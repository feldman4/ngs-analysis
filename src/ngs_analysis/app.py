"""Analyze sequencing reads and compare DNA and protein parts to a reference.

Pattern-matching is used to define parts (e.g., designs, 
variable linkers, fixed sequences, and barcodes). Identified parts can be 
mapped to inferred reference sequences (e.g., design inferred 
from barcode). See README for details.
"""
import os
import re
import subprocess
import tempfile
from glob import glob
from subprocess import DEVNULL
from typing import Optional, Tuple

import fire
import Levenshtein
import numpy as np
import pandas as pd
import pandera as pa
import yaml
from pandera.typing import Series

# import quickdna
from .sequence import (read_fasta, read_fastq, reverse_complement,
                       translate_dna, write_fake_fastq, write_fasta)
from .timer import Timer

# user-provided
config_yaml = 'config.yaml' # how to parse and compare
sample_table = 'samples.csv' # NGS samples
possible_dna_table = 'possible_dna.csv' # designed insert DNA
planned_dna_table = 'planned_dna.csv' # cloning plan

working_directories = '0_paired_reads', '1_reads', '2_parsed', '3_mapped', '3_mapped/map'

# generated or provided
design_table = 'designs.csv' # designed DNA or protein parts

ngmerge = 'NGmerge'
bwa = 'bwa'
mmseqs = 'mmseqs'

MMSEQS_KMER_AA = 6
MMSEQS_KMER_DNA = 12

### TYPES
# TODO: joint uniqueness
class PossibleDNA(pa.DataFrameModel):
    subpool: Optional[Series[str]] = pa.Field(nullable=False, coerce=True)
    name: Series[str] = pa.Field(nullable=False, unique=True)
    vector_dna: Series[str] = pa.Field(nullable=False, unique=True)


class PlannedDNA(pa.DataFrameModel):
    sample: Series[str] = pa.Field(nullable=False, coerce=True)
    vector_dna: Series[str] = pa.Field(nullable=False)


class Samples(pa.DataFrameModel):
    sample: Series[str] = pa.Field(nullable=False, unique=True, 
                                   coerce=True)
    fastq_name: Series[str] = pa.Field(nullable=False, unique=True, 
                                       coerce=True)


# build this by using config.yaml to parse PossibleDNA
class Designs(pa.DataFrameModel):
    subpool: Optional[Series[str]] = pa.Field(nullable=False, coerce=True)
    name: Series[str] = pa.Field(nullable=False, coerce=True)


class Candidates(pa.DataFrameModel):
    read_index: Series[int] = pa.Field(nullable=False)
    query: Series[str] = pa.Field(nullable=False)
    reference: Series[str] = pa.Field(nullable=False)
    reference_name: Series[str] = pa.Field(nullable=False)


### MAIN STEPS
def setup():
    print('Creating directories...')
    [os.makedirs(f'{d}/simulate', exist_ok=True) for d in working_directories]

    print('Checking config...')
    # TODO: actual type check for the config
    load_config()

    n = len(load_possible_dna())
    print(f'Found {n} reference sequences in {possible_dna_table}')

    # TODO: check that the input files exist
    n = len(load_samples())
    print(f'Found {n} samples in {sample_table}')
    
    

# ok
def dna_to_designs():
    """Generate a design table by parsing the possible DNA sequences.
    """
    df_possible = load_possible_dna()
    config = load_config()
    parsed = (parse_sequences(config, df_possible['vector_dna'])
     .drop('read_index', axis=1))
    m, n = len(parsed), len(df_possible)
    if m != n:
        raise ValueError(f'failed to parse {n-m}/{n} sequences')
    
    print(f'Parsed {len(parsed):,} sequences from {possible_dna_table}')
    describe_parsed(parsed)

    pd.concat([
        df_possible.drop('vector_dna', axis=1),
        parsed,
    ], axis=1).to_csv(design_table, index=None)
    print(f'Wrote {n:,} rows to {design_table}')

# ok
def simulate_fastqs(read_lengths: Tuple[int, int]=(300, 300)):
    """Create simulated fastq files based on planned samples.

    :param read_lengths: forward and reverse read lengths
    """
    df_planned = load_planned_dna()
    os.makedirs('1_reads/simulate', exist_ok=True)

    for sample, df in df_planned.groupby('sample', sort=False):
        r1 = df['vector_dna'].str[:read_lengths[0]]
        r2 = df['vector_dna'].apply(reverse_complement).str[:read_lengths[1]]

        write_fake_fastq(f'1_reads/simulate/{sample}_R1.fastq', r1)
        write_fake_fastq(f'1_reads/simulate/{sample}_R2.fastq', r2)
        print(f'Wrote {len(r1):,} paired reads to '
              f'1_reads/simulate/{sample}_R[12].fastq')

# ok
def merge_read_pairs(sample, simulate=False):
    """Use NGmerge to merge read pairs.

    :param sample: sample name
    :param simulate: if True, use ./*/simulate/ subdirectories
    """
    filenames = get_filenames(sample, simulate=simulate)
    
    command = (ngmerge, 
    '-1', filenames['paired_r1'], '-2', filenames['paired_r2'], 
    '-o', filenames['reads'], 
    '-z', # always gzip output
    )
    with Timer(verbose=f'Running ngmerge on sample {sample}...'):
        subprocess.run(command, check=True)
    print(f'Wrote output to {filenames["reads"]}') 

# ok
def parse_reads(sample, simulate=False):
    config = load_config()
    filenames = get_filenames(sample, simulate)
    with Timer('', verbose='Loading reads...'):
        reads = read_fastq(filenames['reads'])
    with Timer('', verbose=f'Parsing {len(reads):,} reads...'):
        os.makedirs(os.path.dirname(filenames['parsed']), exist_ok=True)
        parse_sequences(config, reads).to_parquet(filenames['parsed'])
    
# ok
def map_parsed_reads(sample, simulate=False, mmseqs_max_seqs=10):
    """For each read, find nearest match for the requested fields. 
    
    Candidate matches from the design table are found via mmseqs search, 
    then the nearest match is determined by Levenshtein distance. For 
    inferred fields, the match is set to the inferred value and the 
    distance to the design table entry accordingly. For barcodes, only
    exact matches are considered.

    :param sample: sample name
    :param simulate: if True, use ./*/simulate/ subdirectories
    """
    filenames = get_filenames(sample, simulate)

    # prepare mmseqs db
    mmseqs_make_design_db()
    
    match_to_ngs, secondary = get_comparison_fields()
    arr = []
    for field in match_to_ngs:
        if field == 'barcode': 
            # for barcode, only permit exact matches
            with Timer('', verbose=f'Matching {field}...'):
                arr += [map_barcode(sample, simulate, field)]
        else:
            with Timer('', verbose=f'Searching for {field} candidates...'):
            # fastmap to get candidates
                df_mapped = mmseqs_prefilter(sample, simulate, field,
                 mmseqs_max_seqs=mmseqs_max_seqs)
            # calculate edit distances for candidates
            msg = f'Calculating distances for {len(df_mapped)} {field} candidates...'
            with Timer('', verbose=msg):
                df_match = match_mapped(df_mapped, field)
            arr += [df_match]

    df_match = arr[0]
    for df in arr[1:]:
        df_match = df_match.merge(df, on='read_index')

    if secondary:
        df_designs = load_designs()
        arr = []
        for field_a, field_b in secondary:
            field = f'{field_b}_from_{field_a}'
            if '_aa' in field_b:
                df_designs[field_b] = [quick_translate(x) for x in df_designs[field_b[:-3]]]
            name_to_field = df_designs.set_index('name')[field_b].to_dict()
            # define reference for field_b using field_a
            reference = (df_match
            .rename(columns={f'{field_a}_match': 'reference_name'})
            .assign(reference=lambda x: x['reference_name'].map(name_to_field))
            .dropna()
            )
            # calculate distance between the parsed value for field_b and
            # the reference (where both exist)
            df = (load_seqs_to_map(sample, simulate, field_b)
                .merge(reference, on='read_index')
                .pipe(Candidates)
                .pipe(match_mapped, field)
                .drop(f'{field}_match_equidistant', axis=1)
            )
            df_match = df_match.merge(df, on='read_index')
    
    df_match.to_csv(filenames['mapped'], index=None)
    print(f'Wrote match table ({len(df_match):,} rows) to {filenames["mapped"]}')


### LOAD


def load_config():
    with open(config_yaml, 'r') as fh:
        return yaml.safe_load(fh)


def load_possible_dna() -> PossibleDNA:
    return pd.read_csv(possible_dna_table).pipe(PossibleDNA)


def load_planned_dna() -> PlannedDNA:
    return pd.read_csv(planned_dna_table).pipe(PlannedDNA)


def load_samples() -> Samples:
    return pd.read_csv(sample_table).pipe(Samples)


def load_designs() -> Designs:
    """Add translations that are needed for comparison but not already present.
    """
    df_designs = pd.read_csv(design_table)
    fields, _ = get_comparison_fields()
    for field in fields:
        if field.endswith('_aa') and field[:-3] in df_designs:
            if field in df_designs:
                raise ValueError(f'Cannot provide both DNA and amino acid for {field}')
            df_designs[field] = [quick_translate(x) for x in df_designs[field[:-3]]]
        elif field not in df_designs:
            raise ValueError(f'Cannot derive {field} from design table')
    return df_designs.pipe(Designs)


def load_seqs_to_map(sample, simulate, field):
    filenames = get_filenames(sample, simulate, field)

    aa = field.endswith('_aa')
    field_no_aa = field[:-3] if aa else field

    seqs_to_map = (pd.read_parquet(filenames['parsed'], columns=['read_index', field_no_aa])
     .rename(columns={field_no_aa: 'query'})
     .dropna()
    )
    if aa:
        seqs_to_map['query'] = [quick_translate(x) for x in seqs_to_map['query']]

    return seqs_to_map


### UTILITIES


def describe_parsed(df_parsed):
    df_parsed = df_parsed.copy()
    df_parsed[df_parsed == ''] = np.nan
    dna_length = df_parsed.select_dtypes('object').map(len, na_action='ignore')
    aa_length = dna_length / 3

    null_summary = pd.concat([
        (df_parsed.isnull() * 1).sum().rename('null entries'),
        (df_parsed.isnull() * 1).mean().rename('null fraction'),
    ], axis=1)
    stats = ['count', 'min', 'median', 'max']
    describe = lambda x: x.describe().rename({'50%': 'median'}).loc[stats].astype(int)

    to_string = lambda x: '    ' + x.to_string().replace('\n', '\n    ')
    print('Null entries:')
    print(null_summary.pipe(to_string))
    print('DNA length:')
    print(dna_length.pipe(describe).pipe(to_string))
    print('Protein length:')
    print(aa_length.pipe(describe).pipe(to_string))


def quick_translate(seq):
    # return str(quickdna.DnaSequence(seq).translate())
    n = len(seq) - len(seq) % 3
    return translate_dna(seq[:n])


def format_string_to_capture_regex(x, optional=[], **kwargs):
    formatted_string = re.sub(r'{(.*?)}', r'(?P<\1>.*)', x)
    for key, value in kwargs.items():
        tail = '?' if key in optional else ''
        formatted_string = formatted_string.replace(
            f'(?P<{key}>.*)', f'(?P<{key}>{value}){tail}')
    return formatted_string


def get_comparison_fields():
    """Parse fields for direct matching and field pairs for inferred matching.
    """
    direct_match, inferred_match = set(), set()
    config = load_config()
    for c in config['compare']:
        if '->' in c:
            a, b = [x.strip() for x in c.split('->')]
            direct_match |= {a}
            inferred_match |= {(a, b)}
        else:
            direct_match |= {c.strip()}

    return direct_match, inferred_match


def parse_sequences(config, sequences):
    """Match DNA sequence based on pattern config.
    First, the sequence between left_adapter and right_adapter is extracted.
    Then, any DNA parts are matched.
    Then, a protein is translated starting at the position matched by the pattern 
    in config.protein.start_at using `translate_to_stop`.
    Finally, the protein is matched using the template string protein.pattern,
    and named matches are added to the result dictionary. Only DNA is stored here.
    """
    re_insert = re.compile('{left_adapter}(.*){right_adapter}'.format(**config))
    import os
    os.wtf = re_insert
    capture = config['protein'].get('capture', {})
    optional = config['protein'].get('optional', [])
    re_proteins = [re.compile(format_string_to_capture_regex(
        x, optional, **capture)) for x in config['protein']['patterns']]
    
    arr = []
    for i, dna in enumerate(sequences):
        dna = dna.upper()
        entry = {}
        try:
            insert = re_insert.findall(dna)[0]
        except IndexError:
            continue
        entry['read_index'] = i
        entry['read_length'] = len(dna)
        entry['insert_length'] = len(insert)
        for name, pat in config['dna_parts'].items():
            m = re.findall(pat, insert)
            if m:
                entry[name] = m[0]
        
        protein_start = re.search(config['protein']['start_at'], insert).start()
        protein = quick_translate(insert[protein_start:])
        protein = protein.split('*')[0]
        entry['cds'] = insert[protein_start:][:len(protein) * 3]
        
        for re_protein in re_proteins:
            match = re_protein.match(protein)
            if match:
                for i, name in enumerate(match.groupdict()):
                    name = name.rstrip('_') # collapse redundant names
                    start, end = match.span(i + 1)
                    entry[name] = entry['cds'][start*3:end*3]
                arr += [entry]
                # assert False
                break
    if len(arr) == 0:
        raise Exception('Nothing found')
    return pd.DataFrame(arr)


def map_barcode(sample, simulate, field):
    """Only exact matches, distance and equidistant not returned.
    """
    filenames = get_filenames(sample, simulate, field)

    aa = field.endswith('_aa')
    field_no_aa = field[:-3] if aa else field
    seqs_to_map = pd.read_parquet(filenames['parsed'], columns=['read_index', field_no_aa])
    if aa:
        seqs_to_map[field] = [quick_translate(x) for x in seqs_to_map[field_no_aa]]
        seqs_to_map = seqs_to_map[['read_index', field]]    
    ref_seqs = (read_fasta(filenames['map mmseqs target fa'], as_df=True)
    .rename(columns={'name': 'match', 'seq': 'reference'})
    )

    match = f'{field}_match'

    return (seqs_to_map
     .merge(ref_seqs, left_on=field, right_on='reference', how='left')
     .rename(columns={'match': match})
     [['read_index', match]]
    )


def match_mapped(df_mapped: Candidates, field):
    """Filter mapping table for a specific field.

    What happens if inference leads to multiple sequences?
    """
    match = f'{field}_match'
    distance = f'{field}_distance'
    equidistant = f'{field}_match_equidistant'

    df_mapped = df_mapped.copy()
    df_mapped[distance] = [Levenshtein.distance(a,b) for a,b in df_mapped[['query', 'reference']].values]
    # keep only the minimum values
    min_distance = df_mapped.groupby('read_index')[distance].transform('min')
    df_mapped = df_mapped.loc[lambda x: x[distance] == min_distance]

    A = df_mapped.groupby('read_index').size().rename(equidistant).reset_index()
    B = (df_mapped.sort_values(distance).drop_duplicates('read_index')
     .rename(columns={'reference_name': match})[['read_index', match, distance]])
    return A.merge(B)[['read_index', match, distance, equidistant]]


### FILES


def get_filenames(sample, simulate=False, field=None):
    """Get filenames used throughout pipeline for this sample.

    Fastq files ending with .fastq, .fq., .fastq.gz, or .fq.gz are identified 
    based on the sample table. In simulate mode, the fastq name is identical 
    to the sample name.
    """
    if simulate:
        search = f'1_reads/simulate/{sample}_R1'
    else:
        sample_to_fastq = load_samples().set_index('sample')['fastq_name'].to_dict()
        fastq_search = sample_to_fastq[sample]
        search = f'1_reads/*{fastq_search}*'

    extensions = '.fastq', '.fq', '.fq.gz', '.fastq.gz'
    def search_extensions(search):
        return sum([glob(search + x) for x in extensions], [])

    r = search_extensions(search)
    r1 = search_extensions('paired_' + search + '_R1')
    r2 = search_extensions('paired_' + search + '_R2')
    paired_reads_found = len(r1) == 1 & len(r2) == 1

    if len(r1) > 1 or len(r2) > 1:
        raise ValueError(f'Too many paired read files for sample {sample}, found: {r1 + r2}')
    
    if len(r) > 1:
        raise ValueError(f'Too many read files for sample {sample}, found: {r}')

    if not (paired_reads_found or len(r) == 1):
        raise ValueError(f'Neither reads nor paired reads provided for sample {sample}, searched with: {search}')

    add_simulate = 'simulate/' if simulate else ''
    filenames = {
        'reads': f'1_reads/{add_simulate}{sample}.fastq.gz',
        'parsed': f'2_parsed/{add_simulate}{sample}.parsed.pq',
        'mapped': f'3_mapped/{add_simulate}{sample}.mapped.csv',
    }
    if paired_reads_found:
        filenames |= {'paired_r1': r1[0], 'paired_r2': r2[0]}

    if field:
        filenames['map query'] = f'3_mapped/map/{add_simulate}{sample}_{field}.fa'
        filenames['map mmseqs target db'] = f'3_mapped/map/{field}.target.db'
        filenames['map mmseqs target fa'] = f'3_mapped/map/{field}.fa'
        filenames['map mmseqs query db'] = f'3_mapped/map/{add_simulate}{sample}_{field}.query.db'
        filenames['map mmseqs result db'] = f'3_mapped/map/{add_simulate}{sample}_{field}.result.db'
        filenames['map mmseqs prefilter db'] = f'3_mapped/map/{add_simulate}{sample}_{field}.prefilter.db'
        filenames['map mmseqs prefilter tsv'] = f'3_mapped/map/{add_simulate}{sample}_{field}.prefilter.tsv'
        filenames['map mmseqs result m8'] = f'3_mapped/map/{add_simulate}{sample}_{field}.result.m8'
        filenames['map result'] = f'3_mapped/map/{add_simulate}{sample}_{field}.csv'

    return filenames



### MMSEQS


def mmseqs_make_design_db():
    df_designs = load_designs()
    fields, _ = get_comparison_fields()
    with Timer('', verbose=f'Creating mmseqs index for fields: {", ".join(fields)}'):
        for field in fields:
            aa = field.endswith('_aa')
            f = f'3_mapped/map/{field}.fa' # should live in get_filenames
            seqs_to_index = df_designs[['name', field]]
            write_fasta(f, seqs_to_index.drop_duplicates(field))
            if field != 'barcode':
                make_mmseqs_db(f, aa, is_query=False)


def make_mmseqs_db(fasta_file, dbtype_is_amino_acid, is_query):
    """
    """
    db = fasta_file.replace('.fa', '.query.db' if is_query else '.target.db')
    command = (mmseqs, 'createdb', fasta_file, db,
     '--dbtype', '1' if dbtype_is_amino_acid else '2'
    )
    result = subprocess.run(command, stdout=DEVNULL, stderr=DEVNULL)
    if result.returncode != 0:
        raise ValueError('Non-zero return code in mmseqs createdb')


def mmseqs_prefilter(sample, simulate, field, mmseqs_max_seqs=10, verbose=False):
    filenames = get_filenames(sample, simulate, field)
    
    # make query database
    aa = field.endswith('_aa')
    field_no_aa = field[:-3] if aa else field
    seqs_to_map = (pd.read_parquet(filenames['parsed'], columns=['read_index', field_no_aa])
     .dropna() # remove reads that failed parsing
     .rename(columns={field_no_aa: 'query'})
    )
    if aa:
        seqs_to_map['query'] = [quick_translate(x) for x in seqs_to_map['query']]
        kmer = 6
    else:
        kmer = 12

    write_fasta(filenames['map query'], seqs_to_map.dropna())
    make_mmseqs_db(filenames['map query'], aa, is_query=True)

    # search
    [os.remove(x) for x in glob(filenames['map mmseqs prefilter db'] + '*')]
    with tempfile.TemporaryDirectory() as tmp:
        command = ('mmseqs', 'prefilter', 
        filenames['map mmseqs query db'], # query
        filenames['map mmseqs target db'], # target
        filenames['map mmseqs prefilter db'], 
        '-k', str(kmer), 
        '--max-seqs', str(mmseqs_max_seqs),
        # without these flags, mmseqs will sometimes throw out exact matches
        '--comp-bias-corr', '0',
        '--mask', '0',
        # '--exact-kmer-matching', '1', # prefilter only retains exact kmer matches
        )
        if verbose:
            print(' '.join(command))
        result = subprocess.run(command, stderr=DEVNULL, stdout=DEVNULL)
        # result = subprocess.run(command)
    if result.returncode != 0:
        raise ValueError('Non-zero return code in mmseqs prefilter')

    # convert to tsv
    command = ('mmseqs', 'createtsv', 
    filenames['map mmseqs query db'], # query
    filenames['map mmseqs target db'], # target
    filenames['map mmseqs prefilter db'],
    filenames['map mmseqs prefilter tsv'],
    )
    if verbose:
        print(' '.join(command))
    result = subprocess.run(command, stderr=DEVNULL, stdout=DEVNULL)
    if result.returncode != 0:
        raise ValueError('Non-zero return code in mmseqs prefilter')
    
    ref_seqs = (read_fasta(filenames['map mmseqs target fa'], as_df=True)
     .rename(columns={'name': 'reference_name', 'seq': 'reference'}))

    df_mmseqs = (read_mmseqs_prefilter_tsv(filenames['map mmseqs prefilter tsv'])
    .rename(columns={'query': 'read_index', 'target': 'reference_name'})
    [['read_index', 'reference_name']]
    )
    return df_mmseqs.merge(seqs_to_map).merge(ref_seqs).pipe(Candidates)
    

def read_mmseqs_prefilter_tsv(filename):
    return pd.read_csv(filename, sep='\t', names=('query', 'target', 'score', 'diagonal'))


def read_m8(filename):
    columns = ('query', 'target', 'sequence_identity', 
               'alignment_length', 'num_mismatches', 'num_gaps', 
               'query_start', 'query_end', 'target_start', 
               'target_end', 'e_value', 'bit_score')
    return pd.read_csv(filename, sep=r'\s+', header=None, names=columns)


def mmseqs_search(sample, field, simulate=False, min_seq_id=0.8) -> Candidates:
    """Not in use.
    
    Use mmseqs search to map `field` in parsed reads to the design table.

    :param sample: sample name
    :param field: a field name, if it ends in "_aa" the translation will be mapped
    :param simulate: if True, use ./*/simulate/ subdirectories
    :param min_seq_id: mmseqs2 prefilter parameter
    """
    filenames = get_filenames(sample, simulate, field)
    
    # make query database
    aa = field.endswith('_aa')
    field_no_aa = field[:-3] if aa else field
    seqs_to_map = (pd.read_parquet(filenames['parsed'], columns=['read_index', field_no_aa])
     .dropna() # remove reads that failed parsing
     .rename(columns={field_no_aa: 'query'})
    )
    if aa:
        seqs_to_map['query'] = [quick_translate(x) for x in seqs_to_map['query']]
        kmer = MMSEQS_KMER_AA
    else:
        kmer = MMSEQS_KMER_DNA
    write_fasta(filenames['map query'], seqs_to_map.dropna())
    make_mmseqs_db(filenames['map query'], aa, is_query=True)

    # search
    [os.remove(x) for x in glob(filenames['map mmseqs result db'] + '*')]
    with tempfile.TemporaryDirectory() as tmp:
        command = ('mmseqs', 'search', 
        filenames['map mmseqs query db'], # query
        filenames['map mmseqs target db'], # target
        filenames['map mmseqs result db'], 
        tmp,
        '--min-seq-id', str(min_seq_id),
        '-k', str(kmer), # memory usage is 21**k * 8 bytes (protein) or 5**k * 8 bytes (DNA)
        '-s', '1', # sensitivity, 1 is fastest, 7.5 is most sensitive
        '-c', '0.9', # minimum alignment coverage on query/target
        '--exact-kmer-matching', '1', # prefilter only retains exact kmer matches
        # split query database based on available memory
        # '--split 0 --split-mode 1',
        )
        if not aa:
            command += ('--search-type', '3')
        print(' '.join(command))
        # result = subprocess.run(command, stderr=DEVNULL, stdout=DEVNULL)
        result = subprocess.run(command)
    if result.returncode != 0:
        raise ValueError('Non-zero return code in mmseqs search')

    # convert results
    command = ('mmseqs', 'convertalis', 
    filenames['map mmseqs query db'], # query
    filenames['map mmseqs target db'], # target
    filenames['map mmseqs result db'], # result
    filenames['map mmseqs result m8'], # result
    )
    # print(' '.join(command))
    result = subprocess.run(command, stderr=DEVNULL, stdout=DEVNULL)
    if result.returncode != 0:
        raise ValueError('Non-zero return code in mmseqs convertalis')
    
    ref_seqs = read_fasta(filenames['map mmseqs target fa'], as_df=True).rename(columns={'name': 'reference_name', 'seq': 'reference'})

    df_mmseqs = (read_m8(filenames['map mmseqs result m8'])
    .rename(columns={'query': 'read_index', 'target': 'reference_name'})
    [['read_index', 'reference_name']].drop_duplicates()
    )
    return df_mmseqs.merge(seqs_to_map).merge(ref_seqs).pipe(Candidates)


def load_fastmap(filename):
    """A bit slow, 10% of fastmap runtime
    """
    txt = open(filename).read()
    queries = txt.strip().split('//')
    arr = []
    for q in queries[:-1]:
        if not q:
            continue
        lines = q.strip().split('\n')

        name = lines[0].split()[1]
        for x in set([l.split()[-1].split(':')[0] for l in lines[1:]]):
            arr += [{'query': name, 'match': x}]
    return pd.DataFrame(arr)


if __name__ == '__main__':

    # order is preserved
    commands = [
        'setup',
        'dna_to_designs', 
        'simulate_fastqs', 
        'merge_read_pairs', 
        'parse_reads',
        'map_parsed_reads',
    ]
    # if the command name is different from the function name
    named = {
        # 'search': search_app,
        }

    final = {}
    for k in commands:
        try:
            final[k] = named[k]
        except KeyError:
            final[k] = eval(k)

    try:
        fire.Fire(final)
    except BrokenPipeError:
        pass
    

