"""Analyze sequencing reads and compare DNA and protein parts to a reference.

Pattern-matching is used to define parts (e.g., designs, 
variable linkers, fixed sequences, and barcodes). Identified parts can be 
mapped to inferred reference sequences (e.g., design inferred 
from barcode). See README for details.
"""
import os
import re
import shutil
import subprocess
import tempfile
from glob import glob
from subprocess import DEVNULL
from typing import Tuple

import fire
import Levenshtein
import numpy as np
import pandas as pd

# import quickdna
from .constants import *
from .load import *
from .nanopore import (collect_assemblies, demux_reads,
                       setup_from_nanopore_fastqs, write_flye_commands)
from .sequence import (read_fasta, read_fastq, reverse_complement,
                       write_fake_fastq, write_fasta)
from .timer import Timer
from .types import *
from .utils import nglob


### MAIN STEPS


def setup(clean=False):
    """Setup analysis directory. 
    
    Create working directories; check config, reference_dna.csv, 
    and samples.csv; check NGmerge and mmseqs executables.

    :param clean: whether to remove existing working directories
    """
    
    if clean:
        print('Removing working directories...')
        [shutil.rmtree(d, ignore_errors=True) for d in working_directories]

    print('Creating working directories...')
    create_directories()

    print('Checking config...')
    # TODO: actual type check for the config
    load_config()

    n = len(load_reference_dna())
    print(f'Found {n:,} reference sequences in {reference_dna_table}')

    print(f'Converting reference sequences to design table...')
    dna_to_designs()

    # TODO: check that the input files exist
    n = len(load_samples())
    print(f'Found {n:,} samples in {sample_table}')

    # TODO: check that mmseqs and ngmerge are on path
    if not shutil.which('mmseqs'):
        print('Warning: mmseqs not found')

    if not shutil.which('NGmerge'):
        print('Warning: NGmerge not found')
    

def dna_to_designs():
    """Generate a design table by parsing the reference DNA sequences.
    """
    df_reference = load_reference_dna()
    config = load_config()
    config['left_adapter'] = config.get('left_adapter_reference', config['left_adapter'])
    config['right_adapter'] = config.get('right_adapter_reference', config['right_adapter'])

    parsed = (parse_sequences(config, df_reference['reference_dna'])
     .drop('read_index', axis=1))
    m, n = len(parsed), len(df_reference)
    if m != n:
        raise ValueError(f'failed to parse {n-m}/{n} sequences')
    
    print(f'Parsed {len(parsed):,} sequences from {reference_dna_table}')
    describe_parsed(parsed)

    pd.concat([
        df_reference.drop('reference_dna', axis=1),
        parsed,
    ], axis=1).to_csv(design_table, index=None)
    print(f'Wrote {n:,} rows to {design_table}')


def clear_simulation():
    files = glob('*/simulate/*')
    [os.remove(x) for x in files]


def simulate_single_reads(mutations_per_read: int=0, coverage=None, mutate_to_N=False):
    """Create simulated fastq files based on planned samples.
    """
    df_planned = load_sample_plan()
    df_reference = load_reference_dna()
    df = df_planned.merge(df_reference, how='left')
    if 'coverage' not in df:
        df['coverage'] = coverage
    
    create_directories()

    it = df.groupby(['sample', 'coverage'], sort=False, dropna=False)
    for (sample, coverage), df in it:
        reads = df['reference_dna']

        if not pd.isnull(coverage):
            reads = pd.Series(reads).sample(frac=coverage, replace=True)

        if mutations_per_read != 0:
            reads = mutate_reads(reads, mutations_per_read, mutate_to_N)

        f = f'1_reads/simulate/{sample}.fastq'
        write_fake_fastq(f, reads)
        print(f'Wrote {len(reads):,} single reads to {f}')


def simulate_paired_reads(read_lengths: Tuple[int, int]=(300, 300), 
        mutations_per_read: int=0, coverage=None, mutate_to_N=False):
    """Create simulated fastq files based on planned samples.

    :param read_lengths: forward and reverse read lengths
    """
    df_reference = load_reference_dna()
    df_planned = load_sample_plan().merge(df_reference, how='left')
    if 'coverage' not in df_planned:
        df_planned['coverage'] = coverage

    create_directories()

    it = df_planned.groupby(['sample', 'coverage'], sort=False, dropna=False)
    for (sample, coverage), df in it:
        r1 = df['reference_dna'].str[:read_lengths[0]]
        r2 = df['reference_dna'].apply(reverse_complement).str[:read_lengths[1]]

        if not pd.isnull(coverage):
            r1 = pd.Series(r1).sample(frac=coverage, replace=True)
            r2 = pd.Series(r2).sample(frac=coverage, replace=True)

        if mutations_per_read != 0:
            r1 = mutate_reads(r1, mutations_per_read, mutate_to_N)
            r2 = mutate_reads(r2, mutations_per_read, mutate_to_N)

        write_fake_fastq(f'0_paired_reads/simulate/{sample}_R1.fastq', r1)
        write_fake_fastq(f'0_paired_reads/simulate/{sample}_R2.fastq', r2)
        print(f'Wrote {len(r1):,} paired reads to '
              f'0_paired_reads/simulate/{sample}_R[12].fastq')


def mutate_reads(reads, mutations_per_read, mutate_to_N, seed=0):
    allowed = list('ACGT')
    if mutate_to_N:
        allowed += ['N']
    rs = np.random.RandomState(seed)
    arr = []
    for read in reads:
        mutant = str(read)
        for j in rs.randint(len(read), size=mutations_per_read):
            mutant = mutant[:j - 1] + rs.choice(allowed) + mutant[j:]
        arr += [mutant]
    return arr


def merge_read_pairs(sample, simulate=False):
    """Use NGmerge to merge read pairs.

    :param sample: sample name
    :param simulate: if True, use ./*/simulate/ subdirectories
    """
    filenames = get_filenames(sample, simulate=simulate)
    
    command = [ngmerge, 
    '-1', filenames['paired_r1'], '-2', filenames['paired_r2'], 
    '-o', filenames['reads'], 
    ]
    if filenames['paired_r1'].endswith('.gz'):
        command[-1] += '.gz'
        command += ['-z'] # gzip output

    with Timer(verbose=f'Running ngmerge on sample {sample}...'):
        subprocess.run(command, check=True)
    print(f'Wrote output to {filenames["reads"]}') 


def parse_reads(sample, simulate=False):
    """Find DNA and protein parts within (merged) reads.

    :param sample: sample name
    :param simulate: if True, use ./*/simulate/ subdirectories
    """
    config = load_config()
    filenames = get_filenames(sample, simulate)
    with Timer('', verbose='Loading reads...'):
        reads = read_fastq_or_gz(filenames['reads'])
    with Timer('', verbose=f'Parsing {len(reads):,} reads...'):
        os.makedirs(os.path.dirname(filenames['parsed']), exist_ok=True)
        parse_sequences(config, reads).to_parquet(filenames['parsed'])
    

def map_parsed_reads(sample, simulate=False, mmseqs_max_seqs=10, drop_unmapped=True):
    """For each read, find nearest match for the requested fields. 
    
    Candidate matches from the design table are found via mmseqs search, 
    then the nearest match is determined by Levenshtein distance. For 
    inferred fields, the match is set to the inferred value and the 
    distance to the design table entry accordingly. For barcodes, only
    exact matches are considered.

    :param sample: sample name
    :param simulate: if True, use ./*/simulate/ subdirectories
    :param drop_unmapped: whether to drop reads that fail to map to all
        of the requested fields, if not dropping then the mapped table
        may contain null entries
    """
    filenames = get_filenames(sample, simulate)

    # prepare mmseqs db
    mmseqs_make_design_db()
    
    match_to_ngs, secondary = get_comparison_fields(filter_existing=(sample, simulate))
    if len(match_to_ngs + secondary) == 0:
        print(f'No comparison fields with parsed values available for sample {sample}, skipping')

    arr = []
    for field in match_to_ngs:
        if not is_field_in_parsed(sample, simulate, field):
            print(f'No sequences parsed for field {field}, skipping')
            continue
        with Timer('', verbose=f'Searching for {field} candidates...'):
        # fastmap to get candidates
            df_mapped = mmseqs_prefilter(sample, simulate, field,
                mmseqs_max_seqs=mmseqs_max_seqs)
        # calculate edit distances for candidates
        msg = f'Calculating distances for {len(df_mapped)} {field} candidates...'
        with Timer('', verbose=msg):
            df_match = match_mapped(df_mapped, field)
        arr += [df_match]

    join = 'inner' if drop_unmapped else 'outer'
    df_match = arr[0]
    for df in arr[1:]:
        df_match = df_match.merge(df, on='read_index', how=join)

    if secondary:
        df_designs = load_designs()
        arr = []
        for field_a, field_b in secondary:
            field = f'{field_b}_from_{field_a}'
            if '_aa' in field_b:
                df_designs[field_b] = translate_non_null(df_designs[field_b[:-3]])
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


### UTILITIES


def describe_parsed(df_parsed):
    df_parsed = df_parsed.copy()
    df_parsed[df_parsed == ''] = np.nan
    # pandas applymap was recently deprecated
    if hasattr(df_parsed, 'map'):
        dna_length = df_parsed.select_dtypes('object').map(len, na_action='ignore')
    else:
        dna_length = df_parsed.select_dtypes('object').applymap(len, na_action='ignore')

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


def format_string_to_capture_regex(x, optional=[], **kwargs):
    formatted_string = re.sub(r'{(.*?)}', r'(?P<\1>.*)', x)
    for key, value in kwargs.items():
        tail = '?' if key in optional else ''
        formatted_string = formatted_string.replace(
            f'(?P<{key}>.*)', f'(?P<{key}>{value}){tail}')
    return formatted_string


def parse_sequences(config, sequences):
    """Match DNA sequence based on pattern config.
    First, the sequence between left_adapter and right_adapter is extracted.
    Then, any DNA parts are matched.
    Then, a protein is translated starting at the position matched by the pattern 
    in config.protein.start_at using `translate_to_stop`.
    Finally, the protein is matched using the template string protein.pattern,
    and named matches are added to the result dictionary. Only DNA is stored here.
    """
    re_insert = re.compile('{left_adapter}(.*?){right_adapter}'.format(**config))
    capture = config['protein'].get('capture', {})
    optional = config['protein'].get('optional', [])
    try:
        patterns = config['protein']['patterns']
        re_proteins = [re.compile(format_string_to_capture_regex(
            x, optional, **capture)) for x in patterns]
    except KeyError:
        re_proteins = []

    arr = []
    for i, dna in enumerate(sequences):
        dna = dna.upper()
        entry = {'sense': 1}
        matches = re_insert.findall(dna)
        if len(matches) == 0:
            matches = re_insert.findall(reverse_complement(dna))
            entry['sense'] = -1
        if len(matches) == 0:
            continue
        insert = matches[0]

        entry['read_index'] = i
        entry['read_length'] = len(dna)
        entry['insert_length'] = len(insert)
        for name, pat in config['dna_parts'].items():
            m = re.findall(pat, insert)
            if m:
                entry[name] = m[0]
        if 'protein' not in config:
            continue
        protein_start = re.search(config['protein']['start_at'], insert)
        if protein_start is None:
            continue
        protein_start = protein_start.start()
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
                break
        arr += [entry]
    if len(arr) == 0:
        raise Exception('Nothing found')
    return pd.DataFrame(arr)


def match_mapped(df_mapped: Candidates, field):
    """Filter mapping table for a specific field.

    What happens if inference leads to multiple sequences?
    """
    match = f'{field}_match'
    distance = f'{field}_distance'
    equidistant = f'{field}_match_equidistant'

    # the reference is deduplicated in mmseqs_make_design_db, so 
    # we have to track counts separately
    f = f'3_mapped/map/{field}.counts.csv'
    counts = pd.read_csv(f).rename(columns={'name': 'reference_name'})

    df_mapped = df_mapped.merge(counts)
    df_mapped[distance] = [Levenshtein.distance(a,b) for a,b in df_mapped[['query', 'reference']].values]
    # keep only the minimum values
    min_distance = df_mapped.groupby('read_index')[distance].transform('min')
    df_mapped = df_mapped.loc[lambda x: x[distance] == min_distance]

    A = df_mapped.groupby('read_index')['count'].sum().rename(equidistant).reset_index()
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
    samples = load_samples()['sample']
    if sample not in samples.values:
        raise ValueError(f'Sample {sample} not present in {sample_table}')

    if simulate:
        search = f'1_reads/simulate/{sample}'
    else:
        sample_to_fastq = load_samples().set_index('sample')['fastq_name'].to_dict()
        fastq_search = sample_to_fastq[sample]
        search = f'1_reads/*{fastq_search}*'

    extensions = '.fastq', '.fq', '.fq.gz', '.fastq.gz'
    def search_extensions(search):
        return sum([glob(search + x) for x in extensions], [])

    r = search_extensions(search)
    search_paired = search.replace('1_reads', '0_paired_reads')
    r1 = search_extensions(search_paired + '_R1')
    r2 = search_extensions(search_paired + '_R2')
    paired_reads_found = len(r1) == 1 & len(r2) == 1

    if len(r1) > 1 or len(r2) > 1:
        raise ValueError(f'Too many paired read files for sample {sample}, found: {r1 + r2}')
    
    if len(r) > 1:
        raise ValueError(f'Too many read files for sample {sample}, found: {r}')

    if not (paired_reads_found or len(r) == 1):
        msg = [
            f'Neither paired reads nor single reads provided for '
            f'sample {sample}, searches returned:',
            f'{search_paired}: {nglob(search_paired)}',
            f'{search}: {nglob(search)}',
        ]
        raise ValueError('\n'.join(msg))

    add_simulate = 'simulate/' if simulate else ''
    filenames = {
        'reads': f'1_reads/{add_simulate}{sample}.fastq',
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

def read_fastq_or_gz(filename):
    if os.path.exists(filename):
        return read_fastq(filename)
    else:
        return read_fastq(filename + '.gz')


def mmseqs_make_design_db():
    df_designs = load_designs()
    fields, _ = get_comparison_fields()
    with Timer('', verbose=f'Creating mmseqs index for fields: {", ".join(fields)}'):
        for field in fields:
            aa = field.endswith('_aa')
            f = f'3_mapped/map/{field}.fa' # should live in get_filenames
            seqs_to_index = (df_designs[['name', field]]
             .assign(count=lambda x: x.groupby(field).transform('size'))
             .drop_duplicates(field)
            )
            write_fasta(f, seqs_to_index[['name', field]])
            make_mmseqs_db(f, aa, is_query=False)

            f = f'3_mapped/map/{field}.counts.csv'
            seqs_to_index[['name', 'count']].to_csv(f, index=None)


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
    seqs_to_map = load_seqs_to_map(sample, simulate, field)
    if aa:
        kmer = 6
    else:
        # what should this be?
        kmer = 6

    write_fasta(filenames['map query'], seqs_to_map)
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
        '--spaced-kmer-mode', '0',
        # '--exact-kmer-matching', '1', # prefilter only retains exact kmer matches
        )
        if verbose:
            print(' '.join(command))
        result = subprocess.run(command, capture_output=True)
        # result = subprocess.run(command)
    if result.returncode != 0:
    # if result.returncode != 0 or True:
        msg = '\n'.join([str(x) for x in 
            (result.args,
            result.stdout.decode(),
            result.stderr.decode())
            ])
        raise ValueError(f'Non-zero return code in mmseqs prefilter:\n{msg}')
        
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
    seqs_to_map = load_seqs_to_map(sample, simulate, field)
    kmer = MMSEQS_KMER_AA if aa else MMSEQS_KMER_DNA
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


def create_directories():
    for d in working_directories:
        os.makedirs(f'{d}/simulate', exist_ok=True)


def main():
    # order is preserved
    commands = [
        'setup',
        'dna_to_designs', 
        'simulate_paired_reads', 
        'merge_read_pairs', 
        'parse_reads',
        'map_parsed_reads',
        'setup_from_nanopore_fastqs',
        'demux_reads',
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
    

if __name__ == '__main__':
    main()
