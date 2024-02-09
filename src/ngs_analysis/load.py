import pandas as pd
import pyarrow.parquet as pq
import yaml

from .constants import *
from .sequence import quick_translate, read_fastq
from .types import *
from .utils import *


def load_config():
    with open(config_yaml, 'r') as fh:
        return yaml.safe_load(fh)


def load_reference_dna() -> ReferenceDna:
    """Repeat entries are allowed (same sequence from a different
    source), but `name` must uniquely specify `reference_dna`.
    """
    df = (pd.read_csv(reference_dna_table)
     .pipe(ReferenceDna)
     .assign(reference_dna=lambda x: x['reference_dna'].str.upper())
    )
    (df
     .drop_duplicates(['name', 'reference_dna'])
     .pipe(assert_unique, 'reference_dna')
    )
    return df


def load_sample_plan() -> SamplePlan:
    df = pd.read_csv(sample_plan_table).pipe(SamplePlan)
    if 'source' in df:
        df_ref = load_reference_dna()
        if 'source' not in df_ref:
            raise ValueError('Sample plan has source column but reference does not')
        sample_sources = set(df['source'])
        ref_sources = set(df_ref['source'])
        if missing := sample_sources - ref_sources:
            msg = ('Sample plan refers to sources not in '
            f'reference: {missing}')
            raise ValueError(msg)

    # TODO: check that samples are in samples.csv
    return df


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
            df_designs[field] = translate_non_null(df_designs[field[:-3]])
        elif field not in df_designs:
            raise ValueError(f'Cannot derive {field} from design table')
    return df_designs.pipe(Designs)


def translate_non_null(xs):
    arr = []
    for x in xs:
        if pd.isnull(x):
            arr += [None]
        else:
            arr += [quick_translate(x)]
    return arr


def load_seqs_to_map(sample, simulate, field):
    """Returns dataframe with columns "read_index" and "query". Null values are dropped.
    """
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


def is_field_in_parsed(sample, simulate, field):
    filenames = get_filenames(sample, simulate, field)
    aa = field.endswith('_aa')
    field_no_aa = field[:-3] if aa else field
    return field_no_aa in pq.ParquetFile(filenames['parsed']).schema.names


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


def get_comparison_fields(filter_existing=None):
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

    if filter_existing:
        sample, simulate = filter_existing
        ok = lambda x: is_field_in_parsed(sample, simulate, x)
        direct_match = [x for x in direct_match if ok(x)]
        inferred_match = [(a, b) for a,b in inferred_match if ok(a) and ok(b)]
    return direct_match, inferred_match


def load_data_for_plotting(samples, simulate=False):
    samples = pd.Series(list(samples)).drop_duplicates()

    arr = []
    for sample in samples:
        f = get_filenames(sample, simulate=simulate)['mapped']
        arr += [pd.read_csv(f).assign(sample=sample)]
    df_mapped = pd.concat(arr).reset_index(drop=True)

    direct, indirect = get_comparison_fields()
    fields = direct | set(y for x in indirect for y in x)
    
    return df_mapped, fields


def reference_contains_source():
    return 'source' in pd.read_csv(reference_dna_table, nrows=1)


def load_fastq(samples):
    # load fastq data
    arr = []
    for sample in samples:
        reads, quality, names = read_fastq(
            f'1_reads/{sample}.fastq', include_quality=True, 
            include_name=True, full_name=True)
        (pd.DataFrame({'read': reads, 'name': names, 'quality': quality})
         .assign(read_index=lambda x: np.arange(len(x)))
         .assign(sample=sample)
         .pipe(arr.append)
        )
    return pd.concat(arr)

