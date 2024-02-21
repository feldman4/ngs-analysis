import gzip
import pandas as pd
import os
from .load import load_config, load_fastq, load_samples
from .utils import nglob, csv_frame, assign_format
from .timer import Timer
from .types import Samples
from .sequence import reverse_complement, write_fastq, read_fasta
import re
import shutil


def setup_from_nanopore_fastqs(path_to_fastq_pass, min_reads=10, 
        exclude='unclassified', include=None, write_sample_table=True):
    """Set up analysis from fastq.gz files produced by dorado
    basecalling and barcode demultiplexing.

    Example:
    remote = '/path/to/nanopore/experiment/basecalling/fastq_pass/'
    setup_from_nanopore_fastqs(remote)

    :param path_to_fastq_pass: path to the directory (usually named 
        "fastq_pass") with folders containing .fastq.gz files for each
        barcode
    :param min_reads: skip barcodes with fewer reads
    :param exclude: skip barcodes (folders) matching this regex
    :param write_sample_table: write "samples.csv"
    """
    search = os.path.join(path_to_fastq_pass, '*')
    folders = []
    for x in nglob(search):
        if not os.path.isdir(x) or re.findall(exclude, x):
            continue
        if include is not None and not re.findall(include, x):
            continue
        folders += [x]

    samples = []
    # combine nanopore output fastq files    
    for folder in folders:
        name = os.path.basename(folder)
        files = nglob(f'{folder}/*.fastq.gz')

        num_reads = 0
        arr = []
        for f in files:
            with gzip.open(f, 'rt') as fh:
                try:
                    arr += [fh.read()]
                    num_reads += len(arr[-1].split('\n')) // 4
                except:
                    continue
        if num_reads < min_reads:
            print(f'Skipping {name}, only detected {num_reads} reads')
            continue

        print(f'Wrote 1_reads/{name}.fastq with {num_reads} reads '
              f'from {len(files)} .fastq.gz files')
        with open(f'1_reads/{name}.fastq', 'w') as fh:
            fh.write(''.join(arr))
            
        samples += [name]
    if write_sample_table:
        (pd.DataFrame({'fastq_name': samples, 'sample': samples})
         .pipe(Samples)
         .to_csv('samples.csv', index=None))

    
def demux_reads(flye='flye'):
    shutil.rmtree('4_demux', ignore_errors=True)
    os.makedirs('4_demux')
    print('Cleared files from 4_demux/')

    c = load_config()['demux']
    max_reads_per_assembly = c.get('max_reads_per_assembly', 500)
    
    # filter and sort reads based on this table
    with Timer(verbose='Loading mapped reads...'):
        df_parsed = csv_frame('2_parsed/{sample}.parsed.pq', 
            columns=['read_index', 'read_length', 'insert_length'])
        df_mapped = csv_frame('3_mapped/{sample}.mapped.csv')
        df_fastq = load_fastq(load_samples()['sample'])
    
    key = ['sample', 'read_index']
    df = (df_parsed
     .merge(df_mapped, on=key)
     .merge(df_fastq, on=key)
     .query(c['gate'])
     .pipe(assign_format, fastq_text='{name}\n{read}\n+\n{quality}')
     .groupby(['sample', c['field']])
     .head(max_reads_per_assembly)
    )
    
    num_files = 0
    for (sample, name), df_ in df.groupby(['sample', c['field']]):
        f = f'4_demux/{sample}/{name}.fastq'
        os.makedirs(os.path.dirname(f), exist_ok=True)
        with open(f, 'w') as fh:
            fh.write('\n'.join(df_['fastq_text']))
        num_files += 1
    
    print(f'Wrote {num_files} files to 4_demux/{{sample}}/{{name}}.fastq')

    write_flye_commands(flye=flye)


def write_flye_commands(search='4_demux/*/*.fastq', flags='--nano-hq', flye='flye'):
    """Write commands that can be run with bash, parallel, job submission etc.
    """
    f_out = '4_demux/flye.sh'
    if not shutil.which('flye'):
        print('Warning! flye executable not found')
    arr = []
    for fastq in nglob(search):
        out_dir = fastq.replace(".fastq", "")
        arr += [f'{flye} {flags} {fastq} --out-dir {out_dir}']
    with open(f_out, 'w') as fh:
        fh.write('\n'.join(arr))
    print(f'Wrote {len(arr)} commands to {f_out}')
    print('These can be run in parallel with a command like:')
    print('  cat 4_demux/flye.sh | parallel -j $(nproc)')


def collect_assemblies(backbone_start=None):
    """Collect flye assemblies into a new sample, optionally 
    reverse-complementing/re-indexing.
    
    :param backbone_start: unique sequence at the start of the plasmid, 
        typically ~15 nt
    """
    search = '4_demux/*/*/assembly.fasta'
    files = nglob(search)
    if len(files) == 0:
        print(f'No files matched {search}, aborting')
    
    seqs = [read_fasta(f)[0][1] for f in files]

    df_assembled = (pd.DataFrame({'file': files, 'seq': seqs})
     .assign(sample=lambda x: x['file'].str.split('/').str[1])
     .assign(name=lambda x: x['file'].str.split('/').str[2])
     .assign(length=lambda x: x['seq'].str.len())
    )
    
    if backbone_start == 'pT02':
        backbone_start = 'ATTCTCCTTGGAATT'
    
    if backbone_start is not None:
        reindex = lambda x: reindex_plasmid_sequence(x, [backbone_start])
    else:
        # in one test, flye made an antisense assembly out of sense reads
        reindex = reverse_complement

    df_assembled['read_name'] = df_assembled['sample'] + '_' + df_assembled['name']
    df_assembled['sense'] = df_assembled['seq'].apply(reindex)
    
    f = '1_reads/assembled.fastq'
    write_fastq('1_reads/assembled.fastq', df_assembled['read_name'], df_assembled['sense'])
    print(f'Wrote {len(df_assembled)} assembled sequences to {f}')

    df_samples = load_samples()
    if 'assembled' not in df_samples['sample'].values:
        new_sample = pd.DataFrame({'sample': ['assembled'], 'fastq_name': ['assembled']})
        (pd.concat([df_samples, new_sample])
        .reset_index(drop=True).pipe(Samples)
        .to_csv('samples.csv', index=None)
        )
        print('Added sample "assembled" to samples.csv')
    

def reindex_plasmid_sequence(seq, starting_sequence_candidates, missing='ignore'):
    """Reindex genbank sequence and features.
    """
    seq = seq.upper()
    seq_rc = reverse_complement(seq)
    for start in starting_sequence_candidates:
        start = start.upper()
        if start in seq:
            break
        elif start in seq_rc:
            seq = seq_rc
            break
    else:
        if missing == 'ignore':
            return seq
        else:
            raise ValueError('Plasmid not recognized')
        
    offset = seq.index(start.upper())
    return seq[offset:] + seq[:offset]
