import gzip
from glob import glob

import numpy as np
import pandas as pd
from natsort import natsorted

from .constants import CODONS


watson_crick = {'A': 'T',
                'T': 'A',
                'C': 'G',
                'G': 'C',
                'U': 'A',
                'N': 'N'}

watson_crick.update({k.lower(): v.lower()
                     for k, v in watson_crick.items()})

iupac = {'A': ['A'],
 'C': ['C'],
 'G': ['G'],
 'T': ['T'],
 'M': ['A', 'C'],
 'R': ['A', 'G'],
 'W': ['A', 'T'],
 'S': ['C', 'G'],
 'Y': ['C', 'T'],
 'K': ['G', 'T'],
 'V': ['A', 'C', 'G'],
 'H': ['A', 'C', 'T'],
 'D': ['A', 'G', 'T'],
 'B': ['C', 'G', 'T'],
 'N': ['G', 'A', 'T', 'C']}

codon_maps = {}


def read_fasta(f, as_df=False):
    if f.endswith('.gz'):
        fh = gzip.open(f)
        txt = fh.read().decode()
    else:
        fh = open(f, 'r')
        txt = fh.read()
    fh.close()
    records = parse_fasta(txt)
    if as_df:
        return pd.DataFrame(records, columns=('name', 'seq'))
    else:
        return records


def parse_fasta(txt):
    entries = []
    txt = '\n' + txt.strip()
    for raw in txt.split('\n>'):
        name = raw.split('\n')[0].strip()
        seq = ''.join(raw.split('\n')[1:]).replace(' ', '')
        if name:
            entries += [(name, seq)]
    return entries


def write_fasta(filename, list_or_records):
    if isinstance(list_or_records, pd.DataFrame) and list_or_records.shape[1] == 2:
        list_or_records = list_or_records.values
    list_or_records = list(list_or_records)
    with open(filename, 'w') as fh:
        fh.write(format_fasta(list_or_records))


def write_fake_fastq(filename, list_or_records):
    if filename.endswith('.gz'):
        fh = gzip.open(filename, 'wt')
    else:
        fh = open(filename, 'w')
    fh.write(format_fake_fastq(list_or_records))
    fh.close()


def format_fake_fastq(list_or_records):
    """Generates a fake header for each read that is sufficient to fool bwa/NGmerge.
    """
    fake_header = '@M08044:78:000000000-L568G:1:{tile}:{x}:{y} 1:N:0:AAAAAAAA'
    if isinstance(next(iter(list_or_records)), str):
        records = list_to_records(list_or_records)
    else:
        records = list_or_records

    max_value = 1000
    lines = []
    for i, (_, seq) in enumerate(records):
        tile, rem = divmod(i, max_value**2)
        x, y = divmod(rem, max_value)
        lines.extend([fake_header.format(tile=tile, x=x, y=y), seq.upper(), '+', 'G' * len(seq)])
    return '\n'.join(lines)


def write_fastq(filename, names, sequences, quality_scores):
    with open(filename, 'w') as fh:
        fh.write(format_fastq(names, sequences, quality_scores))


def format_fastq(names, sequences, quality_scores):
    lines = []
    for name, seq, q_score in zip(names, sequences, quality_scores):
        lines.extend([name, seq, '+', q_score])
    return '\n'.join(lines)


def list_to_records(xs):
    n = len(xs)
    width = int(np.ceil(np.log10(n)))
    fmt = '{' + f':0{width}d' + '}'
    records = []
    for i, s in enumerate(xs):
        records += [(fmt.format(i), s)]
    return records


def format_fasta(list_or_records):
    if len(list_or_records) == 0:
        records = []
    elif isinstance(list_or_records[0], str):
        records = list_to_records(list_or_records)
    else:
        records = list_or_records
    
    lines = []
    for name, seq in records:
        lines.extend([f'>{name}', str(seq)])
    return '\n'.join(lines)


def fasta_frame(files_or_search):
    """Convenience function, pass either a list of files or a 
    glob wildcard search term.
    """
    
    if isinstance(files_or_search, str):
        files = natsorted(glob(files_or_search))
    else:
        files = files_or_search

    cols = ['name', 'seq', 'file_ix', 'file']
    records = []
    for f in files:
        for i, (name, seq) in enumerate(read_fasta(f)):
            records += [{
                'name': name, 'seq': seq, 'file_ix': i, 
                'file': f,
            }]

    return pd.DataFrame(records)[cols]


def cast_cols(df, int_cols=tuple(), float_cols=tuple(), str_cols=tuple(), 
              cat_cols=tuple(), uint16_cols=tuple()):
    return (df
           .assign(**{c: df[c].astype(int) for c in int_cols})
           .assign(**{c: df[c].astype(np.uint16) for c in uint16_cols})
           .assign(**{c: df[c].astype(float) for c in float_cols})
           .assign(**{c: df[c].astype(str) for c in str_cols})
           .assign(**{c: df[c].astype('category') for c in cat_cols})
           )


def translate_dna(s):
    assert len(s) % 3 == 0, 'length must be a multiple of 3'
    return ''.join([CODONS.get(s[i*3:(i+1)*3], 'X') for i in range(int(len(s)/3))])


def reverse_complement(seq):
    return ''.join(watson_crick[x] for x in seq)[::-1]


def get_kmers(s, k):
    n = len(s)
    return [s[i:i+k] for i in range(n-k+1)]


def read_fastq(filename, max_reads=1e12, include_quality=False, include_name=False, 
               include_index=False, progress=lambda x: x, full_name=False):
    if max_reads is None:
        max_reads = 1e12
    if filename.endswith('gz'):
        fh = gzip.open(filename, 'rt')
    else:
        fh = open(filename, 'r')
    reads, quality_scores, names, indices = [], [], [], []
    read_count = 0
    for i, line in progress(enumerate(fh)):
        if i % 4 == 1:
            reads.append(line.strip())
            read_count += 1
        if include_quality and i % 4 == 3:
            quality_scores.append(line.strip())
        if include_name and i % 4 == 0:
            if full_name:
                names.append(line.strip())
            else:
                names.append(':'.join(line.split()[0].split(':')[3:7]))
        if include_index and i % 4 == 0:
            indices.append(line.split(':')[-1].strip())
        if i % 4 == 3 and read_count >= max_reads:
            break
        
    fh.close()
    if include_quality or include_name or include_index:
        return_val = (reads,)
        if include_quality:
            return_val += (quality_scores,)
        if include_name:
            return_val += (names,)
        if include_index:
            return_val += (indices,)
        return return_val
    else:
        return reads


def try_translate_dna(s):
    try:
        return translate_dna(s)
    except:
        return None


def translate_to_stop(x):
    if not isinstance(x, str):
        return
    if 'N' in x:
        return
    y = translate_dna(x[:3 * int(len(x)/3)])
    if '*' in y:
        return y.split('*')[0]
    return


def to_codons(dna):
    assert len(dna) % 3 == 0
    return [dna[i * 3:(i + 1) * 3] for i in range(int(len(dna) / 3))]


def translate_to_stop(dna):
    n = len(dna) % 3
    if n != 0:
        dna = dna[:-(len(dna) % 3)]
    dna = dna.upper()
    assert len(dna) % 3 == 0
    aa = translate_dna(dna)
    if '*' in aa:
        return aa.split('*')[0]
    return aa
    # assert '*' in aa


def quick_translate(seq):
    # return str(quickdna.DnaSequence(seq).translate())
    n = len(seq) - len(seq) % 3
    return translate_dna(seq[:n])

    