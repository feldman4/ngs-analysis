import gzip
from collections import defaultdict
from glob import glob

import numpy as np
import pandas as pd
from natsort import natsorted

from .constants import CODONS, CODONS_REVERSE


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
    return ''.join([CODONS[s[i*3:(i+1)*3]] for i in range(int(len(s)/3))])


def make_equivalent_codons(dna_to_aa):
    """Make dictionary from codon to other codons for the same amino acid
    """
    aa_to_dna = defaultdict(list)
    for codon, aa in dna_to_aa.items():
        aa_to_dna[aa] += [codon]
    
    equivalent_codons = {}
    for codon in dna_to_aa:
        aa = dna_to_aa[codon]
        equivalent_codons[codon] = list(set(aa_to_dna[aa]) - {codon})
            
    return equivalent_codons


equivalent_codons = make_equivalent_codons(CODONS)


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


def quality_scores_to_array(quality_scores, baseline=ord('!')):
    """Only works if all quality scores have equal length.
    Expects strings not bytes.
    """
    q = np.array(quality_scores)
    return (np.array(q).astype(f'S{len(q[0])}')
              .view('S1').view(np.uint8)
              .reshape(len(q), len(q[0]))
               - baseline
              )


def aa_to_dna_re(aa_seq):
    return ''.join(f'(?:{CODONS_REVERSE[x]})' for x in aa_seq)


def make_kmer_dict(sequences, k):
    """
    """
    kmers = defaultdict(list)
    for i, seq in enumerate(sequences):
        for kmer in get_kmers(seq, k):
            kmers[kmer].append(i)
    return kmers


def match_nearest(query, sequences, kmers):
    from Levenshtein import distance

    k = len(next(iter(kmers.keys())))

    candidates = []
    for kmer in get_kmers(query, k):
        candidates.extend(kmers[kmer])
    candidates = set(candidates)
    # guess
    candidates = sorted(candidates, key=lambda i: ~
                        sequences[i].startswith(query[:2]))

    matches = []
    for i in candidates:
        d = distance(sequences[i], query)
        matches.append((d, i))
        # exact match
        if d == 0:
            break
    d, i = sorted(matches)[0]
    return d, i


def match_queries(queries, sequences, window, k, progress=lambda x: x):
    """Match queries to reference sequences based on Levenshtein distance between
    prefixes of length `window`. Only pairs with a shared kmer of length `k` are
    checked. For each query, finds the first nearest prefix and returns all sequences 
    that share that prefix.
    """
    query_lookup = {x: x[:window] for x in queries}
    query_prefixes = sorted(set([x[:window] for x in queries]))

    ref_lookup = defaultdict(list)
    for x in sequences:
        ref_lookup[x[:window]].append(x)
    ref_prefixes = sorted(set([x[:window] for x in sequences]))

    kmers = make_kmer_dict(ref_prefixes, k)

    hits = {}
    for q in progress(query_prefixes):
        try:
            hits[q] = match_nearest(q, ref_prefixes, kmers)
        except IndexError:
            pass

    results = []
    for q in queries:
        try:
            d, i = hits[query_lookup[q]]
            results.append(ref_lookup[ref_prefixes[i]])
        except KeyError:
            results.append([])
    return results


def add_design_matches(df_reads, col, reference, window, k):
    """
    `df_reads` is a dataframe containing `col` with sequences
    `reference` is a list of references
    """
    queries = df_reads[col].fillna('').pipe(list)
    queries = [q if '*' not in q else '' for q in queries]
    results = match_queries(queries, reference, window, k)

    df_reads = df_reads.copy()
    design_distance, design_match, design_equidistant = zip(
        *calculate_distance_matches(queries, results))
    return (df_reads
            .assign(design_match=design_match, design_distance=design_distance,
                    design_equidistant=design_equidistant)
            )


def calculate_distance_matches(queries, results):
    """Get columns `design_distance` and `design_match` from results of `match_queries`.
    """
    from Levenshtein import distance

    arr = []
    for q, rs in zip(queries, results):
        if len(rs) == 0:
            arr += [(-1, '', 0)]
        else:
            ds = [(distance(q, r), r) for r in rs]
            d, s = sorted(ds)[0]
            degeneracy = sum([x[0] == d for x in ds])
            arr += [(d, s, degeneracy)]
    return arr


def match_and_check(queries, reference, window, k, ignore_above=40, progress=lambda x: x):
    """Perform fast Levenshtein distance matching of queries to reference
    and check the results by brute force calculation (all pairs). Mismatches
    with an edit distance greater than `ignore_above` are ignored.
    """
    from Levenshtein import distance

    print(f'Matching {len(queries)} queries to {len(reference)} '
          f'reference sequences, window={window} and k={k}')
    df_matched = (pd.DataFrame({'sequence': queries})
                  .pipe(add_design_matches, col='sequence',
                        reference=reference, window=window, k=k))
    it = (df_matched
          [['sequence', 'design_match']].values)
    print('Checking fast matches against brute force matches...')
    different = 0
    for seq, match in progress(it):
        xs = sorted(reference, key=lambda x: distance(x, seq))
        a, b = distance(seq, match), distance(seq, xs[0])
        if b < a and b <= ignore_above:
            different += 1
            print(f'{a},{b} (fast,exact distance); query={seq}')
    print(f'Total mismatches: {different}')
    return df_matched


def try_translate_dna(s):
    try:
        return translate_dna(s)
    except:
        return None


def findone(aa, dna):
    """Simple case of amino acid substring in in-frame DNA.
    """
    aa_ = translate_dna(dna)
    i = aa_.index(aa)
    return dna[i * 3:(i + len(aa)) * 3]


def select_most_different(xs, n):
    """Quickly select sequences with high mutual Levenshtein distance.
    """
    from Levenshtein import distance
    xs = list(xs)
    arr = [xs.pop()]
    for _ in range(n - 1):
        new = sorted(xs, key=lambda x: -min(distance(x, y) for y in arr))[0]
        xs.remove(new)
        arr += [new]
    return arr


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


def pairwise_levenshtein(seqs):
    """Calculate distance matrix for a list of sequences.
    """
    from Levenshtein import distance
    arr = []
    for i, a in enumerate(seqs):
        for j, b in enumerate(seqs[i+1:]):
            arr += [(i, i + j + 1, distance(a,b))]
    n = len(seqs)
    D = np.zeros((n, n), dtype=int)
    i, j, d = zip(*arr)
    D[i, j] = d
    return D + D.T


def kmers_unique(s, k):
    kmers = get_kmers(s, k)
    return len(kmers) == len(set(kmers))


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


def find_aa_in_dna(aa, dna):
    """Simple case of amino acid substring in DNA.
    """
    dna = dna[:-(len(dna) % 3)].upper()
    aa_ = translate_dna(dna)
    i = aa_.index(aa)
    return dna[i * 3:(i + len(aa)) * 3]


def quick_translate(seq):
    # return str(quickdna.DnaSequence(seq).translate())
    n = len(seq) - len(seq) % 3
    return translate_dna(seq[:n])

    