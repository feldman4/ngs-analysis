import gzip
import os
import re
from collections import defaultdict
from glob import glob

import numpy as np
import pandas as pd
from natsort import natsorted

from .constants import RESOURCES, CODONS, CODONS_REVERSE


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
    with open(filename, 'w') as fh:
        fh.write(format_fake_fastq(list_or_records))


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


def load_codons(organism):
    f = os.path.join(RESOURCES, 'codon_usage', 'organisms.csv')
    taxids = pd.read_csv(f).set_index('organism')['taxid'].to_dict()
    
    if organism.lower() == 'yeast':
        organism = 's_cerevisiae'
    organism = organism.lower().replace('.', '').replace(' ', '_')
    
    try:
        table = f'{organism}_{taxids[organism]}.csv'
    except KeyError:
        raise ValueError(f'{organism} must be one of {list(taxids.keys())}')
    f = os.path.join(RESOURCES, 'codon_usage', table)
    return (pd.read_csv(f)
            .assign(codon_dna=lambda x: x['codon'].str.replace('U', 'T')))


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


def sanger_database(drive):
    df_sanger = drive.get_excel('cloning/sanger')

    extra_cols = [x for x in df_sanger.columns
                  if x not in ('identifier', 'search')]

    arr = []
    for _, row in df_sanger.iterrows():
        files = natsorted(glob(row['search']))
        if len(files) == 0:
            print(f'No files found from row {row.to_dict()}')
            continue
        (pd.DataFrame({'file': files})
         .assign(**{x: row[x] for x in extra_cols})
         .assign(name=lambda x: x['file'].str.extract(row['identifier']))
         .assign(seq=lambda x: x['file'].apply(read_ab1))
         .assign(seq_rc=lambda x: x['seq'].apply(reverse_complement))
         .pipe(arr.append)
         )

    cols = extra_cols + ['name', 'file', 'seq', 'seq_rc']
    return pd.concat(arr)[cols]


def read_ab1(f):
    from Bio import SeqIO
    with open(f, 'rb') as fh:
        records = list(SeqIO.parse(fh, 'abi'))
        assert len(records) == 1
        seq = str(records[0].seq)
    return seq


def print_alignment(a, b, width=60, as_string=False):
    """Levenshtein alignment.
    """
    import edlib
    alignment = edlib.align(a, b, task='path')
    d = edlib.getNiceAlignment(alignment, a, b)
    lines = []
    for i in range(0, max(map(len, d.values())), width):
        lines += [str(i)]
        for x in d.values():
            lines += [x[i:i+width]]

    txt = '\n'.join(lines)
    if as_string:
        return txt
    else:
        print(txt)


def reverse_translate_max(aa_seq, organism='e_coli'):
    if organism not in codon_maps:
        codon_maps[organism] = (load_codons(organism)
        .sort_values('relative_frequency', ascending=False)
        .drop_duplicates('amino_acid')
        .set_index('amino_acid')['codon_dna'].to_dict()
        ) 
    codon_map = codon_maps[organism]
    return ''.join([codon_map[x] for x in aa_seq])


def reverse_translate_random(aa_seq, organism='e_coli', rs='input', cutoff=0.12):
    if rs == 'input':
        seed = hash(aa_seq) % 10**8
        rs = np.random.RandomState(seed=seed)
    if (organism, cutoff) not in codon_maps:
        codon_maps[(organism, cutoff)] = (load_codons(organism)
        .query('relative_frequency > @cutoff')
        .groupby('amino_acid')['codon_dna'].apply(list).to_dict()
        )
    codon_map = codon_maps[(organism, cutoff)]
    return ''.join([rs.choice(codon_map[x]) for x in aa_seq])


def get_genbank_features(f, error_on_repeat=True):
    from Bio import SeqIO
    records = list(SeqIO.parse(open(f,'r'), 'genbank'))
    if len(records) != 1:
        raise ValueError(f'found {len(records)} records in genbank {f}')

    features = {}
    for f in records[0].features:
        label = f.qualifiers['label'][0]
        seq = f.extract(records[0].seq)
        if label in features and features[label] != seq and error_on_repeat:
            raise ValueError(f'repeated feature {label}')
        features[label] = str(seq)
    return features


def rolling_gc(seq, window):
    from scipy.ndimage.filters import convolve
    gc_window = 20
    return convolve(np.array([x in 'GC' for x in seq])*1., 
                          np.ones(window)/window, 
                          mode='reflect')


def to_codons(seq):
    assert len(seq) % 3 == 0
    return [seq[i*3:(i+1)*3] for i in range(int(len(seq)/3))]


def codon_adaptation_index(seq, organism='e_coli'):
    return (load_codons(organism)
        .assign(w=lambda x: x.groupby('amino_acid')['relative_frequency']
                .transform(lambda y: y / y.max()))
        .set_index('codon_dna')
        .loc[to_codons(seq)]['w']
        .pipe(lambda x: np.prod(x)**(1/len(x)))
           )


def compare_sequences(sequences, window=25, k=6):
    import matplotlib.pyplot as plt

    fig, (ax0, ax1) = plt.subplots(figsize=(12, 4), ncols=2)
    for name in sequences:
        cai = codon_adaptation_index(sequences[name])
        mean_gc = np.mean([x in 'GC' for x in sequences[name]])
        gc_trace = rolling_gc(sequences[name], window)
        label = f'{name}: avg={mean_gc:.2g} std={np.std(gc_trace):.2g} cai={cai:.2g}'
        ax0.plot(gc_trace, label=label)
    
        (pd.Series(get_kmers(sequences[name], k))
         .value_counts().value_counts().sort_index()
         .pipe(lambda x: x/x.sum())
         .plot(ax=ax1, marker='.', ms=10, label=name))
        
    
    ax0.plot([0, len(gc_trace)], [0.5, 0.5], color='gray', lw=1, ls='--', zorder=-1)
    ax0.legend()
    ax0.set_title(f'average GC content over {window} nt window')
    ax0.set_ylabel('GC fraction')
    ax0.set_xlabel('DNA sequence position')
    ax0.set_ylim([0.25, 0.85])
    ax0.set_xlim([0, len(gc_trace)])
    
    ax1.set_title(f'repeated kmers with k={k}')
    ax1.set_xlabel('kmer count')
    ax1.set_ylabel(f'fraction of kmers')
    ax1.set_xticks([1,2,3,4,5])
    ax1.legend()

    return fig


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


def load_abi_zip(filename):
    """Extract Bio.SeqIO records from sanger zip file.
    """
    import zipfile
    from io import BytesIO

    from Bio import SeqIO
    zh = zipfile.ZipFile(filename, 'r')
    arr = []
    for zi in zh.filelist:
        if not zi.filename.endswith('ab1'):
            print(f'Skipping {zi.filename}')
            continue
        fh = zh.open(zi.filename, 'r')
        buffer = BytesIO(fh.read())
        arr += [SeqIO.read(buffer, 'abi')]
    return arr


def get_abi_traces(abi_record):
    channels = 'DATA9', 'DATA10', 'DATA11', 'DATA12'
    bases = list('GATC')
    traces = np.array([abi_record.annotations['abif_raw'][c]
                       for c in channels])
    df = pd.DataFrame(traces).T
    df.columns = bases
    return df


def try_translate_dna(s):
    try:
        return translate_dna(s)
    except:
        return None


def digest_protein_fasta(filename, digest_pat='[R|K]'):
    """Convert protein fasta into a table of peptides, digesting with provided regex.
    """
    records = read_fasta(filename)

    arr = []
    for name, seq in records:
        parts = re.split(f'({digest_pat})', seq)
        parts = [x for x in parts if x]
        peptides = [a+b for a,b in zip(parts[::2], parts[1::2])]
        for peptide in peptides:
            arr += [{'name': name, 'sequence': peptide}]

    return pd.DataFrame(arr)


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


def edge_primers(dna, melt_temp): 
    from Bio.SeqUtils.MeltingTemp import Tm_NN
    dna_rev = reverse_complement(dna)

    fwd = dna[:10]
    while Tm_NN(fwd) < melt_temp:
        fwd = dna[:len(fwd) + 1]
        
    rev = dna_rev[:10]
    while Tm_NN(rev) < melt_temp:
        rev = dna_rev[:len(rev) + 1]
        
    return fwd, rev


def find_orfs(seq, kozak='GCCACC', return_aa=True):
    """Find open reading frames starting with kozak in circular 
    DNA sequence.
    """
    plasmid = str(seq)
    pat = f'{kozak}(ATG(?:...)*?)(?:TAA|TAG|TGA)'
    ref = plasmid + plasmid
    orfs = re.findall(pat, ref) + re.findall(pat, reverse_complement(ref))
    if return_aa:
        orfs = [translate_dna(x) for x in orfs]
    return sorted(set(orfs), key=lambda x: -1 * len(x))


def find_longest_orf(dna, kozak='GCCACC'):
    pat = f'{kozak}(ATG(?:...)*?)(?:TAA|TGA|TAG)'
    orfs = re.findall(pat, dna.upper())
    if orfs:
        return sorted(orfs, key=len)[-1]


def add_features(record, features):
    """Add features to `record` based on name=>DNA dictionary `features`.
    A feature with "orf" in the name is given a special annotation that 
    shows up as a translation in benchling.
    """

    vector = str(record.seq).upper()

    n = len(vector)

    features = dict(features)
    arr = []
    for name, feature in features.items():
        feature = feature.upper()
        m = len(feature)

        for strand in (-1, 1):
            key = reverse_complement(feature) if strand == -1 else feature
            starts = [x.start() for x in re.finditer(key, vector * 2)] 
            starts = [x for x in starts if x < n]
            for start in starts:
                end = start + m
                qualifiers = dict(label=name)
                feature_type = 'misc'
                if 'orf' in name:
                    feature_type = 'CDS'
                    qualifiers['translation'] = translate_dna(key) + '*'
                    end += strand * 3

                if end < n:
                    location = FeatureLocation(start, end, strand=strand)
                else:
                    f1 = FeatureLocation(start, n, strand=strand)
                    f2 = FeatureLocation(0, end % n, strand=strand)
                    location = f1 + f2
                arr += [SeqFeature(location, type=feature_type, qualifiers=qualifiers)]

    # copies features but not annotations
    new_record = record[:]
    new_record.annotations = record.annotations
    new_record.features += arr
    return new_record


def translate_to_stop(x):
    if not isinstance(x, str):
        return
    if 'N' in x:
        return
    y = translate_dna(x[:3 * int(len(x)/3)])
    if '*' in y:
        return y.split('*')[0]
    return


def ion_plasmid_integration(filename):
    """Generate plasmid map after integration for an iON vector
    :param filename: iON vector genbank
    """
    left_itr = 'CCCTAGAAAGATAGTCTGCGTAAAATTGACGCATG'.upper()
    right_itr = 'ccctagaaagataatcatattgtgacgtacgttaaagataatcatgcgtaaaattgacgcatg'.upper()

    record = list(SeqIO.parse(filename, 'genbank'))[0]
    original_features = get_genbank_features(filename, error_on_repeat=False)
    # remove annoying small features
    original_features = {k: v for k,v in original_features.items() if 6 < len(v)}
    
    dna = str(record.seq)
    start = dna * 2
    # in iON figure 1A notation
    # partial, one, two, one, partial
    _, one, two, _, _ = split_by_regex(f'{left_itr}|{right_itr}', start, keep='prefix')

    # deduplicate
    one, two = ['AA' + x[:-2] for x in (one, two)]
    result = two + reverse_complement(one)
    
    result_record = SeqRecord(Seq(result), name=record.name + '_integrated', 
                              annotations={'molecule_type': 'DNA'})    
    result_record = add_features(result_record, original_features)
    
    return result_record


def remove_restriction_sites(dna, sites, rs):
    codon_sequence = to_codons(dna)
    sites = set(sites)
    for site in list(sites):
        sites.add(reverse_complement(site))

    for site in sites:
        width = len(site)
        dna = ''.join(codon_sequence)
        if site not in dna:
            continue

        for i in range(len(dna)):
            # if we encounter a restriction site
            if dna[i:i + len(site)] == site:
                # change any of these codons
                overlapped_codons = sorted(
                    set([int((i + offset) / 3) for offset in range(width)]))
                # accept first change that removes restriction site
                for j in overlapped_codons:
                    # change this codon
                    new_codon = swap_codon(codon_sequence[j], rs)
                    local_dna = ''.join([new_codon if k == j else codon_sequence[k]
                                         for k in overlapped_codons])
                    # if codon removes this site, keep it
                    if site not in local_dna:
                        codon_sequence = codon_sequence[:j] + \
                            [new_codon] + codon_sequence[j + 1:]
                        break
    dna = ''.join(codon_sequence)
    for site in sites:
        assert site not in dna
            
    return dna


def restriction_sites_in_dna(dna, sites):
    sites = set(list(sites) + [reverse_complement(x) for x in sites])
    for site in sites:
        if site in dna:
            return True
    return False


def to_codons(dna):
    assert len(dna) % 3 == 0
    return [dna[i * 3:(i + 1) * 3] for i in range(int(len(dna) / 3))]


def swap_codon(codon, rs):
    """Swap codon at random, if possible.
    """
    options = equivalent_codons[codon]
    if not options:
        return codon
    return rs.choice(options)


def remove_restriction_sites_dnachisel(cds, avoid, rs=None):
    if not restriction_sites_in_dna(cds, avoid):
        return cds

    import dnachisel as dc
    
    if rs is None:
        np.random.seed(0)
    else:
        np.random.seed(int(1e9*rs.rand()))
        
    n = len(cds)

    constraints = [dc.AvoidPattern(x, location=0) for x in avoid]
    constraints += [dc.EnforceTranslation(location=(0, n))]

    problem = dc.DnaOptimizationProblem(
        sequence=cds,
        constraints=constraints,
        objectives=[dc.AvoidChanges()], 
        logger=None
    )
    problem.max_iters = 50
    problem.resolve_constraints()
    problem.optimize()
    
    assert not restriction_sites_in_dna(problem.sequence, avoid)
    return problem.sequence


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


def generate_adapters(priming_length=12, enzyme='BsmBI', tm_range=2, target_GC=0.5,
                      max_calculate=10_000, num_sample=100_000, distance_threshold=6, seed=0):
    """Tm is calculated for priming site plus restriction site. A target Tm is defined
    based on the average of sequences within 10% of target_GC.
    """
    from Levenshtein import distance
    from scipy.sparse import csr_matrix
    from .utils import maxy_clique_groups
    from Bio.SeqUtils import MeltingTemp as mt

    enzymes = {'BsmBI': 'CGTCTC'}
    site = enzymes[enzyme]
    avoid = 'AAA', 'GGG', 'CCC', 'TTT', site, reverse_complement(site)
    rs = np.random.RandomState(seed)
    candidates = [''.join(x) for x in rs.choice(list('ACTG'), size=(num_sample, priming_length))]
    candidates = [x for x in candidates if not any(y in x for y in avoid)]

    df_adapters = pd.DataFrame({'sequence': candidates})

    df_adapters['Tm'] = (df_adapters['sequence'] + site).apply(mt.Tm_NN)
    df_adapters['GC'] = df_adapters['sequence'].str.count('G|C') / df_adapters['sequence'].str.len()

    mean_tm = df_adapters.query('abs(GC - @target_GC) < 0.1')['Tm'].mean()
    df_adapters = df_adapters.query('abs(Tm - @mean_tm) < @tm_range/2')
    
    seqs = df_adapters['sequence'][:max_calculate]
    print(f'Calculating pairwise distances for {len(seqs):,} priming sites')
    D = pairwise_levenshtein(seqs)
    
    print(f'Selecting priming sites with pairwise distance >= {distance_threshold}')
    
    A = distance_threshold <= D
    cm = csr_matrix(1-A, shape=A.shape)
    centers = np.array(maxy_clique_groups(cm, [1] * A.shape[0]))
    
    return df_adapters['sequence'].iloc[centers].pipe(list)


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


def sequence_sap_score(x):
    """Sum of SAP weights for all residues, regardless of solvent exposure.
    """
    from .ppi.constants import SAP_WEIGHTS
    return sum([SAP_WEIGHTS[aa] for aa in x])