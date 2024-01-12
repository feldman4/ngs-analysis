from collections import defaultdict

# user-provided
config_yaml = 'config.yaml' # how to parse and compare
sample_table = 'samples.csv' # NGS samples
reference_dna_table = 'reference_dna.csv' # designed insert DNA
sample_plan_table = 'sample_plan.csv' # cloning plan

working_directories = '0_paired_reads', '1_reads', '2_parsed', '3_mapped', '3_mapped/map'

# generated or provided
design_table = 'designs.csv' # designed DNA or protein parts

ngmerge = 'NGmerge'
bwa = 'bwa'
mmseqs = 'mmseqs'

MMSEQS_KMER_AA = 6
MMSEQS_KMER_DNA = 6 # TODO: problematic for short sequences?


AA_3 = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
           'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']
AA_1 = list('ARNDCQEGHILKMFPSTWYV')
CANONICAL_AA = AA_1
AA_3_1 = dict(zip(AA_3, AA_1))
AA_1_3 = dict(zip(AA_1, AA_3))

CODONS = {
    'TAA': '*',
    'TAG': '*',
    'TGA': '*',
    'GCA': 'A',
    'GCC': 'A',
    'GCG': 'A',
    'GCT': 'A',
    'TGC': 'C',
    'TGT': 'C',
    'GAC': 'D',
    'GAT': 'D',
    'GAA': 'E',
    'GAG': 'E',
    'TTC': 'F',
    'TTT': 'F',
    'GGA': 'G',
    'GGC': 'G',
    'GGG': 'G',
    'GGT': 'G',
    'CAC': 'H',
    'CAT': 'H',
    'ATA': 'I',
    'ATC': 'I',
    'ATT': 'I',
    'AAA': 'K',
    'AAG': 'K',
    'CTA': 'L',
    'CTC': 'L',
    'CTG': 'L',
    'CTT': 'L',
    'TTA': 'L',
    'TTG': 'L',
    'ATG': 'M',
    'AAC': 'N',
    'AAT': 'N',
    'CCA': 'P',
    'CCC': 'P',
    'CCG': 'P',
    'CCT': 'P',
    'CAA': 'Q',
    'CAG': 'Q',
    'AGA': 'R',
    'AGG': 'R',
    'CGA': 'R',
    'CGC': 'R',
    'CGG': 'R',
    'CGT': 'R',
    'AGC': 'S',
    'AGT': 'S',
    'TCA': 'S',
    'TCC': 'S',
    'TCG': 'S',
    'TCT': 'S',
    'ACA': 'T',
    'ACC': 'T',
    'ACG': 'T',
    'ACT': 'T',
    'GTA': 'V',
    'GTC': 'V',
    'GTG': 'V',
    'GTT': 'V',
    'TGG': 'W',
    'TAC': 'Y',
    'TAT': 'Y',
}

CODONS_REVERSE = defaultdict(list)
[CODONS_REVERSE[v].append(k) for k, v in CODONS.items()]
CODONS_REVERSE = {k: '|'.join(v) for k, v in CODONS_REVERSE.items()}

SLUGIFY_REPLACEMENTS = [
    ('<', 'lt'),
    ('<=', 'lte'),
    ('>', 'gt'),
    ('>=', 'gte'),
    ('&', 'and'),
    ('|', 'or'),
]