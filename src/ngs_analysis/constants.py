from collections import defaultdict
import os
from pathlib import Path
import string

try:
    HOME = Path(os.environ['HOME'])
except KeyError:
    HOME = Path('/home/dfeldman')

JOBLIB_CACHE = HOME / '.joblib'

RESOURCES = Path(__file__).parents[0] / 'resources'
RULE_SETS = RESOURCES / 'rule_sets.csv'
VISITOR_FONT_PATH = RESOURCES / 'visitor1.ttf'

BQ_PROJECT_ID = 'bilf-350112'

GO_TERM = 'GO_term'
GO = 'GO ID'
GO_SYNONYM = 'DB Object Synonym (|Synonym)'
GO_SYMBOL = 'DB Object Symbol'
GO_TERM_COUNTS = 'GO_term_counts'
SEARCH_KEYWORD = 'keyword'
GO_SYNONYM = 'DB Object Synonym (|Synonym)'

GENE_ID = 'gene_id'
GENE_SYMBOL = 'gene_symbol'

HGNC = 'HGNC'
UNIPROTKB = 'UniProtKB'
ENSG = 'ENSG'
GENE_ALIAS = 'gene_alias'
RCSB = 'RCSB'

biomart_columns = {'Gene stable ID': ENSG,
                   'HGNC ID': HGNC,
                   'NCBI gene ID': GENE_ID,
                   'NCBI gene (formerly Entrezgene) ID': GENE_ID,
                   'HGNC symbol': GENE_SYMBOL,
                   'UniProtKB Gene Name ID': UNIPROTKB,
                  }

PTHR = 'PTHR'
PTHR_SF = 'PTHR_SF'
PTHR_FAMILY = 'PTHR_family'
PTHR_SUBFAMILY = 'PTHR_subfamily'
PTHR_CLASS_LIST = 'PC_list'

pthr_columns = {
    0: 'identifier',
    2: PTHR_SF,
    3: PTHR_FAMILY,
    4: PTHR_SUBFAMILY,
    8: PTHR_CLASS_LIST,
}

MZ_DOUBLE_SPACING = 0.5001917279701898

AA_3 = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
           'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']
AA_1 = list('ARNDCQEGHILKMFPSTWYV')
CANONICAL_AA = AA_1
AA_3_1 = dict(zip(AA_3, AA_1))
AA_1_3 = dict(zip(AA_1, AA_3))

skyline_columns = {
    'Replicate': 'sample', 
    'Replicate Name': 'sample',
    'Protein': 'short_name', 
    'Peptide': 'sequence',
    'Peptide Retention Time': 'RTime',
    'Normalized Area': 'peak_area',
    'Total Area MS1': 'ms1_area',
    'Best Retention Time': 'RTime',
    'Min Start Time': 'RTime_start',
    'Max End Time': 'RTime_end',
    'Average Mass Error PPM': 'mass_error_ppm',
    'Isotope Dot Product': 'idotp',
    }


PT02_BACKBONE_START = 'ATTCTCCTTGGAATTTGCCCTTTTTGAGTTTGGATCTTGGTTCAT'
iON_BACKBONE_START = 'GACATTGATTATTGACTAGTTATTAATAGTAATCAAT'
pET_BACKBONE_START = 'TAATACGACTCACTATAGGGGAATTGTGAGCGGATAACAATTCC'

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

CHAIN_ALPHABET = string.ascii_uppercase + string.ascii_lowercase

SLUGIFY_REPLACEMENTS = [
    ('<', 'lt'),
    ('<=', 'lte'),
    ('>', 'gt'),
    ('>=', 'gte'),
    ('&', 'and'),
    ('|', 'or'),
]