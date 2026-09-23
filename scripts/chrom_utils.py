
import re


MITO_PATTERNS = frozenset({'chrM', 'MT', 'chrMT', 'chrm', 'M'})
CHRX_PATTERNS = frozenset({'chrX', 'X', 'chrx', 'ChrX', 'CM039033.1'})
CHRY_PATTERNS = frozenset({'chrY', 'Y', 'chry', 'ChrY'})
SEX_CHROM_PATTERNS = CHRX_PATTERNS | CHRY_PATTERNS


_CHR_TO_CM039 = {
    'chr1':  'CM039011.1', 'chr2':  'CM039012.1', 'chr3':  'CM039013.1',
    'chr4':  'CM039014.1', 'chr5':  'CM039015.1', 'chr6':  'CM039016.1',
    'chr7':  'CM039017.1', 'chr8':  'CM039018.1', 'chr9':  'CM039019.1',
    'chr10': 'CM039020.1', 'chr11': 'CM039021.1', 'chr12': 'CM039022.1',
    'chr13': 'CM039023.1', 'chr14': 'CM039024.1', 'chr15': 'CM039025.1',
    'chr16': 'CM039026.1', 'chr17': 'CM039027.1', 'chr18': 'CM039028.1',
    'chr19': 'CM039029.1', 'chr20': 'CM039030.1', 'chr21': 'CM039031.1',
    'chr22': 'CM039032.1', 'chrX':  'CM039033.1',
}
_CM039_TO_CHR = {v: k for k, v in _CHR_TO_CM039.items()}
_CHR_ALIAS = {**_CHR_TO_CM039, **_CM039_TO_CHR}



def is_mito(name: str) -> bool:
    return name in MITO_PATTERNS


def is_chrx(name: str) -> bool:
    return name in CHRX_PATTERNS


def is_chry(name: str) -> bool:
    return name in CHRY_PATTERNS


def is_sex_chrom(name: str) -> bool:
    return name in SEX_CHROM_PATTERNS



def resolve_chrom(name: str, target_chroms: set) -> str:
    if name in target_chroms:
        return name
    alias = _CHR_ALIAS.get(name)
    if alias and alias in target_chroms:
        return alias
    return None


def display_label(name: str) -> str:
    return _CM039_TO_CHR.get(name, name)



def natural_chrom_key(name: str) -> tuple:
    n = name.lower().replace('chr', '')
    if n in ('x',):
        return (23, 0, name)
    if n in ('y',):
        return (24, 0, name)
    if n in ('m', 'mt'):
        return (25, 0, name)
    m = re.match(r'(\d+)', n)
    if m:
        return (int(m.group(1)), 0, name)
    m = re.match(r'cm0*(\d+)', n)
    if m:
        return (int(m.group(1)), 0, name)
    return (100, 0, name)
