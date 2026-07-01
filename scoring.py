from constants import GROUPS, BLOSUM62_DIAG

DNA_LETTERS = set('ACGTURYKMSWBDHVN')


def conservative(a, b, sequence_type='protein'):
    a = a.upper()
    b = b.upper()
    if a == b or a == '-' or b == '-':
        return False
    if sequence_type != 'protein':
        return False
    return any(a in g and b in g for g in GROUPS)


def sub_score(a, b, matrix, match, mismatch, sequence_type='protein'):
    if matrix == 'simple':
        return match if a == b else -mismatch
    if a == b:
        return BLOSUM62_DIAG.get(a, match) if matrix == 'blosum62' else match
    if conservative(a, b, sequence_type=sequence_type):
        return 2 if matrix in ['blosum_like', 'blosum62'] else 1
    return -mismatch


def clean(seq):
    return ''.join(ch for ch in seq.upper() if ch.isalpha())
