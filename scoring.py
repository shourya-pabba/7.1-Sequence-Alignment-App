from constants import GROUPS, BLOSUM62_DIAG


def conservative(a, b):
    return a != b and a != '-' and b != '-' and any(a in g and b in g for g in GROUPS)


def sub_score(a, b, matrix, match, mismatch):
    if matrix == 'simple':
        return match if a == b else -mismatch
    if a == b:
        return BLOSUM62_DIAG.get(a, match) if matrix == 'blosum62' else match
    if conservative(a, b):
        return 2 if matrix in ['blosum_like', 'blosum62'] else 1
    return -mismatch


def clean(seq):
    return ''.join(ch for ch in seq.upper() if ch.isalpha())
