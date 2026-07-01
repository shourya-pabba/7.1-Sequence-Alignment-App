import ipywidgets as widgets
from IPython.display import display, HTML, Markdown, clear_output

DNA_SIMPLE = 'simple'
PROTEIN_MATRICES = ['simple', 'blosum_like', 'pam_like', 'blosum62']

GROUPS = [set('AVLIM'), set('FYW'), set('STNQ'), set('KRH'), set('DE'), set('CGP')]

BLOSUM62_DIAG = dict(A=4,R=5,N=6,D=6,C=9,Q=5,E=5,G=6,H=8,I=4,L=4,K=5,M=5,F=6,P=7,S=4,T=5,W=11,Y=7,V=4)

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

def global_alignment(s1, s2, match=2, mismatch=3, gap=4, matrix='simple'):
    s1, s2 = clean(s1), clean(s2)
    rows, cols = len(s1)+1, len(s2)+1
    score = [[0]*cols for _ in range(rows)]
    ptr = [['stop']*cols for _ in range(rows)]
    for i in range(1, rows):
        score[i][0] = score[i-1][0] - gap; ptr[i][0] = 'up'
    for j in range(1, cols):
        score[0][j] = score[0][j-1] - gap; ptr[0][j] = 'left'
    for i in range(1, rows):
        for j in range(1, cols):
            diag = score[i-1][j-1] + sub_score(s1[i-1], s2[j-1], matrix, match, mismatch)
            up = score[i-1][j] - gap
            left = score[i][j-1] - gap
            best = max(diag, up, left)
            score[i][j] = best
            ptr[i][j] = 'diag' if best == diag else 'up' if best == up else 'left'
    return traceback(s1, s2, score, ptr, len(s1), len(s2), 'global')

def local_alignment(s1, s2, match=2, mismatch=3, gap=4, matrix='simple'):
    s1, s2 = clean(s1), clean(s2)
    rows, cols = len(s1)+1, len(s2)+1
    score = [[0]*cols for _ in range(rows)]
    ptr = [['stop']*cols for _ in range(rows)]
    best_i = best_j = best_score = 0
    for i in range(1, rows):
        for j in range(1, cols):
            diag = score[i-1][j-1] + sub_score(s1[i-1], s2[j-1], matrix, match, mismatch)
            up = score[i-1][j] - gap
            left = score[i][j-1] - gap
            best = max(0, diag, up, left)
            score[i][j] = best
            ptr[i][j] = 'stop' if best == 0 else 'diag' if best == diag else 'up' if best == up else 'left'
            if best > best_score:
                best_i, best_j, best_score = i, j, best
    return traceback(s1, s2, score, ptr, best_i, best_j, 'local')

def traceback(s1, s2, score, ptr, i, j, mode):
    end_i, end_j = i, j
    a, b = '', ''
    while i > 0 or j > 0:
        move = ptr[i][j]
        if mode == 'local' and move == 'stop': break
        if move == 'diag':
            a = s1[i-1] + a; b = s2[j-1] + b; i -= 1; j -= 1
        elif move == 'up':
            a = s1[i-1] + a; b = '-' + b; i -= 1
        elif move == 'left':
            a = '-' + a; b = s2[j-1] + b; j -= 1
        else: break
    return {'a': a, 'b': b, 'score': score[end_i][end_j], 'start1': i, 'end1': end_i, 'start2': j, 'end2': end_j, 'mode': mode}

def metrics(aln, original_a='', original_b=''):
    a, b = aln['a'], aln['b']
    matches = mismatches = gaps = gap_runs = longest = current = conservative_count = 0
    for i, (x, y) in enumerate(zip(a, b)):
        if x == '-' or y == '-':
            gaps += 1; current += 1; longest = max(longest, current)
            if i == 0 or (a[i-1] != '-' and b[i-1] != '-'): gap_runs += 1
        else:
            current = 0
            if x == y: matches += 1
            elif conservative(x, y): conservative_count += 1
            else: mismatches += 1
    length = max(1, len(a))
    return dict(score=aln['score'], identity=100*matches/length, length=len(a), matches=matches, conservative=conservative_count, mismatches=mismatches, gaps=gaps, gap_runs=gap_runs, longest_gap=longest)

def diagnose(m, mode, mismatch, gap):
    labels, warnings, story = [], [], []
    if m['identity'] >= 75: labels.append('Highly conserved')
    if m['identity'] >= 55 and m['gaps']/max(1,m['length']) < 0.18: labels.append('Biologically plausible')
    if m['gap_runs']/max(1,m['length']) > 0.09: labels.append('Fragmented')
    if m['gaps']/max(1,m['length']) > 0.25: labels.append('Gap-heavy')
    if mismatch <= 1.5: labels.append('Mismatch-tolerant')
    if mode == 'local': labels.append('Local-hit style')
    if m['gap_runs'] >= 8 and m['longest_gap'] <= 3: warnings.append('Many short scattered gaps: mathematically legal, but biologically suspicious.')
    if mismatch <= 1.5 and m['mismatches']/max(1,m['length']) > 0.22: warnings.append('High mismatch tolerance may force unrelated regions to align.')
    if m['identity'] >= 70: story.append('High conservation suggests evolutionary similarity or functional constraint.')
    if m['gap_runs'] > 6: story.append('Many isolated gaps may be less plausible than one insertion/deletion event.')
    if mode == 'local': story.append('Local alignment captures the strongest conserved region rather than forcing end-to-end alignment.')
    if not story: story.append('This alignment needs biological judgment beyond its score.')
    return labels, warnings, story
