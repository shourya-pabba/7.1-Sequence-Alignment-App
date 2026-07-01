from scoring import sub_score, clean


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
