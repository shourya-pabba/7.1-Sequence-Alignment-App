from scoring import sub_score, clean


def global_alignment(s1, s2, match=2, mismatch=3, gap=4, matrix='simple', gap_open=None, gap_extend=None):
    s1, s2 = clean(s1), clean(s2)
    if gap_open is None:
        gap_open = gap
    if gap_extend is None:
        gap_extend = gap

    rows, cols = len(s1)+1, len(s2)+1
    match_scores = [[float('-inf')]*cols for _ in range(rows)]
    gap_a_scores = [[float('-inf')]*cols for _ in range(rows)]
    gap_b_scores = [[float('-inf')]*cols for _ in range(rows)]
    ptr = [['stop']*cols for _ in range(rows)]

    match_scores[0][0] = 0
    for j in range(1, cols):
        gap_a_scores[0][j] = -gap_open - (j-1)*gap_extend
        ptr[0][j] = 'gap_a'
    for i in range(1, rows):
        gap_b_scores[i][0] = -gap_open - (i-1)*gap_extend
        ptr[i][0] = 'gap_b'

    for i in range(1, rows):
        for j in range(1, cols):
            diag = match_scores[i-1][j-1] + sub_score(s1[i-1], s2[j-1], matrix, match, mismatch)
            gap_a = max(match_scores[i][j-1] - gap_open, gap_a_scores[i][j-1] - gap_extend)
            gap_b = max(match_scores[i-1][j] - gap_open, gap_b_scores[i-1][j] - gap_extend)
            best = max(diag, gap_a, gap_b)
            if best == diag:
                match_scores[i][j] = best
                ptr[i][j] = 'diag'
            elif best == gap_a:
                gap_a_scores[i][j] = best
                ptr[i][j] = 'gap_a'
            else:
                gap_b_scores[i][j] = best
                ptr[i][j] = 'gap_b'

    end_state = 'diag'
    if gap_a_scores[-1][-1] > match_scores[-1][-1] and gap_a_scores[-1][-1] >= gap_b_scores[-1][-1]:
        end_state = 'gap_a'
    elif gap_b_scores[-1][-1] > match_scores[-1][-1] and gap_b_scores[-1][-1] >= gap_a_scores[-1][-1]:
        end_state = 'gap_b'

    return traceback(s1, s2, match_scores, gap_a_scores, gap_b_scores, ptr, len(s1), len(s2), 'global', end_state)


def local_alignment(s1, s2, match=2, mismatch=3, gap=4, matrix='simple', gap_open=None, gap_extend=None):
    s1, s2 = clean(s1), clean(s2)
    if gap_open is None:
        gap_open = gap
    if gap_extend is None:
        gap_extend = gap

    rows, cols = len(s1)+1, len(s2)+1
    match_scores = [[0]*cols for _ in range(rows)]
    gap_a_scores = [[0]*cols for _ in range(rows)]
    gap_b_scores = [[0]*cols for _ in range(rows)]
    ptr = [['stop']*cols for _ in range(rows)]
    best_i = best_j = best_score = 0

    for i in range(1, rows):
        for j in range(1, cols):
            diag = match_scores[i-1][j-1] + sub_score(s1[i-1], s2[j-1], matrix, match, mismatch)
            gap_a = max(match_scores[i][j-1] - gap_open, gap_a_scores[i][j-1] - gap_extend)
            gap_b = max(match_scores[i-1][j] - gap_open, gap_b_scores[i-1][j] - gap_extend)
            best = max(0, diag, gap_a, gap_b)
            if best == 0:
                match_scores[i][j] = 0
                gap_a_scores[i][j] = 0
                gap_b_scores[i][j] = 0
                ptr[i][j] = 'stop'
            elif best == diag:
                match_scores[i][j] = best
                ptr[i][j] = 'diag'
            elif best == gap_a:
                gap_a_scores[i][j] = best
                ptr[i][j] = 'gap_a'
            else:
                gap_b_scores[i][j] = best
                ptr[i][j] = 'gap_b'
            if best > best_score:
                best_i, best_j, best_score = i, j, best

    return traceback(s1, s2, match_scores, gap_a_scores, gap_b_scores, ptr, best_i, best_j, 'local', ptr[best_i][best_j])


def traceback(s1, s2, match_scores, gap_a_scores, gap_b_scores, ptr, i, j, mode, end_state):
    end_i, end_j = i, j
    a, b = '', ''
    state = end_state
    while i > 0 or j > 0:
        if mode == 'local' and state == 'stop':
            break
        if state == 'diag':
            a = s1[i-1] + a
            b = s2[j-1] + b
            i -= 1
            j -= 1
        elif state == 'gap_a':
            a = '-' + a
            b = s2[j-1] + b
            j -= 1
        elif state == 'gap_b':
            a = s1[i-1] + a
            b = '-' + b
            i -= 1
        else:
            break
        if i > 0 or j > 0:
            state = ptr[i][j]
    return {'a': a, 'b': b, 'score': match_scores[end_i][end_j] if end_state == 'diag' else gap_a_scores[end_i][end_j] if end_state == 'gap_a' else gap_b_scores[end_i][end_j], 'start1': i, 'end1': end_i, 'start2': j, 'end2': end_j, 'mode': mode}
