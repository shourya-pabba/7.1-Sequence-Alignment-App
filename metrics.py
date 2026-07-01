from scoring import conservative


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
