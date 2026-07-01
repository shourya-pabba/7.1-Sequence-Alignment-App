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
