import os
import ssl
import urllib.request

from flask import Flask, request, render_template
from align_algorithms import global_alignment, local_alignment
from metrics import metrics
from diagnostics import diagnose

app = Flask(__name__, template_folder='templates', static_folder='static')

DATASET_LIBRARY = {
    'custom': {'label': 'Custom input', 'sequence': None, 'sequence_type': 'dna', 'scoring_profile': 'balanced'},
    'sars_cov_spike_dna': {'label': 'SARS-CoV spike DNA', 'db': 'nuccore', 'accession': 'AY278741.1', 'start': 21492, 'end': 25259, 'sequence_type': 'dna', 'scoring_profile': 'balanced'},
    'sars_cov2_spike_dna': {'label': 'SARS-CoV-2 spike DNA', 'db': 'nuccore', 'accession': 'NC_045512.2', 'start': 21563, 'end': 25384, 'sequence_type': 'dna', 'scoring_profile': 'balanced'},
    'sars_cov_spike_protein': {'label': 'SARS-CoV spike protein', 'db': 'protein', 'accession': 'NP_828851.1', 'sequence_type': 'protein', 'scoring_profile': 'blosum62'},
    'sars_cov2_spike_protein': {'label': 'SARS-CoV-2 spike protein', 'db': 'protein', 'accession': 'YP_009724390.1', 'sequence_type': 'protein', 'scoring_profile': 'blosum62'},
}


def fetch_ncbi_sequence(dataset_entry):
    if dataset_entry.get('sequence') is not None:
        return dataset_entry['sequence']

    db = dataset_entry['db']
    accession = dataset_entry['accession']
    url = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db={db}&id={accession}&rettype=fasta&retmode=text'
    if dataset_entry.get('start') is not None and dataset_entry.get('end') is not None:
        url += f'&from={dataset_entry["start"]}&to={dataset_entry["end"]}'

    request_obj = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
    context = ssl._create_unverified_context()
    with urllib.request.urlopen(request_obj, context=context, timeout=30) as response:
        text = response.read().decode('utf-8', 'ignore')

    lines = [line.strip() for line in text.splitlines() if line.strip()]
    seq = ''.join(line for line in lines if not line.startswith('>'))
    return ''.join(seq.split()).upper()


def dataset_context(preset_name):
    entry = DATASET_LIBRARY.get(preset_name, DATASET_LIBRARY['custom'])
    if entry.get('sequence') is not None:
        return {
            'label': entry['label'],
            'sequence': entry['sequence'],
            'sequence_type': entry['sequence_type'],
            'scoring_profile': entry['scoring_profile'],
        }

    sequence = fetch_ncbi_sequence(entry)
    return {
        'label': entry['label'],
        'sequence': sequence,
        'sequence_type': entry['sequence_type'],
        'scoring_profile': entry['scoring_profile'],
    }


def resolve_scoring(sequence_type, scoring_profile):
    if sequence_type == 'protein':
        if scoring_profile == 'blosum62':
            return {'match': 2, 'mismatch': 3, 'gap': 4, 'matrix': 'blosum62', 'name': 'BLOSUM62 protein profile'}
        if scoring_profile == 'blosum_like':
            return {'match': 2, 'mismatch': 3, 'gap': 4, 'matrix': 'blosum_like', 'name': 'Blosum-like protein profile'}
        if scoring_profile == 'pam_like':
            return {'match': 2, 'mismatch': 3, 'gap': 4, 'matrix': 'pam_like', 'name': 'Pam-like protein profile'}
        return {'match': 2, 'mismatch': 3, 'gap': 4, 'matrix': 'simple', 'name': 'Simple protein profile'}

    if scoring_profile == 'strict':
        return {'match': 3, 'mismatch': 4, 'gap': 2, 'matrix': 'simple', 'name': 'Strict DNA profile'}
    if scoring_profile == 'permissive':
        return {'match': 1, 'mismatch': 1, 'gap': 3, 'matrix': 'simple', 'name': 'Permissive DNA profile'}
    return {'match': 2, 'mismatch': 3, 'gap': 4, 'matrix': 'simple', 'name': 'Balanced DNA profile'}


@app.route('/', methods=['GET', 'POST'])
def index():
    result = None
    seq1_value = ''
    seq2_value = ''
    seq1_preset = request.form.get('seq1_dataset', 'custom')
    seq2_preset = request.form.get('seq2_dataset', 'custom')

    if request.method == 'POST':
        seq1_value = request.form.get('seq1', '').strip()
        seq2_value = request.form.get('seq2', '').strip()

        seq1_context = dataset_context(seq1_preset) if seq1_preset != 'custom' else None
        seq2_context = dataset_context(seq2_preset) if seq2_preset != 'custom' else None

        if seq1_context:
            seq1_value = seq1_context['sequence']
        if seq2_context:
            seq2_value = seq2_context['sequence']

        sequence_type = request.form.get('sequence_type', 'dna')
        scoring_profile = request.form.get('scoring_profile', 'balanced')
        if seq1_context or seq2_context:
            if (seq1_context and seq1_context['sequence_type'] == 'protein') or (seq2_context and seq2_context['sequence_type'] == 'protein'):
                sequence_type = 'protein'
                scoring_profile = 'blosum62'
            else:
                sequence_type = 'dna'
                scoring_profile = 'balanced'

        mode = request.form.get('mode', 'global')
        scoring = resolve_scoring(sequence_type, scoring_profile)
        match = scoring['match']
        mismatch = scoring['mismatch']
        gap = scoring['gap']
        matrix = scoring['matrix']

        if seq1_value and seq2_value:
            aligner = global_alignment if mode == 'global' else local_alignment
            alignment = aligner(seq1_value, seq2_value, match=match, mismatch=mismatch, gap=gap, matrix=matrix)
            metrics_result = metrics(alignment, sequence_type=sequence_type)
            labels, warnings, story = diagnose(metrics_result, mode, mismatch, gap)
            result = {
                'mode': mode,
                'score': alignment['score'],
                'alignment_a': alignment['a'],
                'alignment_b': alignment['b'],
                'metrics': metrics_result,
                'labels': labels,
                'warnings': warnings,
                'story': story,
                'scoring': scoring,
                'seq1_label': seq1_context['label'] if seq1_context else 'Custom input',
                'seq2_label': seq2_context['label'] if seq2_context else 'Custom input',
            }

    return render_template('index.html', result=result, seq1_value=seq1_value, seq2_value=seq2_value, seq1_preset=seq1_preset, seq2_preset=seq2_preset)


if __name__ == '__main__':
    app.run(host='0.0.0.0', port=int(os.environ.get('PORT', '5000')), debug=True)
