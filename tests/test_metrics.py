import unittest

from align_algorithms import global_alignment
from metrics import metrics


class MetricsTests(unittest.TestCase):
    def test_dna_sequences_count_all_mismatches(self):
        aln = global_alignment('CTCTCT', 'CTCAGA', match=2, mismatch=3, gap=4, matrix='simple')
        result = metrics(aln, sequence_type='dna')

        self.assertEqual(result['matches'], 3)
        self.assertEqual(result['mismatches'], 3)
        self.assertEqual(result['conservative'], 0)

    def test_protein_sequences_can_still_track_conservative_substitutions(self):
        aln = global_alignment('ACDE', 'AGDE', match=2, mismatch=3, gap=4, matrix='simple')
        result = metrics(aln, sequence_type='protein')

        self.assertEqual(result['matches'], 3)
        self.assertEqual(result['mismatches'], 0)
        self.assertEqual(result['conservative'], 1)


if __name__ == '__main__':
    unittest.main()
