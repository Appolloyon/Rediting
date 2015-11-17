from unittest import TestCase
from unittest import skip

from rediting.util import sequence

class TestSeqs(TestCase):

    def test_compare_seqs_nonequal(self):
        seq1 = 'abcdefghi'
        seq2 = 'jklmnopqr'
        # 7 of 9 resiues used as standard setting in main program
        self.assertFalse(sequence.compare_seqs(seq1,seq2,7))

    def test_compare_seqs_allequal(self):
        seq1 = 'abcdefghi'
        seq2 = 'abcdefghi'
        self.assertTrue(sequence.compare_seqs(seq1,seq2,7))

    def test_compare_seqs_below_threshold(self):
        seq1 = 'abcdefghi'
        seq2 = 'abcdefxyz'
        # only six characters are the same above
        # therefore the test should fail
        self.assertFalse(sequence.compare_seqs(seq1,seq2,7))

    def test_compare_seqs_opening_gap1(self):
        seq1 = '-bcdefghi'
        seq2 = 'abcdefghi'
        # first character is a gap, therefore should evaluate False
        self.assertFalse(sequence.compare_seqs(seq1,seq2,7))

    def test_compare_seqs_opening_gap2(self):
        seq1 = '-bcdefghi'
        seq2 = '-bcdefghi'
        self.assertFalse(sequence.compare_seqs(seq1,seq2,7))

    def test_compare_seqs_internal_gap1(self):
        seq1 = 'a-cdefghi'
        seq2 = 'abcdefghi'
        self.assertTrue(sequence.compare_seqs(seq1,seq2,7))

    def test_compare_seqs_internal_gap2(self):
        seq1 = 'a-cdefghi'
        seq2 = 'abc-efghi'
        self.assertTrue(sequence.compare_seqs(seq1,seq2,7))

    # GC function returns floats, therefore must use decimals
    def test_GC_100(self):
        seq = 'GCGCGCGCGCGCGCGC'
        self.assertEqual(sequence.calc_gc(seq), 100.0)

    def test_GC_50(self):
        seq = 'AGTCAGTCAGTC'
        self.assertEqual(sequence.calc_gc(seq), 50.0)

    def test_GC_0(self):
        seq = 'ATATATATATATATATATATATATA'
        self.assertEqual(sequence.calc_gc(seq), 0.0)

    # calc_percent takes window size, use 20 to test
    def test_calc_percent_100(self):
        seq = '11111111111111111111'
        result = sequence.calc_percent(seq,0,20,20)
        self.assertEqual(result, 100.0)

    def test_calc_percent_50(self):
        seq = '11111111110000000000'
        result = sequence.calc_percent(seq,0,20,20)
        self.assertEqual(result, 50.0)

    def test_calc_percent_0(self):
        seq = '00000000000000000000'
        result = sequence.calc_percent(seq,0,20,20)
        self.assertEqual(result, 0.0)

    # get_indices takes window size, here use only 2 for test purposes
    def test_get_indices(self):
        seq = 'abcdefg'
        indices = sequence.get_indices(seq,2)
        self.assertEqual(indices, [[0,2],[1,3],[2,4],[3,5],[4,6],[5,7]])

    # again, use a value of 2 for the 'size' variable
    def test_get_non_overlapping_indices(self):
        seq = 'abcdefg'
        indices = sequence.get_non_overlapping_indices(seq,0,2,indices=[])
        self.assertEqual(indices, [[0,2],[2,4],[4,6]])

    @skip("skipped for now")
    def test_polyT(self):
        pass

    @skip("skipped for now")
    def test_polyTpercent(self):
        pass

    def test_codon_position(self):
        results = []
        for i in range(1,5):
            results.append(sequence.incr_codon_position(i))
        self.assertEqual(results, [2,3,1,1])

    def test_calculate_codons1(self):
        seq = 'ATGCTTATTCGG'
        self.assertEqual(sequence.calculate_codons(seq,1), ['ATG','CTT','ATT','CGG'])

    def test_calculate_codons2(self):
        seq = 'ATGCTTATTCGG'
        self.assertEqual(sequence.calculate_codons(seq,2), ['GCT','TAT','TCG'])

    def test_calculate_codons3(self):
        seq = 'ATGCTTATTCGG'
        self.assertEqual(sequence.calculate_codons(seq,3), ['TGC','TTA','TTC'])
