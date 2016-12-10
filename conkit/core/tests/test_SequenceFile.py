"""Testing facility for conkit.core.sequencesequence_file"""

__author__ = "Felix Simkovic"
__date__ = "07 Sep 2016"

from conkit.core.Sequence import Sequence
from conkit.core.SequenceFile import SequenceFile

import unittest


class Test(unittest.TestCase):

    def test_is_alignment(self):
        # ======================================================
        # Test Case 1
        sequence_file = SequenceFile('test')
        sequence_file.add(Sequence('foo', 'AAAAA'))
        sequence_file.add(Sequence('bar', 'BBBBB'))
        self.assertTrue(sequence_file.is_alignment)
        # ======================================================
        # Test Case 2
        sequence_file = SequenceFile('test')
        sequence_file.add(Sequence('foo', 'AAAAA'))
        sequence_file.add(Sequence('bar', 'BBBB'))
        self.assertFalse(sequence_file.is_alignment)

    def test_nseqs(self):
        # ======================================================
        # Test Case 1
        sequence_file = SequenceFile('test')
        self.assertEqual(0, sequence_file.nseqs)
        # ======================================================
        # Test Case 2
        sequence_file = SequenceFile('test')
        sequence_file.add(Sequence('foo', 'AAAAA'))
        self.assertEqual(1, sequence_file.nseqs)
        # ======================================================
        # Test Case 3
        sequence_file = SequenceFile('test')
        sequence_file.add(Sequence('foo', 'AAAAA'))
        sequence_file.add(Sequence('bar', 'BBBBB'))
        self.assertEqual(2, sequence_file.nseqs)

    def test_remark(self):
        # ======================================================
        # Test Case 1
        sequence_file = SequenceFile('test')
        sequence_file.remark = 'Hello'
        self.assertEqual(['Hello'], sequence_file.remark)
        # ======================================================
        # Test Case 2
        sequence_file = SequenceFile('test')
        sequence_file.remark = 'Hello'
        sequence_file.remark = 'World'
        self.assertEqual(['Hello', 'World'], sequence_file.remark)
        # ======================================================
        # Test Case 3
        sequence_file = SequenceFile('test')
        sequence_file.remark = 'Hello'
        sequence_file.remark = '5'
        sequence_file.remark = 'World'
        sequence_file.remark = '!'
        self.assertEqual(['Hello', '5', 'World', '!'], sequence_file.remark)
        # ======================================================
        # Test Case 4
        sequence_file = SequenceFile('test')
        self.assertEqual([], sequence_file.remark)
        # ======================================================
        # Test Case 5
        sequence_file = SequenceFile('test')
        sequence_file.remark = 'hello'
        sequence = Sequence('foo', 'GSMFTPK')
        sequence.remark = 'bar'
        sequence_file.add(sequence)
        self.assertEqual(['hello'], sequence_file.remark)
        self.assertEqual(['bar'], sequence_file[0].remark)

    def test_top_sequence(self):
        # ======================================================
        # Test Case 1
        sequence_file = SequenceFile('test')
        self.assertEqual(None, sequence_file.top_sequence)
        # ======================================================
        # Test Case 2
        sequence_file = SequenceFile('test')
        sequence1 = Sequence('foo', 'AAAAA')
        sequence_file.add(sequence1)
        self.assertEqual(sequence1, sequence_file.top_sequence)
        # ======================================================
        # Test Case 3
        sequence_file = SequenceFile('test')
        sequence1 = Sequence('foo', 'AAAAA')
        sequence2 = Sequence('bar', 'BBBBB')
        sequence_file.add(sequence1)
        sequence_file.add(sequence2)
        self.assertEqual(sequence1, sequence_file.top_sequence)

    def test_calculate_meff(self):
        # ======================================================
        # Test Case 1
        sequence_file = SequenceFile('test')
        for s in [Sequence('foo', 'AAAAAAA'), Sequence('bar', 'AAAAAAA'),
                  Sequence('cho', 'AAAAAAA'), Sequence('baz', 'AAAAAAA')]:
            sequence_file.add(s)
        m_eff = sequence_file.calculate_meff(identity=0.7)
        self.assertTrue(isinstance(m_eff, int))
        self.assertEqual(1, m_eff)
        # ======================================================
        # Test Case 2
        sequence_file = SequenceFile('test')
        for s in [Sequence('foo', 'AAAAAAA'), Sequence('bar', 'AAAAAAA'),
                  Sequence('cho', 'AAAAAAA'), Sequence('baz', 'BBBBBBB')]:
            sequence_file.add(s)
        m_eff = sequence_file.calculate_meff(identity=0.7)
        self.assertTrue(isinstance(m_eff, int))
        self.assertEqual(2, m_eff)
        # ======================================================
        # Test Case 3
        sequence_file = SequenceFile('test')
        for s in [Sequence('foo', 'AAAAAAA'), Sequence('bar', 'A-AABA-'),
                  Sequence('cho', 'B-BAA--'), Sequence('baz', 'BBBBBBB')]:
            sequence_file.add(s)
        m_eff = sequence_file.calculate_meff(identity=0.7)
        self.assertTrue(isinstance(m_eff, int))
        self.assertEqual(4, m_eff)
        # ======================================================
        # Test Case 4
        sequence_file = SequenceFile('test')
        for s in [Sequence('foo', 'AAAAAAA'), Sequence('bar', 'AAAABA-'),
                  Sequence('cho', 'B-BAA--'), Sequence('baz', 'BBBBBBB')]:
            sequence_file.add(s)
        m_eff = sequence_file.calculate_meff(identity=0.7)
        self.assertTrue(isinstance(m_eff, int))
        self.assertEqual(3, m_eff)
        # ======================================================
        # Test Case 5
        sequence_file = SequenceFile('test')
        for s in [Sequence('foo', 'AAAAAAA'), Sequence('bar', 'AA-ABA-'),
                  Sequence('cho', 'B-BAA--'), Sequence('baz', 'BBBBBBB')]:
            sequence_file.add(s)
        m_eff = sequence_file.calculate_meff(identity=0.6)
        self.assertTrue(isinstance(m_eff, int))
        self.assertEqual(4, m_eff)
        self.assertNotEqual(3, m_eff)
        self.assertNotEqual(3, m_eff)
        # ======================================================
        # Test Case 5
        sequence_file = SequenceFile('test')
        for s in [Sequence('foo', 'AAAAAAA'), Sequence('bar', 'AA-ABA-'),
                  Sequence('cho', 'AAACBAA'), Sequence('doo', 'B-BAA--'),
                  Sequence('miu', 'BBBBBBB'), Sequence('nop', 'AAAAAAB')]:
            sequence_file.add(s)
        m_eff = sequence_file.calculate_meff(identity=0.6)
        self.assertTrue(isinstance(m_eff, int))
        self.assertEqual(4, m_eff)

    def test_calculate_freq(self):
        # ======================================================
        # Test Case 1
        sequence_file = SequenceFile('test')
        for s in [Sequence('foo', 'AAAAAAA'), Sequence('bar', 'A-AAAA-'), Sequence('cho', '--AAA--')]:
            sequence_file.add(s)
        calculated_freqs = [round(i, 6) for i in sequence_file.calculate_freq()]
        self.assertEqual([0.666667, 0.333333, 1.0, 1.0, 1.0, 0.666667, 0.333333], calculated_freqs)
        # ======================================================
        # Test Case 2
        sequence_file = SequenceFile('test')
        for s in [Sequence('foo', '-------'), Sequence('bar', '-------'), Sequence('cho', '-------')]:
            sequence_file.add(s)
        calculated_freqs = [round(i, 6) for i in sequence_file.calculate_freq()]
        self.assertEqual([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], calculated_freqs)
        # ======================================================
        # Test Case 3
        sequence_file = SequenceFile('test')
        for s in [Sequence('foo', 'AAAAAAA'), Sequence('bar', 'AAAAAAA'), Sequence('cho', 'AAAAAAA')]:
            sequence_file.add(s)
        calculated_freqs = [round(i, 6) for i in sequence_file.calculate_freq()]
        self.assertEqual([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0], calculated_freqs)

    def test_sort(self):
        # ======================================================
        # Test Case 1
        sequence_file = SequenceFile('test')
        for seq in [Sequence('foo', 'AAAAA'), Sequence('bar', 'BBBBB'), Sequence('doe', 'CCCCC')]:
            sequence_file.add(seq)
        sequence_file_sorted = sequence_file.sort('id', reverse=False, inplace=False)
        self.assertEqual(['bar', 'doe', 'foo'], [s.id for s in sequence_file_sorted])
        self.assertEqual(['BBBBB', 'CCCCC', 'AAAAA'], [s.seq for s in sequence_file_sorted])
        self.assertNotEqual(sequence_file, sequence_file_sorted)
        # ======================================================
        # Test Case 2
        sequence_file = SequenceFile('test')
        for seq in [Sequence('foo', 'AAAAA'), Sequence('bar', 'BBBBB'), Sequence('doe', 'CCCCC')]:
            sequence_file.add(seq)
        sequence_file_sorted = sequence_file.sort('id', reverse=True, inplace=False)
        self.assertEqual(['foo', 'doe', 'bar'], [s.id for s in sequence_file_sorted])
        self.assertEqual(['AAAAA', 'CCCCC', 'BBBBB'], [s.seq for s in sequence_file_sorted])
        self.assertNotEqual(sequence_file, sequence_file_sorted)
        # ======================================================
        # Test Case 3
        sequence_file = SequenceFile('test')
        for seq in [Sequence('foo', 'AAAAA'), Sequence('bar', 'BBBBB'), Sequence('doe', 'CCCCC')]:
            sequence_file.add(seq)
        sequence_file_sorted = sequence_file.sort('seq', reverse=False, inplace=True)
        self.assertEqual(['foo', 'bar', 'doe'], [s.id for s in sequence_file_sorted])
        self.assertEqual(['AAAAA', 'BBBBB', 'CCCCC'], [s.seq for s in sequence_file_sorted])
        self.assertEqual(sequence_file, sequence_file_sorted)
        # ======================================================
        # Test Case 4
        sequence_file = SequenceFile('test')
        for seq in [Sequence('foo', 'AAAAA'), Sequence('bar', 'BBBBB'), Sequence('doe', 'CCCCC')]:
            sequence_file.add(seq)
        sequence_file_sorted = sequence_file.sort('seq', reverse=True, inplace=True)
        self.assertEqual(['doe', 'bar', 'foo'], [s.id for s in sequence_file_sorted])
        self.assertEqual(['CCCCC', 'BBBBB', 'AAAAA'], [s.seq for s in sequence_file_sorted])
        self.assertEqual(sequence_file, sequence_file_sorted)


if __name__ == "__main__":
    unittest.main(verbosity=2)
