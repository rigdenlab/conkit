"""Testing facility for conkit.core.sequence"""

__author__ = "Felix Simkovic"
__date__ = "07 Sep 2016"

from conkit.core.Sequence import Sequence

import unittest


class Test(unittest.TestCase):

    def test_remark_1(self):
        sequence = Sequence('foo', 'GSMFTPK')
        sequence.remark = 'bar'
        self.assertEqual(['bar'], sequence.remark)

    def test_remark_2(self):
        sequence = Sequence('foo', 'GSMFTPK')
        sequence.remark = 'bar'
        sequence.remark = 'baz'
        self.assertEqual(['bar', 'baz'], sequence.remark)

    def test_seq_1(self):
        sequence = Sequence('foo', 'GSMFTPK')
        self.assertEqual('foo', sequence.id)
        self.assertEqual('GSMFTPK', sequence.seq)

    def test_seq_2(self):
        sequence = Sequence('foo', 'GSMFTPK')
        sequence.seq = 'AAAAAA'
        self.assertEqual('foo', sequence.id)
        self.assertEqual('AAAAAA', sequence.seq)

    def test_seq_3(self):
        sequence = Sequence('foo', 'GSMFTPK')
        with self.assertRaises(ValueError):
            sequence.seq = 'A2A'

    def test_seq_4(self):
        sequence = Sequence('foo', 'GSMFTPK')
        sequence.seq = '-------'

    def test_seq_len_1(self):
        sequence = Sequence('foo', 'GSMFTPK')
        self.assertEqual('foo', sequence.id)
        self.assertEqual('GSMFTPK', sequence.seq)
        self.assertEqual(7, sequence.seq_len)

    def test_seq_len_2(self):
        sequence = Sequence('foo', 'GSMFTPK')
        self.assertEqual('foo', sequence.id)
        self.assertEqual('GSMFTPK', sequence.seq)
        self.assertEqual(7, sequence.seq_len)
        sequence.seq = 'AAAAAAAAAA'
        self.assertEqual('foo', sequence.id)
        self.assertEqual('AAAAAAAAAA', sequence.seq)
        self.assertEqual(10, sequence.seq_len)

    def test_align_local_1(self):
        sequence1 = Sequence('foo', 'GSMFTPKPPQDSAVIKAGYCVKQGAVMKNWKRRYFQLDENTI'
                                    'GYFKSELEKEPLRVIPLKEVHKVQECKQSDIMMRDNLFEIVT'
                                    'TSRTFYVQADSPEEMHSWIKAVSGAIVAQRGPGRSASSEHP')
        sequence2 = Sequence('bar', 'Q-------YF-------P------------------------'
                                    '--F----------VQADSPEEMHSWIKAVSGAIVAQR')
        sequence1.align_local(sequence2, id_chars=2, nonid_chars=1, gap_open_pen=-0.5, gap_ext_pen=-0.2, inplace=True)
        aligned1 = "GSMFTPKPPQDSAVIKAGYCVKQGAVMKNWKRRYFQLDENTIGYFKSELEKEPLRVIPLKEVHKVQECKQSDIMM" \
                   "RDNLFEIVTTSRTFYVQADSPEEMHSWIKAVSGAIVAQRGPGRSASSEHP"
        aligned2 = "-----------------------------------Q-------YF-------P----------------------" \
                   "----F----------VQADSPEEMHSWIKAVSGAIVAQR-----------"
        self.assertEqual(aligned1, sequence1.seq)
        self.assertEqual(aligned2, sequence2.seq)

    def test_align_local_2(self):
        sequence1 = Sequence('foo', 'DDLTISSLAKGETTKAAFNQMVQGHKLPAWVMKGGTYTPAQTV'
                                    'TLGDETYQVMSACKPHDCGSQRIAVMWSEKSNQMTGLFSTIDE'
                                    'KTSQEKLTWLNVNDALSIDGKTVLFAALTGSLENHPDGFNFKV'
                                    'FGRCELAAAMKRHGLDNYRGYSLGNWVCAAKFESNFNTQATNR'
                                    'NTDGSTDYGILQINSRWWCNDGRTPGSRNLCNIPCSALLSSDI'
                                    'TASVNCAKKIVSDGNGMNAWVAWRNRCKGTDVQAWIRGCR')
        sequence2 = Sequence('bar', '-------------------------------------------'
                                    '-------------------------------------------'
                                    '--------W------------TV--------------------'
                                    'F--C----AM---GLD-----------C--KFE-NF-------'
                                    'N-D-----G---------C-D----G--NLC-IP--------I'
                                    '--------------NG--------------D----IRGC-')
        sequence1.align_local(sequence2, id_chars=2, nonid_chars=1, gap_open_pen=-0.5, gap_ext_pen=-0.2, inplace=True)
        aligned1 = "DDLTISSLAKGETTKAAFNQMVQGHKLPAWVMKGGTYTPAQTVTLGDETYQVMSACKPHDCGSQRIAVMWSEKSN" \
                   "QMTGLFSTIDEKTSQEKLTWLNVNDALSIDGKTVLFAALTGSLENHPDGFNFKVFGRCELAAAMKRHGLDNYRGY" \
                   "SLGNWVCAAKFESNFNTQATNRNTDGSTDYGILQINSRWWCNDGRTPGSRNLCNIPCSALLSSDITASVNCAKKI" \
                   "VSDGNGMNAWVAWRNRCKGTDVQAWIRGCR"
        aligned2 = "---------------------------------------------------------------------------" \
                   "-------------------W------------TV--------------------F--C----AM---GLD-----" \
                   "------C--KFE-NF-------N-D-----G---------C-D----G--NLC-IP--------I----------" \
                   "----NG--------------D----IRGC-"
        self.assertEqual(aligned1, sequence1.seq)
        self.assertEqual(aligned2, sequence2.seq)

    def test_align_local_3(self):
        sequence1 = Sequence('foo', '------------------------------------------'
                                    '------------------------------------------'
                                    '----------W------------TV-----------------'
                                    '---F--C----AM---GLD-----------C--KFE-NF---'
                                    '----N-D-----G---------C-D----G--NLC-IP----'
                                    '----I--------------NG--------------D----IR'
                                    'GC-')
        sequence2 = Sequence('bar', '-D-------------------------------GGTYTP---'
                                    '-------------C-PHDCGS-R-------------------'
                                    '------------------------------TG--EN------'
                                    '-KV------------------------------KFESN-N-Q'
                                    'ATNR------D----Q--------------------------'
                                    '------------------------WVA--NR-----------'
                                    '---')
        sequence1.align_local(sequence2, id_chars=2, nonid_chars=1, gap_open_pen=-1.0, gap_ext_pen=-0.5, inplace=True)
        aligned1 = "---------------------------------------------------------------------------" \
                   "-------------------W------------TV--------------------F--C----AM---GLD-----" \
                   "------C--KFE-NF-------N-D-----G---------C-D----G--NLC-IP--------I----------" \
                   "----NG--------------D----IRGC-"
        aligned2 = "-D-------------------------------GGTYTP----------------C-PHDCGS-R----------" \
                   "---------------------------------------TG--EN-------KV---------------------" \
                   "---------KFESN-N-QATNR------D----Q-----------------------------------------" \
                   "---------WVA--NR--------------"
        self.assertEqual(aligned1, sequence1.seq)
        self.assertEqual(aligned2, sequence2.seq)


if __name__ == "__main__":
    unittest.main(verbosity=2)
