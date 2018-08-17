"""Testing facility for conkit.io.ncont"""

__author__ = "Felix Simkovic"
__date__ = "26 Oct 2016"

import os
import unittest

from conkit.io.ncont import NcontParser
from conkit.io._iotools import create_tmp_f


class TestNcontParser(unittest.TestCase):
    def test_read_1(self):
        content = """
 ###############################################################
 ###############################################################
 ###############################################################
 ### CCP4 7.0.045: NCONT            version 7.0.045 :         ##
 ###############################################################
 User: unknown  Run date: 16/11/2017 Run time: 09:51:16


 Please reference: Collaborative Computational Project, Number 4. 2011.
 "Overview of the CCP4 suite and current developments". Acta Cryst. D67, 235-242.
 as well as any specific reference in the program write-up.


 ------------------------------------------------------------------------------

 PDB file /home/felix/Dropbox/subprojects/ample-rio/tmppyndBZ.pdb has been read in.

 ------------------------------------------------------------------------------
  Input cards

 Data line--- source A//CA
 Data line--- target B//CA
 Data line--- mindist 0.0
 Data line--- maxdist 1.5
 Data line--- sort target inc

 ------------------------------------------------------------------------------

 Selected  15  source atoms
 Selected  8154  target atoms

 ------------------------------------------------------------------------------

  12 contacts found:

      SOURCE ATOMS               TARGET ATOMS         DISTANCE

 /1/A/ 127(ALA). / CA [ C]:  /1/B/ 129(GLY). / CA [ C]:   0.41
 /1/A/ 128(GLN). / CA [ C]:  /1/B/ 130(ALA). / CA [ C]:   0.37
 /1/A/ 129(GLY). / CA [ C]:  /1/B/ 131(MET). / CA [ C]:   0.19
 /1/A/ 130(ALA). / CA [ C]:  /1/B/ 132(ASN). / CA [ C]:   0.10
 /1/A/ 131(MET). / CA [ C]:  /1/B/ 133(LYS). / CA [ C]:   0.35
 /1/A/ 132(ASN). / CA [ C]:  /1/B/ 134(ALA). / CA [ C]:   0.20
 /1/A/ 133(LYS). / CA [ C]:  /1/B/ 135(LEU). / CA [ C]:   0.39
 /1/A/ 134(ALA). / CA [ C]:  /1/B/ 136(GLU). / CA [ C]:   0.34
 /1/A/ 135(LEU). / CA [ C]:  /1/B/ 137(LEU). / CA [ C]:   0.48
 /1/A/ 136(GLU). / CA [ C]:  /1/B/ 138(PHE). / CA [ C]:   0.46
 /1/A/ 137(LEU). / CA [ C]:  /1/B/ 139(ARG). / CA [ C]:   0.47
 /1/A/ 138(PHE). / CA [ C]:  /1/B/ 140(LYS). / CA [ C]:   0.89


  Total 12 contacts

--------------------------------------------------------------------------

 NCONT:  Normal termination
Times: User:       0.1s System:    0.0s Elapsed:     0:00
"""
        f_name = create_tmp_f(content=content)
        with open(f_name, 'r') as f_in:
            contact_file = NcontParser().read(f_in)
        os.unlink(f_name)
        cmap = contact_file.top
        self.assertEqual(cmap.ncontacts, 12)
        self.assertEqual([c.res1_seq for c in cmap], list(range(127, 139)))
        self.assertEqual([c.res2_seq for c in cmap], list(range(129, 141)))
        self.assertEqual(set([c.res1_chain for c in cmap]), set(["A"]))
        self.assertEqual(set([c.res2_chain for c in cmap]), set(["B"]))

    def test_read_2(self):
        content = """
 ###############################################################
 ###############################################################
 ###############################################################
 ### CCP4 7.0.045: NCONT            version 7.0.045 :         ##
 ###############################################################
 User: unknown  Run date: 16/11/2017 Run time: 10:04:06


 Please reference: Collaborative Computational Project, Number 4. 2011.
 "Overview of the CCP4 suite and current developments". Acta Cryst. D67, 235-242.
 as well as any specific reference in the program write-up.


 ------------------------------------------------------------------------------

 PDB file /home/felix/Dropbox/subprojects/ample-rio/tmp5x3eJ1.pdb has been read in.

 ------------------------------------------------------------------------------
  Input cards

 Data line--- source A//CA
 Data line--- target B//CA
 Data line--- mindist 0.0
 Data line--- maxdist 0.01
 Data line--- sort target inc

 ------------------------------------------------------------------------------

 Selected  15  source atoms
 Selected  8154  target atoms

 ------------------------------------------------------------------------------

  NO CONTACTS FOUND.



  Total 0 contacts

--------------------------------------------------------------------------

 NCONT:  Normal termination
Times: User:       0.1s System:    0.0s Elapsed:     0:00
"""
        f_name = create_tmp_f(content=content)
        with open(f_name, 'r') as f_in:
            contact_file = NcontParser().read(f_in)
        os.unlink(f_name)
        cmap = contact_file.top
        self.assertEqual(cmap.ncontacts, 0)

    def test_read_3(self):
        content = """
 ###############################################################
 ###############################################################
 ###############################################################
 ### CCP4 7.0.045: NCONT            version 7.0.045 :         ##
 ###############################################################
 User: unknown  Run date: 16/11/2017 Run time: 10:07:41


 Please reference: Collaborative Computational Project, Number 4. 2011.
 "Overview of the CCP4 suite and current developments". Acta Cryst. D67, 235-242.
 as well as any specific reference in the program write-up.


 ------------------------------------------------------------------------------

 PDB file /home/felix/Dropbox/subprojects/ample-rio/tmpftYndd.pdb has been read in.

 ------------------------------------------------------------------------------
  Input cards

 Data line--- source A//*
 Data line--- target B//*
 Data line--- mindist 0.0
 Data line--- maxdist 0.2
 Data line--- sort target inc

 ------------------------------------------------------------------------------

 Selected  74  source atoms
 Selected  57888  target atoms

 ------------------------------------------------------------------------------

  9 contacts found:

      SOURCE ATOMS               TARGET ATOMS         DISTANCE

 /1/A/ 129(GLY). / CA [ C]:  /1/B/ 131(MET). / CA [ C]:   0.19
 /1/A/ 129(GLY). / C  [ C]:  /1/B/ 131(MET). / C  [ C]:   0.13
 /1/A/ 129(GLY). / O  [ O]:  /1/B/ 131(MET). / O  [ O]:   0.15
 /1/A/ 130(ALA). / N  [ N]:  /1/B/ 132(ASN). / N  [ N]:   0.10
 /1/A/ 130(ALA). / CA [ C]:  /1/B/ 132(ASN). / CA [ C]:   0.10
 /1/A/ 130(ALA). / CB [ C]:  /1/B/ 132(ASN). / CB [ C]:   0.14
 /1/A/ 131(MET). / O  [ O]:  /1/B/ 133(LYS). / O  [ O]:   0.19
 /1/A/ 132(ASN). / C  [ C]:  /1/B/ 134(ALA). / C  [ C]:   0.18
 /1/A/ 135(LEU). / N  [ N]:  /1/B/ 137(LEU). / N  [ N]:   0.20


  Total 9 contacts

--------------------------------------------------------------------------

 NCONT:  Normal termination
Times: User:       0.1s System:    0.0s Elapsed:     0:00
"""
        f_name = create_tmp_f(content=content)
        with open(f_name, 'r') as f_in:
            contact_file = NcontParser().read(f_in)
        os.unlink(f_name)
        cmap = contact_file.top
        self.assertEqual(cmap.ncontacts, 5)
        self.assertEqual([c.res1_seq for c in cmap], [129, 130, 131, 132, 135])
        self.assertEqual([c.res2_seq for c in cmap], [131, 132, 133, 134, 137])
        self.assertEqual(set([c.res1_chain for c in cmap]), set(["A"]))
        self.assertEqual(set([c.res2_chain for c in cmap]), set(["B"]))


if __name__ == "__main__":
    unittest.main(verbosity=2)
