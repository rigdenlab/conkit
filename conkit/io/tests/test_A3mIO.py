"""Testing facility for conkit.io.A3mIO"""

__author__ = "Felix Simkovic"
__date__ = "11 Sep 2016"

from conkit.io.A3mIO import A3mIO
from conkit.io._iotools import create_tmp_f

import os
import unittest


class Test(unittest.TestCase):

    def test_read_1(self):
        msa = """>d1a1x__ b.63.1.1 (-) p13-MTCP1 {Human (Homo sapiens)}
PPDHLWVHQEGIYRDEYQRTWVAVVEEETSFLRARVQQIQVPLGDAARPSHLLTSQL
>gi|6678257|ref|NP_033363.1|:(7-103) T-cell lymphoma breakpoint 1 [Mus musculus]
HPNRLWIWEKHVYLDEFRRSWLPVVIKSNEKFQVILRQEDVTLGEAMSPSQLVPYEL
>gi|7305557|ref|NP_038800.1|:(8-103) T-cell leukemia/lymphoma 1B, 3 [Mus musculus]
PPRFLVCTRDDIYEDENGRQWVVAKVETSRSpygsrietcITVHLQHMTTIPQEPTPQQPINNNSL
>gi|11415028|ref|NP_068801.1|:(2-106) T-cell lymphoma-1; T-cell lymphoma-1A [Homo sapiens]
HPDRLWAWEKFVYLDEKQHAWLPLTIEikDRLQLRVLLRREDVVLGRPMTPTQIGPSLL
>gi|7305561|ref|NP_038804.1|:(7-103) T-cell leukemia/lymphoma 1B, 5 [Mus musculus]
----------GIYEDEHHRVWIAVNVETSHSSHgnrietcvtVHLQHMTTLPQEPTPQQPINNNSL
>gi|7305553|ref|NP_038801.1|:(5-103) T-cell leukemia/lymphoma 1B, 1 [Mus musculus]
LPVYLVSVRLGIYEDEHHRVWIVANVETshSSHGNRRRTHVTVHLWKLIPQQVIPFNplnydFL
>gi|27668591|ref|XP_234504.1|:(7-103) similar to Chain A, Crystal Structure Of Murine Tcl1
-PDRLWLWEKHVYLDEFRRSWLPIVIKSNGKFQVIMRQKDVILGDSMTPSQLVPYEL
>gi|27668589|ref|XP_234503.1|:(9-91) similar to T-cell leukemia/lymphoma 1B, 5;
-PHILTLRTHGIYEDEHHRLWVVLDLQAShlSFSNRLLIYLTVYLQqgvafplESTPPSPMNLNGL
>gi|7305559|ref|NP_038802.1|:(8-102) T-cell leukemia/lymphoma 1B, 4 [Mus musculus]
PPCFLVCTRDDIYEDEHGRQWVAAKVETSSHSPycskietcvtVHLWQMTTLFQEPSPDSLKTFNFL
>gi|7305555|ref|NP_038803.1|:(9-102) T-cell leukemia/lymphoma 1B, 2 [Mus musculus]
---------PGFYEDEHHRLWMVAKLETCSHSPycnkietcvtVHLWQMTRYPQEPAPYNPMNYNFL
"""
        f_name = create_tmp_f(content=msa)
        parser = A3mIO()
        with open(f_name, 'r') as f_in:
            sequence_file = parser.read(f_in, remove_insert=True)     # <------------
        for i, sequence_entry in enumerate(sequence_file):
            if i == 0:
                self.assertEqual('d1a1x__ b.63.1.1 (-) p13-MTCP1 {Human (Homo sapiens)}',
                                 sequence_entry.id)
                self.assertEqual('PPDHLWVHQEGIYRDEYQRTWVAVVEEETSFLRARVQQIQVPLGDAARPSHLLTSQL',
                                 sequence_entry.seq)
            elif i == 1:
                self.assertEqual('gi|6678257|ref|NP_033363.1|:(7-103) T-cell lymphoma breakpoint 1 [Mus musculus]',
                                 sequence_entry.id)
                self.assertEqual('HPNRLWIWEKHVYLDEFRRSWLPVVIKSNEKFQVILRQEDVTLGEAMSPSQLVPYEL',
                                 sequence_entry.seq)
            elif i == 2:
                self.assertEqual('gi|7305557|ref|NP_038800.1|:(8-103) T-cell leukemia/lymphoma 1B, 3 [Mus musculus]',
                                 sequence_entry.id)
                self.assertEqual('PPRFLVCTRDDIYEDENGRQWVVAKVETSRSITVHLQHMTTIPQEPTPQQPINNNSL',
                                 sequence_entry.seq)
            elif i == 3:
                self.assertEqual('gi|11415028|ref|NP_068801.1|:(2-106) T-cell lymphoma-1; T-cell lymphoma-1A [Homo sapiens]',
                                 sequence_entry.id)
                self.assertEqual('HPDRLWAWEKFVYLDEKQHAWLPLTIEDRLQLRVLLRREDVVLGRPMTPTQIGPSLL',
                                 sequence_entry.seq)
            elif i == 4:
                self.assertEqual('gi|7305561|ref|NP_038804.1|:(7-103) T-cell leukemia/lymphoma 1B, 5 [Mus musculus]',
                                 sequence_entry.id)
                self.assertEqual('----------GIYEDEHHRVWIAVNVETSHSSHVHLQHMTTLPQEPTPQQPINNNSL',
                                 sequence_entry.seq)
            elif i == 5:
                self.assertEqual('gi|7305553|ref|NP_038801.1|:(5-103) T-cell leukemia/lymphoma 1B, 1 [Mus musculus]',
                                 sequence_entry.id)
                self.assertEqual('LPVYLVSVRLGIYEDEHHRVWIVANVETSSHGNRRRTHVTVHLWKLIPQQVIPFNFL',
                                 sequence_entry.seq)
            elif i == 6:
                self.assertEqual('gi|27668591|ref|XP_234504.1|:(7-103) similar to Chain A, Crystal Structure Of Murine Tcl1',
                                 sequence_entry.id)
                self.assertEqual('-PDRLWLWEKHVYLDEFRRSWLPIVIKSNGKFQVIMRQKDVILGDSMTPSQLVPYEL',
                                 sequence_entry.seq)
            elif i == 7:
                self.assertEqual('gi|27668589|ref|XP_234503.1|:(9-91) similar to T-cell leukemia/lymphoma 1B, 5;',
                                 sequence_entry.id)
                self.assertEqual('-PHILTLRTHGIYEDEHHRLWVVLDLQASSFSNRLLIYLTVYLQESTPPSPMNLNGL',
                                 sequence_entry.seq)
            elif i == 8:
                self.assertEqual('gi|7305559|ref|NP_038802.1|:(8-102) T-cell leukemia/lymphoma 1B, 4 [Mus musculus]',
                                 sequence_entry.id)
                self.assertEqual('PPCFLVCTRDDIYEDEHGRQWVAAKVETSSHSPVHLWQMTTLFQEPSPDSLKTFNFL',
                                 sequence_entry.seq)
            elif i == 9:
                self.assertEqual('gi|7305555|ref|NP_038803.1|:(9-102) T-cell leukemia/lymphoma 1B, 2 [Mus musculus]',
                                 sequence_entry.id)
                self.assertEqual('---------PGFYEDEHHRLWMVAKLETCSHSPVHLWQMTRYPQEPAPYNPMNYNFL',
                                 sequence_entry.seq)
        os.unlink(f_name)

    def test_read_2(self):
        msa = """>d1a1x__ b.63.1.1 (-) p13-MTCP1 {Human (Homo sapiens)}
PPDHLWVHQEGIYRDEYQRTWVAVVEEETSFLRARVQQIQVPLGDAARPSHLLTSQL
>gi|6678257|ref|NP_033363.1|:(7-103) T-cell lymphoma breakpoint 1 [Mus musculus]
HPNRLWIWEKHVYLDEFRRSWLPVVIKSNEKFQVILRQEDVTLGEAMSPSQLVPYEL
>gi|7305557|ref|NP_038800.1|:(8-103) T-cell leukemia/lymphoma 1B, 3 [Mus musculus]
PPRFLVCTRDDIYEDENGRQWVVAKVETSRSpygsrietcITVHLQHMTTIPQEPTPQQPINNNSL
>gi|11415028|ref|NP_068801.1|:(2-106) T-cell lymphoma-1; T-cell lymphoma-1A [Homo sapiens]
HPDRLWAWEKFVYLDEKQHAWLPLTIEikDRLQLRVLLRREDVVLGRPMTPTQIGPSLL
>gi|7305561|ref|NP_038804.1|:(7-103) T-cell leukemia/lymphoma 1B, 5 [Mus musculus]
----------GIYEDEHHRVWIAVNVETSHSSHgnrietcvtVHLQHMTTLPQEPTPQQPINNNSL
>gi|7305553|ref|NP_038801.1|:(5-103) T-cell leukemia/lymphoma 1B, 1 [Mus musculus]
LPVYLVSVRLGIYEDEHHRVWIVANVETshSSHGNRRRTHVTVHLWKLIPQQVIPFNplnydFL
>gi|27668591|ref|XP_234504.1|:(7-103) similar to Chain A, Crystal Structure Of Murine Tcl1
-PDRLWLWEKHVYLDEFRRSWLPIVIKSNGKFQVIMRQKDVILGDSMTPSQLVPYEL
>gi|27668589|ref|XP_234503.1|:(9-91) similar to T-cell leukemia/lymphoma 1B, 5;
-PHILTLRTHGIYEDEHHRLWVVLDLQAShlSFSNRLLIYLTVYLQqgvafplESTPPSPMNLNGL
>gi|7305559|ref|NP_038802.1|:(8-102) T-cell leukemia/lymphoma 1B, 4 [Mus musculus]
PPCFLVCTRDDIYEDEHGRQWVAAKVETSSHSPycskietcvtVHLWQMTTLFQEPSPDSLKTFNFL
>gi|7305555|ref|NP_038803.1|:(9-102) T-cell leukemia/lymphoma 1B, 2 [Mus musculus]
---------PGFYEDEHHRLWMVAKLETCSHSPycnkietcvtVHLWQMTRYPQEPAPYNPMNYNFL
"""
        f_name = create_tmp_f(content=msa)
        parser = A3mIO()
        with open(f_name, 'r') as f_in:
            sequence_file = parser.read(f_in, remove_insert=False)     # <------------
        for i, sequence_entry in enumerate(sequence_file):
            if i == 0:
                self.assertEqual('d1a1x__ b.63.1.1 (-) p13-MTCP1 {Human (Homo sapiens)}',
                                 sequence_entry.id)
                self.assertEqual('PPDHLWVHQEGIYRDEYQRTWVAVVEE--E--T--SF---------LR----------ARVQQIQVPLG-------DAARPSHLLTS-----QL',
                                 sequence_entry.seq)
            elif i == 1:
                self.assertEqual('gi|6678257|ref|NP_033363.1|:(7-103) T-cell lymphoma breakpoint 1 [Mus musculus]',
                                 sequence_entry.id)
                self.assertEqual('HPNRLWIWEKHVYLDEFRRSWLPVVIK--S--N--EK---------FQ----------VILRQEDVTLG-------EAMSPSQLVPY-----EL',
                                 sequence_entry.seq)
            elif i == 2:
                self.assertEqual('gi|7305557|ref|NP_038800.1|:(8-103) T-cell leukemia/lymphoma 1B, 3 [Mus musculus]',
                                 sequence_entry.id)
                self.assertEqual('PPRFLVCTRDDIYEDENGRQWVVAKVE--T--S--RSpygsrietcIT----------VHLQHMTTIPQ-------EPTPQQPINNN-----SL',
                                 sequence_entry.seq)
            elif i == 3:
                self.assertEqual(
                    'gi|11415028|ref|NP_068801.1|:(2-106) T-cell lymphoma-1; T-cell lymphoma-1A [Homo sapiens]',
                    sequence_entry.id)
                self.assertEqual('HPDRLWAWEKFVYLDEKQHAWLPLTIEikD--R--LQ---------LR----------VLLRREDVVLG-------RPMTPTQIGPS-----LL',
                                 sequence_entry.seq)
            elif i == 4:
                self.assertEqual('gi|7305561|ref|NP_038804.1|:(7-103) T-cell leukemia/lymphoma 1B, 5 [Mus musculus]',
                                 sequence_entry.id)
                self.assertEqual('----------GIYEDEHHRVWIAVNVE--T--S--HS---------SHgnrietcvt-VHLQHMTTLPQ-------EPTPQQPINNN-----SL',
                                 sequence_entry.seq)
            elif i == 5:
                self.assertEqual('gi|7305553|ref|NP_038801.1|:(5-103) T-cell leukemia/lymphoma 1B, 1 [Mus musculus]',
                                 sequence_entry.id)
                self.assertEqual('LPVYLVSVRLGIYEDEHHRVWIVANVE--TshS--SH---------GN----------RRRTHVTVHLW-------KLIPQQVIPFNplnydFL',
                                 sequence_entry.seq)
            elif i == 6:
                self.assertEqual(
                    'gi|27668591|ref|XP_234504.1|:(7-103) similar to Chain A, Crystal Structure Of Murine Tcl1',
                    sequence_entry.id)
                self.assertEqual('-PDRLWLWEKHVYLDEFRRSWLPIVIK--S--N--GK---------FQ----------VIMRQKDVILG-------DSMTPSQLVPY-----EL',
                                 sequence_entry.seq)
            elif i == 7:
                self.assertEqual('gi|27668589|ref|XP_234503.1|:(9-91) similar to T-cell leukemia/lymphoma 1B, 5;',
                                 sequence_entry.id)
                self.assertEqual('-PHILTLRTHGIYEDEHHRLWVVLDLQ--A--ShlSF---------SN----------RLLIYLTVYLQqgvafplESTPPSPMNLN-----GL',
                                 sequence_entry.seq)
            elif i == 8:
                self.assertEqual('gi|7305559|ref|NP_038802.1|:(8-102) T-cell leukemia/lymphoma 1B, 4 [Mus musculus]',
                                 sequence_entry.id)
                self.assertEqual('PPCFLVCTRDDIYEDEHGRQWVAAKVE--T--S--SH---------SPycskietcvtVHLWQMTTLFQ-------EPSPDSLKTFN-----FL',
                                 sequence_entry.seq)
            elif i == 9:
                self.assertEqual('gi|7305555|ref|NP_038803.1|:(9-102) T-cell leukemia/lymphoma 1B, 2 [Mus musculus]',
                                 sequence_entry.id)
                self.assertEqual('---------PGFYEDEHHRLWMVAKLE--T--C--SH---------SPycnkietcvtVHLWQMTRYPQ-------EPAPYNPMNYN-----FL',
                                 sequence_entry.seq)
        os.unlink(f_name)

    def test_read_3(self):
        msa = """>d1a1x__ b.63.1.1 (-) p13-MTCP1 {Human (Homo sapiens)}
PPDHLWVHQEGIYRDEYQRTWVAVVEEETSFLRARVQQIQVPLGDAARPSHLLTSQL
>gi|6678257|ref|NP_033363.1|:(7-103) T-cell lymphoma breakpoint 1 [Mus musculus]
HPNRLWIWEKHVYLDEFRRSWLPVVIKSNEKFQVILRQEDVTLGEAMSPSQLVPYEL
>gi|6678257|ref|NP_033363.1|:(7-103) T-cell lymphoma breakpoint 1 [Mus musculus]
HPNRLWIWEKHVYLDEFRRSWLPVVIKSNEKFQVILRQEDVTLGEAMSPSQLVPYEL
>gi|6678257|ref|NP_033363.1|:(7-103) T-cell lymphoma breakpoint 1 [Mus musculus]
HPNRLWIWEKHVYLDEFRRSWLPVVIKSNEKFQVILRQEDVTLGEAMSPSQLVPYEL
"""
        f_name = create_tmp_f(content=msa)
        parser = A3mIO()
        with open(f_name, 'r') as f_in:
            sequence_file = parser.read(f_in, remove_insert=False)  # <------------
        for i, sequence_entry in enumerate(sequence_file):
            if i == 0:
                self.assertEqual('d1a1x__ b.63.1.1 (-) p13-MTCP1 {Human (Homo sapiens)}',
                                 sequence_entry.id)
                self.assertEqual('PPDHLWVHQEGIYRDEYQRTWVAVVEEETSFLRARVQQIQVPLGDAARPSHLLTSQL',
                                 sequence_entry.seq)
            elif i == 1:
                self.assertEqual('gi|6678257|ref|NP_033363.1|:(7-103) T-cell lymphoma breakpoint 1 [Mus musculus]',
                                 sequence_entry.id[:80])
                self.assertEqual(79, len(sequence_entry.id))
                self.assertEqual('HPNRLWIWEKHVYLDEFRRSWLPVVIKSNEKFQVILRQEDVTLGEAMSPSQLVPYEL',
                                 sequence_entry.seq)
            elif i == 1:
                self.assertEqual('gi|6678257|ref|NP_033363.1|:(7-103) T-cell lymphoma breakpoint 1 [Mus musculus]',
                                 sequence_entry.id[:80])
                self.assertGreater(79, len(sequence_entry.id))
                self.assertEqual('HPNRLWIWEKHVYLDEFRRSWLPVVIKSNEKFQVILRQEDVTLGEAMSPSQLVPYEL',
                                 sequence_entry.seq)
            elif i == 1:
                self.assertEqual('gi|6678257|ref|NP_033363.1|:(7-103) T-cell lymphoma breakpoint 1 [Mus musculus]',
                                 sequence_entry.id[:80])
                self.assertGreater(79, len(sequence_entry.id))
                self.assertEqual('HPNRLWIWEKHVYLDEFRRSWLPVVIKSNEKFQVILRQEDVTLGEAMSPSQLVPYEL',
                                 sequence_entry.seq)
        os.unlink(f_name)

    def test_write_1(self):
        msa = [
            ">d1a1x__ b.63.1.1 (-) p13-MTCP1 {Human (Homo sapiens)}",
            "PPDHLWVHQEGIYRDEYQRTWVAVVEEETSFLRARVQQIQVPLGDAARPSHLLTSQL",
            ">gi|6678257|ref|NP_033363.1|:(7-103) T-cell lymphoma breakpoint 1 [Mus musculus]",
            "HPNRLWIWEKHVYLDEFRRSWLPVVIKSNEKFQVILRQEDVTLGEAMSPSQLVPYEL",
            ">gi|7305557|ref|NP_038800.1|:(8-103) T-cell leukemia/lymphoma 1B, 3 [Mus musculus]",
            "PPRFLVCTRDDIYEDENGRQWVVAKVETSRSpygsrietcITVHLQHMTTIPQEPTPQQPINNNSL",
            ">gi|11415028|ref|NP_068801.1|:(2-106) T-cell lymphoma-1; T-cell lymphoma-1A [Homo sapiens]",
            "HPDRLWAWEKFVYLDEKQHAWLPLTIEikDRLQLRVLLRREDVVLGRPMTPTQIGPSLL",
            ">gi|7305561|ref|NP_038804.1|:(7-103) T-cell leukemia/lymphoma 1B, 5 [Mus musculus]",
            "----------GIYEDEHHRVWIAVNVETSHSSHgnrietcvtVHLQHMTTLPQEPTPQQPINNNSL",
            ">gi|7305553|ref|NP_038801.1|:(5-103) T-cell leukemia/lymphoma 1B, 1 [Mus musculus]",
            "LPVYLVSVRLGIYEDEHHRVWIVANVETshSSHGNRRRTHVTVHLWKLIPQQVIPFNplnydFL",
            ">gi|27668591|ref|XP_234504.1|:(7-103) similar to Chain A, Crystal Structure Of Murine Tcl1",
            "-PDRLWLWEKHVYLDEFRRSWLPIVIKSNGKFQVIMRQKDVILGDSMTPSQLVPYEL",
            ">gi|27668589|ref|XP_234503.1|:(9-91) similar to T-cell leukemia/lymphoma 1B, 5;",
            "-PHILTLRTHGIYEDEHHRLWVVLDLQAShlSFSNRLLIYLTVYLQqgvafplESTPPSPMNLNGL",
            ">gi|7305559|ref|NP_038802.1|:(8-102) T-cell leukemia/lymphoma 1B, 4 [Mus musculus]",
            "PPCFLVCTRDDIYEDEHGRQWVAAKVETSSHSPycskietcvtVHLWQMTTLFQEPSPDSLKTFNFL",
            ">gi|7305555|ref|NP_038803.1|:(9-102) T-cell leukemia/lymphoma 1B, 2 [Mus musculus]",
            "---------PGFYEDEHHRLWMVAKLETCSHSPycnkietcvtVHLWQMTRYPQEPAPYNPMNYNFL",
        ]

        f_name_in = create_tmp_f(content=os.linesep.join(msa))
        f_name_out = create_tmp_f()
        parser = A3mIO()
        with open(f_name_in, 'r') as f_in, open(f_name_out, 'w') as f_out:
            sequence_file = parser.read(f_in, remove_insert=True)
            parser.write(f_out, sequence_file)
        ref = [
            ">d1a1x__ b.63.1.1 (-) p13-MTCP1 {Human (Homo sapiens)}",
            "PPDHLWVHQEGIYRDEYQRTWVAVVEEETSFLRARVQQIQVPLGDAARPSHLLTSQL",
            ">gi|6678257|ref|NP_033363.1|:(7-103) T-cell lymphoma breakpoint 1 [Mus musculus]",
            "HPNRLWIWEKHVYLDEFRRSWLPVVIKSNEKFQVILRQEDVTLGEAMSPSQLVPYEL",
            ">gi|7305557|ref|NP_038800.1|:(8-103) T-cell leukemia/lymphoma 1B, 3 [Mus musculus]",
            "PPRFLVCTRDDIYEDENGRQWVVAKVETSRSITVHLQHMTTIPQEPTPQQPINNNSL",
            ">gi|11415028|ref|NP_068801.1|:(2-106) T-cell lymphoma-1; T-cell lymphoma-1A [Homo sapiens]",
            "HPDRLWAWEKFVYLDEKQHAWLPLTIEDRLQLRVLLRREDVVLGRPMTPTQIGPSLL",
            ">gi|7305561|ref|NP_038804.1|:(7-103) T-cell leukemia/lymphoma 1B, 5 [Mus musculus]",
            "----------GIYEDEHHRVWIAVNVETSHSSHVHLQHMTTLPQEPTPQQPINNNSL",
            ">gi|7305553|ref|NP_038801.1|:(5-103) T-cell leukemia/lymphoma 1B, 1 [Mus musculus]",
            "LPVYLVSVRLGIYEDEHHRVWIVANVETSSHGNRRRTHVTVHLWKLIPQQVIPFNFL",
            ">gi|27668591|ref|XP_234504.1|:(7-103) similar to Chain A, Crystal Structure Of Murine Tcl1",
            "-PDRLWLWEKHVYLDEFRRSWLPIVIKSNGKFQVIMRQKDVILGDSMTPSQLVPYEL",
            ">gi|27668589|ref|XP_234503.1|:(9-91) similar to T-cell leukemia/lymphoma 1B, 5;",
            "-PHILTLRTHGIYEDEHHRLWVVLDLQASSFSNRLLIYLTVYLQESTPPSPMNLNGL",
            ">gi|7305559|ref|NP_038802.1|:(8-102) T-cell leukemia/lymphoma 1B, 4 [Mus musculus]",
            "PPCFLVCTRDDIYEDEHGRQWVAAKVETSSHSPVHLWQMTTLFQEPSPDSLKTFNFL",
            ">gi|7305555|ref|NP_038803.1|:(9-102) T-cell leukemia/lymphoma 1B, 2 [Mus musculus]",
            "---------PGFYEDEHHRLWMVAKLETCSHSPVHLWQMTRYPQEPAPYNPMNYNFL",
            "",
        ]
        ref = os.linesep.join(ref)
        with open(f_name_out, 'r') as f_in:
            output = "".join(f_in.readlines())
        self.assertEqual(ref, output)
        os.unlink(f_name_in)
        os.unlink(f_name_out)

    def test_write_2(self):
        msa = [
            ">d1a1x__ b.63.1.1 (-) p13-MTCP1 {Human (Homo sapiens)}",
            "PPDHLWVHQEGIYRDEYQRTWVAVVEEETSFLRARVQQIQVPLGDAARPSHLLTSQL",
            ">gi|6678257|ref|NP_033363.1|:(7-103) T-cell lymphoma breakpoint 1 [Mus musculus]",
            "HPNRLWIWEKHVYLDEFRRSWLPVVIKSNEKFQVILRQEDVTLGEAMSPSQLVPYEL",
            ">gi|7305557|ref|NP_038800.1|:(8-103) T-cell leukemia/lymphoma 1B, 3 [Mus musculus]",
            "PPRFLVCTRDDIYEDENGRQWVVAKVETSRSpygsrietcITVHLQHMTTIPQEPTPQQPINNNSL",
            ">gi|11415028|ref|NP_068801.1|:(2-106) T-cell lymphoma-1; T-cell lymphoma-1A [Homo sapiens]",
            "HPDRLWAWEKFVYLDEKQHAWLPLTIEikDRLQLRVLLRREDVVLGRPMTPTQIGPSLL",
            ">gi|7305561|ref|NP_038804.1|:(7-103) T-cell leukemia/lymphoma 1B, 5 [Mus musculus]",
            "----------GIYEDEHHRVWIAVNVETSHSSHgnrietcvtVHLQHMTTLPQEPTPQQPINNNSL",
            ">gi|7305553|ref|NP_038801.1|:(5-103) T-cell leukemia/lymphoma 1B, 1 [Mus musculus]",
            "LPVYLVSVRLGIYEDEHHRVWIVANVETshSSHGNRRRTHVTVHLWKLIPQQVIPFNplnydFL",
            ">gi|27668591|ref|XP_234504.1|:(7-103) similar to Chain A, Crystal Structure Of Murine Tcl1",
            "-PDRLWLWEKHVYLDEFRRSWLPIVIKSNGKFQVIMRQKDVILGDSMTPSQLVPYEL",
            ">gi|27668589|ref|XP_234503.1|:(9-91) similar to T-cell leukemia/lymphoma 1B, 5;",
            "-PHILTLRTHGIYEDEHHRLWVVLDLQAShlSFSNRLLIYLTVYLQqgvafplESTPPSPMNLNGL",
            ">gi|7305559|ref|NP_038802.1|:(8-102) T-cell leukemia/lymphoma 1B, 4 [Mus musculus]",
            "PPCFLVCTRDDIYEDEHGRQWVAAKVETSSHSPycskietcvtVHLWQMTTLFQEPSPDSLKTFNFL",
            ">gi|7305555|ref|NP_038803.1|:(9-102) T-cell leukemia/lymphoma 1B, 2 [Mus musculus]",
            "---------PGFYEDEHHRLWMVAKLETCSHSPycnkietcvtVHLWQMTRYPQEPAPYNPMNYNFL",
        ]

        f_name_in = create_tmp_f(content=os.linesep.join(msa))
        f_name_out = create_tmp_f()
        parser = A3mIO()
        with open(f_name_in, 'r') as f_in, open(f_name_out, 'w') as f_out:
            sequence_file = parser.read(f_in, remove_insert=False)
            parser.write(f_out, sequence_file)
        ref = [
            ">d1a1x__ b.63.1.1 (-) p13-MTCP1 {Human (Homo sapiens)}",
            "PPDHLWVHQEGIYRDEYQRTWVAVVEE--E--T--SF---------LR----------ARVQQIQVPLG-------DAARPSHLLTS-----QL",
            ">gi|6678257|ref|NP_033363.1|:(7-103) T-cell lymphoma breakpoint 1 [Mus musculus]",
            "HPNRLWIWEKHVYLDEFRRSWLPVVIK--S--N--EK---------FQ----------VILRQEDVTLG-------EAMSPSQLVPY-----EL",
            ">gi|7305557|ref|NP_038800.1|:(8-103) T-cell leukemia/lymphoma 1B, 3 [Mus musculus]",
            "PPRFLVCTRDDIYEDENGRQWVVAKVE--T--S--RSpygsrietcIT----------VHLQHMTTIPQ-------EPTPQQPINNN-----SL",
            ">gi|11415028|ref|NP_068801.1|:(2-106) T-cell lymphoma-1; T-cell lymphoma-1A [Homo sapiens]",
            "HPDRLWAWEKFVYLDEKQHAWLPLTIEikD--R--LQ---------LR----------VLLRREDVVLG-------RPMTPTQIGPS-----LL",
            ">gi|7305561|ref|NP_038804.1|:(7-103) T-cell leukemia/lymphoma 1B, 5 [Mus musculus]",
            "----------GIYEDEHHRVWIAVNVE--T--S--HS---------SHgnrietcvt-VHLQHMTTLPQ-------EPTPQQPINNN-----SL",
            ">gi|7305553|ref|NP_038801.1|:(5-103) T-cell leukemia/lymphoma 1B, 1 [Mus musculus]",
            "LPVYLVSVRLGIYEDEHHRVWIVANVE--TshS--SH---------GN----------RRRTHVTVHLW-------KLIPQQVIPFNplnydFL",
            ">gi|27668591|ref|XP_234504.1|:(7-103) similar to Chain A, Crystal Structure Of Murine Tcl1",
            "-PDRLWLWEKHVYLDEFRRSWLPIVIK--S--N--GK---------FQ----------VIMRQKDVILG-------DSMTPSQLVPY-----EL",
            ">gi|27668589|ref|XP_234503.1|:(9-91) similar to T-cell leukemia/lymphoma 1B, 5;",
            "-PHILTLRTHGIYEDEHHRLWVVLDLQ--A--ShlSF---------SN----------RLLIYLTVYLQqgvafplESTPPSPMNLN-----GL",
            ">gi|7305559|ref|NP_038802.1|:(8-102) T-cell leukemia/lymphoma 1B, 4 [Mus musculus]",
            "PPCFLVCTRDDIYEDEHGRQWVAAKVE--T--S--SH---------SPycskietcvtVHLWQMTTLFQ-------EPSPDSLKTFN-----FL",
            ">gi|7305555|ref|NP_038803.1|:(9-102) T-cell leukemia/lymphoma 1B, 2 [Mus musculus]",
            "---------PGFYEDEHHRLWMVAKLE--T--C--SH---------SPycnkietcvtVHLWQMTRYPQ-------EPAPYNPMNYN-----FL",
            "",
        ]
        ref = os.linesep.join(ref)
        with open(f_name_out, 'r') as f_in:
            output = "".join(f_in.readlines())
        self.assertEqual(ref, output)
        os.unlink(f_name_in)
        os.unlink(f_name_out)


if __name__ == "__main__":
    unittest.main(verbosity=2)
