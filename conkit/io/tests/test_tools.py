"""Testing facility for conkit.io.tools"""

import unittest

from conkit.io.tools import set_contact_definition


class TestSetContactDefinition(unittest.TestCase):

    def test_protein_default(self):
        atom, cutoff = set_contact_definition('Protein')
        self.assertEqual('CB', atom)
        self.assertEqual(8, cutoff)

    def test_protein_custom_cutoff(self):
        atom, cutoff = set_contact_definition('Protein', cutoff=10)
        self.assertEqual('CB', atom)
        self.assertEqual(10, cutoff)

    def test_rna_default(self):
        atom, cutoff = set_contact_definition('RNA')
        self.assertEqual("C1'", atom)
        self.assertEqual(12.5, cutoff)

    def test_rna_custom_cutoff(self):
        atom, cutoff = set_contact_definition('RNA', cutoff=10)
        self.assertEqual("C1'", atom)
        self.assertEqual(10, cutoff)

    def test_dna_default(self):
        atom, cutoff = set_contact_definition('DNA')
        self.assertEqual("C1'", atom)
        self.assertEqual(12.5, cutoff)

    def test_dna_custom_cutoff(self):
        atom, cutoff = set_contact_definition('DNA', cutoff=10)
        self.assertEqual("C1'", atom)
        self.assertEqual(10, cutoff)

    def test_custom_rep_atom_and_cutoff(self):
        atom, cutoff = set_contact_definition('Custom', rep_atom='CA', cutoff=10)
        self.assertEqual('CA', atom)
        self.assertEqual(10, cutoff)

    def test_custom_rep_atom_default_cutoff(self):
        # When rep_atom is given but no cutoff, the function falls back to 10
        atom, cutoff = set_contact_definition('Custom', rep_atom='CA')
        self.assertEqual('CA', atom)
        self.assertEqual(10, cutoff)

    def test_unsupported_moltype_raises(self):
        with self.assertRaises(ValueError):
            set_contact_definition('Ligand')


if __name__ == "__main__":
    unittest.main(verbosity=2)
