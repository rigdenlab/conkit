"""
Parser module specific to Pdb files
"""

from __future__ import division

__author__ = "Felix Simkovic"
__date__ = "27 Sep 2016"
__version__ = 0.1

import collections
import itertools
import warnings

from Bio.PDB import PDBParser

from conkit import constants
from conkit.core import Contact
from conkit.core import ContactFile
from conkit.core import ContactMap
from conkit.core import Sequence
from conkit.io._ParserIO import _ContactFileParser


class PdbParser(_ContactFileParser):
    """
    Class to parse a PDB file and extract distance restraints
    as residue-residue contacts
    """
    def __init__(self):
        super(PdbParser, self).__init__()

    def _build_sequence(self, chain):
        """Build a peptide using Biopython to extract the sequence"""
        return Sequence(chain.id + '_seq', ''.join(constants.THREE_TO_ONE[residue.resname] for residue in chain))

    def _chain_contacts(self, chain1, chain2):
        """Determine the contact pairs intra- or inter-molecular

        Parameters
        ----------
        chain1 : :obj:`Bio.PDB.Chain`
           A first chain object
        chain2 : :obj:`Bio.PDB.Chain`
           A second chain object

        Yields
        ------
        atom_comb : list
           A list of tuples containing the contact information

        """
        Atom = collections.namedtuple('Atom', 'resname resseq resseq_alt reschain')

        if chain1 == chain2:
            range1 = range2 = list(range(1, len(chain1) + 1))
        else:
            range1 = list(range(1, len(chain1) + 1))
            range2 = list(range(len(chain1) + 1, len(chain1)+len(chain2) + 1))

        assert len(range1) == len(chain1)
        assert len(range2) == len(chain2)

        for (resseq1_alt, residue1), (resseq2_alt, residue2) in itertools.product(list(zip(range1, chain1)),
                                                                                  list(zip(range2, chain2))):
            for atom1, atom2 in itertools.product(residue1, residue2):

                # Ignore duplicates
                if int(residue1.id[1]) > int(residue2.id[1]):
                    continue
                elif int(residue1.id[1]) == int(residue2.id[1]):
                    continue

                # Biopython implementation to calculate distance between atoms
                distance = atom1 - atom2

                construct1 = Atom(
                    resname=residue1.resname,
                    resseq=int(residue1.id[1]),
                    resseq_alt=resseq1_alt,
                    reschain=chain1.id,
                )
                construct2 = Atom(
                    resname=residue2.resname,
                    resseq=int(residue2.id[1]),
                    resseq_alt=resseq2_alt,
                    reschain=chain2.id,
                )
                yield (construct1, construct2, distance)

    def _remove_atom(self, chain, type):
        """Tidy up a chain removing all HETATM entries"""
        for residue in chain.copy():
            for atom in residue.copy():
                if atom.is_disordered():
                    chain[residue.id].detach_child(atom.id)
                elif residue.resname == 'GLY' and type == 'CB' and atom.id == 'CA':
                    continue
                elif atom.id != type:
                    chain[residue.id].detach_child(atom.id)

    def _remove_hetatm(self, chain):
        """Tidy up a chain removing all HETATM entries"""
        for residue in chain.copy():
            if residue.id[0].strip() and residue.resname not in constants.THREE_TO_ONE:
                chain.detach_child(residue.id)

    def read(self, f_handle, f_id="pdb", distance_cutoff=8, atom_type='CB'):
        """Read a contact file

        Parameters
        ----------
        f_handle
           Open file handle [read permissions]
        f_id : str, optional
           Unique contact file identifier
        distance_cutoff : int, optional
           Distance cutoff for which to determine contacts [default: 8]
        atom_type : str, optional
           Atom type between which distances are calculated [default: CB]

        Returns
        -------
        :obj:`ContactFile <conkit.core.ContactFile>`

        """
        structure = PDBParser(QUIET=True).get_structure("pdb", f_handle)

        hierarchies = []
        for model in structure:

            # Create a per-model hierarchy
            hierarchy = ContactFile(f_id + '_' + str(model.id))
            chains = list(chain for chain in model)

            # Tidy the chains of this model
            for chain in chains:
                self._remove_hetatm(chain)
                self._remove_atom(chain, atom_type)

            for chain1, chain2 in itertools.product(chains, chains):

                # Define the ContactMap Id
                if chain1.id == chain2.id:                          # intra
                    contact_map = ContactMap(chain1.id)
                else:                                               # inter
                    contact_map = ContactMap(chain1.id + chain2.id)

                # Find all contacts between the two chains
                for (atom1, atom2, distance) in self._chain_contacts(chain1, chain2):
                    contact = Contact(
                        atom1.resseq,
                        atom2.resseq,
                        round(1.0-(distance/100), 6),
                        distance_bound=(0, distance_cutoff)
                    )

                    contact.res1_altseq = atom1.resseq_alt
                    contact.res2_altseq = atom2.resseq_alt
                    contact.res1 = atom1.resname
                    contact.res2 = atom2.resname
                    contact.res1_chain = atom1.reschain
                    contact.res2_chain = atom2.reschain

                    if distance < distance_cutoff:
                        contact.define_true_positive()
                    else:
                        contact.define_false_positive()

                    contact_map.add(contact)

                # Tidy up empty maps
                if len(contact_map) > 0:
                    # Get the full length of the peptide sequence and store it
                    if len(contact_map.id) == 1:                                                # INTRA !!!
                        contact_map.sequence = self._build_sequence(chain1)
                        assert len(contact_map.sequence.seq) == len(chain1)                     # Check that ncon analyzed == len_chains
                    else:                                                                       # INTER !!!
                        contact_map.sequence = self._build_sequence(chain1) + self._build_sequence(chain2)
                        assert len(contact_map.sequence.seq) == len(chain1) + len(chain2)       # Check that ncon analyzed == len_chains
                    hierarchy.add(contact_map)                                                  # Save the map into the hierarchy
                else:
                    del contact_map                                                             # Delete the empty contact map

            hierarchy.method = 'Contact map extracted from PDB {0}'.format(model.id)
            hierarchy.remark = ['The model id is the chain identifier, i.e XY equates to chain X and chain Y.',
                                'Residue numbers in column 1 are chain X, and numbers in column 2 are chain Y.']
            hierarchies.append(hierarchy)

        if len(hierarchies) > 1:
            msg = "Super-level to contact file not yet implemented. " \
                  "Parser returns hierarchy for top model only!"
            warnings.warn(msg, FutureWarning)
        return hierarchies[0]

    def write(self, f_handle, hierarchy):
        """Write a contact file instance to to file

        Parameters
        ----------
        f_handle
           Open file handle [write permissions]
        hierarchy : :obj:`ContactFile <conkit.core.ContactFile>`, :obj:`ContactMap <conkit.core.ContactMap>`
                    or :obj:`ContactMap <conkit.core.Contact>`

        Raises
        ------
        RuntimeError
           Not available

        """
        raise RuntimeError("Not available")
