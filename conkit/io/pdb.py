# BSD 3-Clause License
#
# Copyright (c) 2016-18, University of Liverpool
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# * Redistributions of source code must retain the above copyright notice, this
#   list of conditions and the following disclaimer.
#
# * Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.
#
# * Neither the name of the copyright holder nor the names of its
#   contributors may be used to endorse or promote products derived from
#   this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""
Parser module specific to Pdb files
"""

from __future__ import division

__author__ = "Felix Simkovic"
__date__ = "27 Sep 2016"
__version__ = "0.1"

import collections
import itertools
import sys
import warnings

if sys.version_info.major < 3:
    from itertools import izip as zip

from Bio.PDB import MMCIFParser
from Bio.PDB import PDBParser

from conkit.io._parser import ContactFileParser
from conkit.core.contact import Contact
from conkit.core.contactmap import ContactMap
from conkit.core.contactfile import ContactFile
from conkit.core.sequence import Sequence
from conkit.core.mappings import AminoAcidThreeToOne

ATOM = collections.namedtuple('Atom', 'resname resseq resseq_alt reschain')


class GenericStructureParser(ContactFileParser):
    """
    Parent class to parse a PDB file and extract distance restraints
    as residue-residue contacts
    """

    def _build_sequence(self, chain):
        """Build a peptide using :mod:`biopython` to extract the sequence"""
        return Sequence(chain.id + '_seq', ''.join(AminoAcidThreeToOne[residue.resname].value for residue in chain))

    def _chain_contacts(self, chain1, chain2):
        """Determine the contact pairs intra- or inter-molecular

        Parameters
        ----------
        chain1 : :obj:`~Bio.PDB.Chain`
           A first chain object
        chain2 : :obj:`~Bio.PDB.Chain`
           A second chain object

        Yields
        ------
        list
           A list of tuples containing the contact information

        """
        if chain1 == chain2:
            range1 = range2 = list(range(1, len(chain1) + 1))
        else:
            range1 = list(range(1, len(chain1) + 1))
            range2 = list(range(len(chain1) + 1, len(chain1) + len(chain2) + 1))
        assert len(range1) == len(chain1) and len(range2) == len(chain2)

        iterator = itertools.product(list(zip(range1, chain1)), list(zip(range2, chain2)))
        for (resseq1_alt, residue1), (resseq2_alt, residue2) in iterator:
            for atom1, atom2 in itertools.product(residue1, residue2):
                if chain1 == chain2 and int(residue1.id[1]) >= int(residue2.id[1]):
                    continue
                construct1 = ATOM(
                    resname=residue1.resname,
                    resseq=int(residue1.id[1]),
                    resseq_alt=resseq1_alt,
                    reschain=chain1.id,
                )
                construct2 = ATOM(
                    resname=residue2.resname,
                    resseq=int(residue2.id[1]),
                    resseq_alt=resseq2_alt,
                    reschain=chain2.id,
                )
                yield (construct1, construct2, atom1 - atom2)

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
            if residue.id[0].strip() and residue.resname not in AminoAcidThreeToOne.__members__:
                chain.detach_child(residue.id)

    def _read(self, structure, f_id, distance_cutoff, atom_type):
        """Read a contact file

        Parameters
        ----------
        structure
           A :obj:`~Bio.PDB.Structure.Structure>` instance
        f_id : str
           Unique contact file identifier
        distance_cutoff : int
           Distance cutoff for which to determine contacts
        atom_type : str
           Atom type between which distances are calculated

        Returns
        -------
        :obj:`~conkit.core.contactfile.ContactFile~`

        """
        hierarchies = []
        for model in structure:
            hierarchy = ContactFile(f_id + '_' + str(model.id))
            chains = list(chain for chain in model)

            for chain in chains:
                self._remove_hetatm(chain)
                self._remove_atom(chain, atom_type)

            for chain1, chain2 in itertools.product(chains, chains):
                if chain1.id == chain2.id:  # intra
                    contact_map = ContactMap(chain1.id)
                else:  # inter
                    contact_map = ContactMap(chain1.id + chain2.id)

                for (atom1, atom2, distance) in self._chain_contacts(chain1, chain2):
                    contact = Contact(
                        atom1.resseq,
                        atom2.resseq,
                        round(1.0 - (distance / 100), 6),
                        distance_bound=(0., float(distance_cutoff)))

                    contact.res1_altseq = atom1.resseq_alt
                    contact.res2_altseq = atom2.resseq_alt
                    contact.res1 = atom1.resname
                    contact.res2 = atom2.resname
                    contact.res1_chain = atom1.reschain
                    contact.res2_chain = atom2.reschain

                    if distance_cutoff == 0 or distance < distance_cutoff:
                        contact.true_positive = True
                        contact_map.add(contact)

                if contact_map.empty:
                    del contact_map
                else:
                    if len(contact_map.id) == 1:
                        contact_map.sequence = self._build_sequence(chain1)
                        assert len(contact_map.sequence.seq) == len(chain1)
                    else:
                        contact_map.sequence = self._build_sequence(chain1) \
                            + self._build_sequence(chain2)
                        assert len(contact_map.sequence.seq) \
                            == len(chain1) + len(chain2)
                    hierarchy.add(contact_map)

            hierarchy.method = 'Contact map extracted from PDB ' + str(model.id)
            hierarchy.remark = [
                'The model id is the chain identifier, i.e XY equates to chain X and chain Y.',
                'Residue numbers in column 1 are chain X, and numbers in column 2 are chain Y.'
            ]
            hierarchies.append(hierarchy)

        if len(hierarchies) > 1:
            msg = "Super-level to contact file not yet implemented. " \
                  "Parser returns hierarchy for top model only!"
            warnings.warn(msg, FutureWarning)
        return hierarchies[0]

    def _write(self, f_handle, hierarchy):
        """Write a contact file instance to to file

        Raises
        ------
        :exc:`NotImplementedError`
           Write function not available

        """
        raise NotImplementedError("Write function not available")


class MmCifParser(GenericStructureParser):
    """
    Class to parse a mmCIF file and extract distance restraints
    as residue-residue contacts
    """

    def __init__(self):
        super(MmCifParser, self).__init__()

    def read(self, f_handle, f_id="mmcif", distance_cutoff=8, atom_type='CB'):
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
        :obj:`~conkit.core.contactfile.ContactFile`

        """
        structure = MMCIFParser(QUIET=True).get_structure("mmcif", f_handle)
        return self._read(structure, f_id, distance_cutoff, atom_type)

    def write(self, f_handle, hierarchy):
        """Write a contact file instance to to file

        Parameters
        ----------
        f_handle
           Open file handle [write permissions]
        hierarchy : :obj:`~conkit.core.contactfile.ContactFile`, :obj:`~conkit.core.contactmap.ContactMap`
                    or :obj:`~conkit.core.contact.Contact`

        Raises
        ------
        :exc:`NotImplementedError`
           Write function not available

        """
        self._write(f_handle, hierarchy)


class PdbParser(GenericStructureParser):
    """
    Class to parse a PDB file and extract distance restraints
    as residue-residue contacts
    """

    def __init__(self):
        super(PdbParser, self).__init__()

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
        :obj:`~conkit.core.contactfile.ContactFile`

        """
        structure = PDBParser(QUIET=True).get_structure("pdb", f_handle)
        return self._read(structure, f_id, distance_cutoff, atom_type)

    def write(self, f_handle, hierarchy):
        """Write a contact file instance to to file

        Parameters
        ----------
        f_handle
           Open file handle [write permissions]
        hierarchy : :obj:`~conkit.core.contactfile.ContactFile`, :obj:`~conkit.core.contactmap.ContactMap`
                    or :obj:`~conkit.core.contact.Contact`

        Raises
        ------
        :exc:`NotImplementedError`
           Write function not available

        """
        self._write(f_handle, hierarchy)
