# MIT License
#
# Copyright (c) 2017-18 Felix Simkovic
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

__author__ = "Felix Simkovic"
__date__ = "13 Aug 2018"
__version__ = "1.0"

from conkit.io._parser import ContactFileParser
from conkit.misc.distances import DynamicDistances
from conkit.misc.energyfunction import RosettaFunctionConstructs


class RosettaParser(ContactFileParser):
    """Implementation of a ROSETTA restraint file parser"""

    def read(self, f_handle, f_id="rosetta"):
        """Read a contact file into a :obj:`~conkit.core.contactfile.ContactFile` instance

        Parameters
        ----------
        f_handle
           Open file handle [read permissions]
        f_id : str, optional
           Unique contact file identifier

        Returns
        -------
        :obj:`~conkit.core.contactfile.ContactFile`

        Raises
        ------
        :exc:`NotImplementedError`

        """
        raise NotImplementedError

    def write(self, f_handle, hierarchy, efunc="FADE"):
        """Write a contact file instance to to file

        Parameters
        ----------
        f_handle
           Open file handle [write permissions]
        hierarchy : :obj:`~conkit.core.contactfile.ContactFile`, :obj:`~conkit.core.contactmap.ContactMap`
                    or :obj:`~conkit.core.contact.Contact`
        efunc : str, optional
           The output format

        """
        if not hasattr(RosettaFunctionConstructs, efunc):
            raise ValueError("Unknown Rosetta energy function: {}".format(efunc))
        contact_file = self._reconstruct(hierarchy)
        construct = getattr(RosettaFunctionConstructs, efunc).fget(RosettaFunctionConstructs)
        for contact in contact_file.top:
            contact_dict = contact._to_dict()
            contact_dict["atom1"] = "CA" if contact.res1 == "G" else "CB"
            contact_dict["atom2"] = "CA" if contact.res2 == "G" else "CB"
            contact_dict["energy_bonus"] = contact.weight * 15.00
            contact_dict["scalar_score"] = contact.scalar_score * contact.weight
            contact_dict["sigmoid_cutoff"] = DynamicDistances.cutoff(contact.res1, contact.res2)
            contact_dict["sigmoid_slope"] = 1.0 / DynamicDistances.percentile(contact.res1, contact.res2)
            f_handle.write(construct.format(**contact_dict) + "\n")
