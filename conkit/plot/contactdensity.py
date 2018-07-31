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
"""A module to produce a domain boundary plot"""

from __future__ import division
from __future__ import print_function

__author__ = "Felix Simkovic"
__date__ = "23 Feb 2017"
__version__ = "0.1"

import matplotlib.pyplot as plt
import numpy as np

from conkit.misc import deprecate
from conkit.plot.figure import Figure
from conkit.plot.tools import ColorDefinitions
from conkit.plot.tools import find_minima
from conkit.plot.tools import _isinstance


class ContactDensityFigure(Figure):
    """A Figure object specifically for a contact density illustration.

    This figure is an adaptation of the algorithm published by Sadowski
    (2013) [#]_.

    .. [#] Sadowski M. (2013). Prediction of protein domain boundaries
       from inverse covariances. Proteins 81(2), 253-260.

    Attributes
    ----------
    hierarchy : :obj:`~conkit.core.contactmap.ContactMap`
       The default contact map hierarchy
    bw_method : str
       The method to estimate the bandwidth

    Examples
    --------
    >>> import conkit
    >>> cmap = conkit.io.read('toxd/toxd.mat', 'ccmpred').top_map
    >>> conkit.plot.ContactDensityFigure(cmap)

    """

    def __init__(self, hierarchy, bw_method='bowman', **kwargs):
        """A new contact density plot

        Parameters
        ----------
        hierarchy : :obj:`~conkit.core.contactmap.ContactMap`
           The default contact map hierarchy
        bw_method : str, optional
           The method to estimate the bandwidth [default: bowman]
        **kwargs
           General :obj:`~conkit.plot.figure.Figure` keyword arguments

        """
        super(ContactDensityFigure, self).__init__(**kwargs)
        self._bw_method = None
        self._hierarchy = None

        self.bw_method = bw_method
        self.hierarchy = hierarchy

        self.minima_ = None

        self.draw()

    def __repr__(self):
        return self.__class__.__name__

    @property
    def bw_method(self):
        """The method to estimate the bandwidth
        
        For a full list of options, please refer to 
        :meth:`~conkit.core.contactmap.ContactMap.get_contact_density`

        """
        return self._bw_method

    @bw_method.setter
    def bw_method(self, bw_method):
        """Define the method to estimate the bandwidth"""
        self._bw_method = bw_method

    @property
    def hierarchy(self):
        """A :obj:`~conkit.core.contactmap.ContactMap`"""
        return self._hierarchy

    @hierarchy.setter
    def hierarchy(self, hierarchy):
        """Define the ConKit :obj:`ContactMap <conkit.core.contactmap.ContactMap>`

        Raises
        ------
        :exc:`TypeError`
           The hierarchy is not a :obj:`~conkit.core.contactmap.ContactMap`

        """
        if hierarchy and _isinstance(hierarchy, "ContactMap"):
            self._hierarchy = hierarchy
        else:
            raise TypeError("The hierarchy is not an contact map")

    @deprecate('0.11', msg='Use draw instead')
    def redraw(self):
        self.draw()

    def draw(self):
        x, y = self.get_xy_data()
        self.ax.plot(x, y, linestyle="solid", color=ColorDefinitions.GENERAL, label="Contact Density", zorder=2)
        line_kwargs = dict(linestyle="--", linewidth=1.0, alpha=0.5, color=ColorDefinitions.MISMATCH, zorder=1)
        self.minima_ = []
        for minimum in find_minima(y, order=1):
            self.minima_.append(x[minimum])
            self.ax.axvline(x[minimum], **line_kwargs)
        self.ax.axvline(0, ymin=0, ymax=0, label="Domain Boundary", **line_kwargs)
        self.ax.set_xlim(x.min(), x.max())
        self.ax.set_ylim(0., y.max())
        self.ax.set_xlabel('Residue number')
        self.ax.set_ylabel('Density Estimate')
        if self.legend:
            self.ax.legend(
                bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=3, mode="expand", borderaxespad=0., scatterpoints=1)
        # TODO: deprecate this in 0.10
        if self._file_name:
            self.savefig(self._file_name, dpi=self._dpi)

    def get_xy_data(self):
        residues = np.asarray(self.hierarchy.as_list()).flatten()
        x = np.arange(residues.min(), residues.max() + 1)
        y = np.asarray(self.hierarchy.get_contact_density(self.bw_method))
        return x, y
