# BSD 3-Clause License
#
# Copyright (c) 2016-17, University of Liverpool
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
"""A module containing generic Figure related objects"""

__author__ = "Felix Simkovic"
__date__ = "08 Jan 2018"
__version__ = "0.2"

import os


class Figure(object):
    """A Figure class to store common features"""

    def __init__(self, ax=None, legend=False):
        """Initialise a new :obj:`conkit.plot.Figure` object
            
        Parameters
        ----------
        ax : 
        legend : 
        
        """
        if ax is None:
            import matplotlib.pyplot as plt
            self._fig, self._ax = plt.subplots()
        else:
            self._fig = None
            self._ax = ax
        self.legend = legend

    def __repr__(self):
        return self.__class__.__name__

    def savefig(self, filename, dpi=300, overwrite=False):
        if os.path.isfile(filename) and not overwrite:
            raise RuntimeError("File exists: %s! Please rename or remove." % filename)
        else:
            self._fig.savefig(filename, dpi=dpi, bbox_inches="tight")

    @property
    def fig(self):
        return self._fig

    @property
    def ax(self):
        return self._ax
