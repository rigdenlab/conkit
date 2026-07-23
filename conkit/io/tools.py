# BSD 3-Clause License
#
# Copyright (c) 2016-21, University of Liverpool
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
"""Internal utility functions"""

__author__ = "Felix Simkovic"
__date__ = "16 Feb 2017"
__version__ = "0.14.1"

import numpy as np
import os
import subprocess

from conkit.core.contact import Contact
from conkit.core.contactmap import ContactMap
from conkit.core.contactfile import ContactFile
from conkit.core.sequence import Sequence
from conkit.core.sequencefile import SequenceFile
from conkit.core.distogram import Distogram
from conkit.core.distancefile import DistanceFile

HierarchyIndex = {
    "Contact": Contact,
    "ContactMap": ContactMap,
    "ContactFile": ContactFile,
    "Sequence": Sequence,
    "SequenceFile": SequenceFile,
    "Distogram": Distogram,
    "DistanceFile": DistanceFile
}

def set_contact_definition(moltype,rep_atom=None,cutoff=None):

    if rep_atom!=None:
        if cutoff!=None:
            return rep_atom, cutoff
        else: return rep_atom, 10

    elif moltype=='Protein':
        rep_atom = "CB"
        if cutoff==None: cutoff = 8
        return rep_atom, cutoff
    
    elif moltype=='RNA' or moltype=='DNA':
        rep_atom = "C1'"
        if cutoff==None: cutoff = 12.5
        return rep_atom, cutoff

    else:
        raise ValueError('Molecule type not supported without explicit contact definition, set atleast the --rep_atom flag and considder setting --contact_dist')