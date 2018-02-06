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
Parent classes for all parser classes
"""

__author__ = "Felix Simkovic"
__date__ = "04 Oct 2016"
__version__ = "0.1"

import abc

ABC = abc.ABCMeta('ABC', (object,), {})

from conkit.core.contact import Contact
from conkit.core.contactmap import ContactMap
from conkit.core.contactfile import ContactFile
from conkit.core.sequence import Sequence
from conkit.core.sequencefile import SequenceFile


class Parser(ABC):
    """Abstract class for all parsers

    """
    @abc.abstractmethod
    def read(self):
        pass

    @abc.abstractmethod
    def write(self):
        pass

    @classmethod
    def _reconstruct(cls, hierarchy):
        """Wrapper to re-construct full hierarchy when parts are provided"""
        if isinstance(hierarchy, ContactFile):
            h = hierarchy
        elif isinstance(hierarchy, ContactMap):
            h = ContactFile('conkit')
            h.add(hierarchy)
        elif isinstance(hierarchy, Contact):
            h = ContactFile('conkit')
            m = ContactMap('1')
            m.add(hierarchy)
            h.add(m)
        elif isinstance(hierarchy, SequenceFile):
            h = hierarchy
        elif isinstance(hierarchy, Sequence):
            h = SequenceFile('conkit')
            h.add(hierarchy)
        return h 
    

class ContactFileParser(Parser):
    """General purpose class for all contact file parsers"""
    pass
    

class SequenceFileParser(Parser):
    """General purpose class for all sequence file parsers"""
    pass
