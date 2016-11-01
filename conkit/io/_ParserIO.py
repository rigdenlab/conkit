"""
Parent classes for all parser classes
"""

__author__ = "Felix Simkovic"
__date__ = "04 Oct 2016"
__version__ = 0.1

from conkit.core import Contact
from conkit.core import ContactFile
from conkit.core import ContactMap
from conkit.core import Sequence
from conkit.core import SequenceFile


class _Parser(object):
    """General purpose class for all parsers

    Notes
    -----
    Do no instantiate this class directly

    """
    def __init__(self):
        pass


class _ContactFileParser(_Parser):
    """General purpose class for all contact file parsers

    """
    def __init__(self):
        super(_ContactFileParser, self).__init__()

    def _reconstruct(self, hierarchy):
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
        return h 
    

class _SequenceFileParser(_Parser):
    """General purpose class for all sequence file parsers

    """
    def __init__(self):
        super(_SequenceFileParser, self).__init__()

    def _reconstruct(self, hierarchy):
        """Wrapper to re-construct full hierarchy when parts are provided"""
        if isinstance(hierarchy, SequenceFile):
            h = hierarchy
        if isinstance(hierarchy, Sequence):
            h = SequenceFile('conkit')
            h.add(hierarchy)
        return h
