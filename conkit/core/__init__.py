"""Core modules for hierarchy construction"""

__author__ = "Felix Simkovic"
__date__ = "03 Aug 2016"
__version__ = 0.1

from conkit.core.ContactCore import Contact
from conkit.core.ContactMapCore import ContactMap
from conkit.core.ContactFileCore import ContactFile
from conkit.core.SequenceCore import Sequence
from conkit.core.SequenceFileCore import SequenceFile

__all__ = [
    'Contact',
    'ContactMap',
    'ContactFile',
    'Sequence',
    'SequenceFile',
]
