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
"""Module containing code for caching all parsers"""

__author__ = 'Felix Simkovic'
__date__ = '19 Jun 2017'
__version__ = "1.0"


class _CacheObj(object):
    """Container for individual module"""

    __slots__ = ["id", "module", "object", "group"]

    def __init__(self, id, module, object, group):
        self.id = id
        self.module = module
        self.object = object
        self.group = group


class _ParserCache(object):
    """Cache to hold handlers to each file parser"""

    def __init__(self):
        self._cfile_parsers = []    # Contact file parsers
        self._sfile_parsers = []    # Sequence file parsers

        self.__construct()

    def __contains__(self, item):
        return item in self.file_parsers

    def __getitem__(self, item):
        return self.file_parsers.get(item, None)

    def __repr__(self):
        return "{0}(nparsers={1})".format(
            self.__class__.__name__, len(self._cfile_parsers) + len(self._sfile_parsers)
        )

    @property
    def contact_file_parsers(self):
        """A dict of contact file parsers"""
        # Mask as dictionaries to not break existing code
        return {c.id: c for c in self._cfile_parsers}

    @property
    def sequence_file_parsers(self):
        """A dict of sequence file parsers"""
        # Mask as dictionaries to not break existing code
        return {c.id: c for c in self._sfile_parsers}

    @property
    def file_parsers(self):
        """A dict of all file parsers"""
        # Mask as dictionaries to not break existing code
        return {c.id: c for c in self._sfile_parsers + self._cfile_parsers}

    def __construct(self):
        """Create the handles"""
        base = "conkit.io."

        self._cfile_parsers += [_CacheObj("bclcontact", base + "BCLContactIO", "BCLContactParser", "contact")]
        self._cfile_parsers += [_CacheObj("bbcontacts", base + "BbcontactsIO", "BbcontactsParser", "contact")]
        self._cfile_parsers += [_CacheObj("ccmpred", base + "CCMpredIO", "CCMpredParser", "contact")]
        self._cfile_parsers += [_CacheObj("casprr", base + "CaspIO", "CaspParser", "contact")]
        self._cfile_parsers += [_CacheObj("comsat", base + "ComsatIO", "ComsatParser", "contact")]
        self._cfile_parsers += [_CacheObj("epcmap", base + "EPCMapIO", "EPCMapParser", "contact")]
        self._cfile_parsers += [_CacheObj("evfold", base + "EVfoldIO", "EVfoldParser", "contact")]
        self._cfile_parsers += [_CacheObj("freecontact", base + "FreeContactIO", "FreeContactParser", "contact")]
        self._cfile_parsers += [_CacheObj("gremlin", base + "GremlinIO", "GremlinParser", "contact")]
        self._cfile_parsers += [_CacheObj("membrain", base + "MemBrainIO", "MemBrainParser", "contact")]
        self._cfile_parsers += [_CacheObj("pconsc", base + "PconsIO", "PconsParser", "contact")]
        self._cfile_parsers += [_CacheObj("pconsc2", base + "PconsIO", "PconsParser", "contact")]
        self._cfile_parsers += [_CacheObj("pconsc3", base + "PconsIO", "PconsParser", "contact")]
        self._cfile_parsers += [_CacheObj("pdb", base + "PdbIO", "PdbParser", "contact")]
        self._cfile_parsers += [_CacheObj("mmcif", base + "PdbIO", "MmCifParser", "contact")]
        self._cfile_parsers += [_CacheObj("plmdca", base + "PlmDCAIO", "PlmDCAParser", "contact")]
        self._cfile_parsers += [_CacheObj("psicov", base + "PsicovIO", "PsicovParser", "contact")]
        self._cfile_parsers += [_CacheObj("metapsicov", base + "PsicovIO", "PsicovParser", "contact")]

        self._sfile_parsers += [_CacheObj("a3m", base + "A3mIO", "A3mParser", "sequence")]
        self._sfile_parsers += [_CacheObj("a3m-inserts", base + "A3mIO", "A3mParser", "sequence")]
        self._sfile_parsers += [_CacheObj("fasta", base + "FastaIO", "FastaParser", "sequence")]
        self._sfile_parsers += [_CacheObj("jones", base + "JonesIO", "JonesParser", "sequence")]
        self._sfile_parsers += [_CacheObj("stockholm", base + "StockholmIO", "StockholmParser", "sequence")]


PARSER_CACHE = _ParserCache()
CONTACT_FILE_PARSERS = PARSER_CACHE.contact_file_parsers
SEQUENCE_FILE_PARSERS = PARSER_CACHE.sequence_file_parsers
