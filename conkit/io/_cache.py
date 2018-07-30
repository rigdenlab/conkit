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
"""Module containing code for caching all parsers

Description
-----------
This module is an import caching facility used for the I/O package in ConKit.

To allow fast access to individual modules required for :func:`read <conkit.io.read>`, :func:`write <conkit.io.write>`
and :func:`convert <conkit.io.convert>` functions, we don't want to import everything every time. 
Thus, we only cache the results to ultimately import the bits we really require.

"""

__author__ = 'Felix Simkovic'
__date__ = '19 Jun 2017'
__version__ = "1.0"

import copy
import glob
import importlib
import os
import re

RE_CLASS_DECLARATION = re.compile("^class\s+([A-Za-z0-9]+)\s*(.*):$")


class _CacheObj(object):
    """Container for individual module"""

    __slots__ = ["id", "module", "object", "group"]

    def __init__(self, id, module, object, group):
        self.id = id
        self.module = module
        self.object = object
        self.group = group

    def __repr__(self):
        return "{name}(id={id} module={module} object={object} group={group}".format(
            name=self.__class__.__name__, **{k: getattr(self, k)
                                             for k in self.__class__.__slots__})


class _ParserCache(object):
    """Cache to hold handlers to each file parser"""

    # This mask is for backward-compatibility and extensions to avoid re-writing the same algorithms
    MASKS = {
        'a2m': ['a2m', 'jones'],
        'a3m': ['a3m', 'a3m-inserts'],
        'casp': ['casp', 'casprr'],
        'pcons': ['flib', 'pconsc', 'pconsc2', 'pconsc3'],
        'psicov': ['psicov', 'metapsicov', 'nebcon'],
    }

    BLINDFOLD = set(['ContactFileParser', 'GenericStructureParser', 'SequenceFileParser'])

    def __init__(self):
        self._parsers = []

        self.__construct()

    def __contains__(self, item):
        return item in self.file_parsers

    def __getitem__(self, item):
        return self.file_parsers.get(item, None)

    def __repr__(self):
        return "{0}(nparsers={1})".format(self.__class__.__name__, len(self._cfile_parsers) + len(self._sfile_parsers))

    @property
    def contact_file_parsers(self):
        return {c.id: c for c in self._parsers if c.group in ["ContactFileParser", "GenericStructureParser"]}

    @property
    def sequence_file_parsers(self):
        return {c.id: c for c in self._parsers if c.group in ["SequenceFileParser"]}

    @property
    def file_parsers(self):
        return {c.id: c for c in self._parsers}

    def __construct(self):
        path = os.path.abspath(os.path.dirname(__file__))
        for m in glob.glob(os.path.join(path, "[!_]*.py")):
            with open(m, "r") as f_in:
                lines = [RE_CLASS_DECLARATION.match(l.strip()) for l in f_in if RE_CLASS_DECLARATION.match(l.strip())]
            for match in lines:
                decl = _CacheObj(None, None, match.group(1), None)
                ggroup = match.group(2)
                if ggroup and ggroup.startswith("(") and ggroup.endswith(")"):
                    decl.group = ggroup.replace("(", "").replace(")", "")
                else:
                    decl.group = "Ungrouped"
                name = os.path.basename(m).replace(".py", "")
                decl.module = "conkit.io." + name
                objname = decl.object.lower().replace("parser", "")
                if decl.object in _ParserCache.BLINDFOLD:
                    continue
                elif objname in _ParserCache.MASKS:
                    for extra in _ParserCache.MASKS[objname]:
                        decl_ = copy.copy(decl)
                        decl_.id = extra
                        self._parsers += [decl_]
                else:
                    decl.id = objname
                    self._parsers += [decl]

    def import_module(self, format):
        return importlib.import_module(self[format].module)

    def import_class(self, format):
        return getattr(self.import_module(format), PARSER_CACHE[format].object)


# Only allow this to be seen from outside
PARSER_CACHE = _ParserCache()
