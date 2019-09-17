# BSD 3-Clause License
#
# Copyright (c) 2016-19, University of Liverpool
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

__author__ = "Felix Simkovic"
__date__ = "19 Jun 2017"
__version__ = "1.0"

import collections
import copy
import glob
import importlib
import os
import re

RE_CLASS_DECLARATION = re.compile(r"^class\s+([A-Za-z0-9]+)\s*(.*):$")

CacheObj = collections.namedtuple("CacheObj", ["id", "module", "object", "group"])


class ParserCache(object):
    """Cache to hold handlers to each file parser"""

    # This mask is for backward-compatibility and extensions to avoid re-writing the same algorithms
    MASKS = {
        "a2m": ["a2m", "jones"],
        "a3m": ["a3m", "a3m-inserts"],
        "aleigen": ["aleigen"],
        "casp": ["casp", "casprr"],
        "mapalign": ["mapalign"],
        "pcons": ["flib", "pconsc", "pconsc2", "pconsc3", "saint2"],
        "psicov": ["psicov", "metapsicov", "nebcon"],
    }

    BLINDFOLD = set(["ContactFileParser", "GenericStructureParser", "SequenceFileParser"])

    def __init__(self):
        self._parsers = []

        self.__construct()

    def __contains__(self, item):
        return item in self.file_parsers

    def __getitem__(self, item):
        return self.file_parsers.get(item, None)

    def __repr__(self):
        nparsers = len(self.contact_file_parsers) + len(self.sequence_file_parsers)
        return "{}(nparsers={})".format(self.__class__.__name__, nparsers)

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
                object_ = match.group(1)
                ggroup = match.group(2)
                if ggroup and ggroup.startswith("(") and ggroup.endswith(")"):
                    group = ggroup.replace("(", "").replace(")", "")
                else:
                    group = "Ungrouped"

                name = os.path.basename(m).replace(".py", "")
                module = "conkit.io." + name
                objname = object_.lower().replace("parser", "")

                if object_ in ParserCache.BLINDFOLD:
                    continue
                elif objname in ParserCache.MASKS:
                    self._parsers.extend(
                        [CacheObj(extra, module, object_, group) for extra in ParserCache.MASKS[objname]]
                    )
                else:
                    self._parsers.append(CacheObj(objname, module, object_, group))

    def import_module(self, format):
        return importlib.import_module(self[format].module)

    def import_class(self, format):
        return getattr(self.import_module(format), PARSER_CACHE[format].object)


# Only allow this to be seen from outside
PARSER_CACHE = ParserCache()
