# coding=utf-8
#
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
"""Core entities for hierarchy construction"""

__author__ = "Felix Simkovic"
__contributing_authors__ = "Jens Thomas"
__date__ = "03 Aug 2016"
__version__ = "1.0"


def Contact(*args, **kwargs):
    """:obj:`Contact <conkit.core.Contact.Contact>` instance"""
    from conkit.core.contact import Contact
    return Contact(*args, **kwargs)


def ContactMap(*args, **kwargs):
    """:obj:`ContactMap <conkit.core.ContactMap.ContactMap>` instance"""
    from conkit.core.contactmap import ContactMap
    return ContactMap(*args, **kwargs)


def ContactFile(*args, **kwargs):
    """:obj:`ContactFile <conkit.core.ContactFile.ContactFile>` instance"""
    from conkit.core.contactfile import ContactFile
    return ContactFile(*args, **kwargs)


def Sequence(*args, **kwargs):
    """:obj:`Sequence <conkit.core.Sequence.Sequence>` instance"""
    from conkit.core.sequence import Sequence
    return Sequence(*args, **kwargs)


def SequenceFile(*args, **kwargs):
    """:obj:`SequenceFile <conkit.core.SequenceFile.SequenceFile>` instance"""
    from conkit.core.sequencefile import SequenceFile
    return SequenceFile(*args, **kwargs)
