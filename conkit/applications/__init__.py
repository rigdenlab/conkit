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
"""Python wrappers for common command line applications"""

__author__ = "Felix Simkovic"
__date__ = "20 Oct 2016"
__version__ = "0.1"


def BbcontactsCommandline(*args, **kwargs):
    from conkit.applications.bbcontacts import BbcontactsCommandline
    return BbcontactsCommandline(*args, **kwargs)


def CCMpredCommandline(*args, **kwargs):
    from conkit.applications.ccmpred import CCMpredCommandline
    return CCMpredCommandline(*args, **kwargs)


def CdhitCommandline(*args, **kwargs):
    from conkit.applications.cdhit import CdhitCommandline
    return CdhitCommandline(*args, **kwargs)


def HHblitsCommandline(*args, **kwargs):
    from conkit.applications.hhblits import HHblitsCommandline
    return HHblitsCommandline(*args, **kwargs)


def HHfilterCommandline(*args, **kwargs):
    from conkit.applications.hhfilter import HHfilterCommandline
    return HHfilterCommandline(*args, **kwargs)


def JackhmmerCommandline(*args, **kwargs):
    from conkit.applications.jackhmmer import JackhmmerCommandline
    return JackhmmerCommandline(*args, **kwargs)


def PsicovCommandline(*args, **kwargs):
    from conkit.applications.psicov import PsicovCommandline
    return PsicovCommandline(*args, **kwargs)
