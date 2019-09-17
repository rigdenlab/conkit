# coding=utf-8
#
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
"""Mappings required for the core functionality"""

from __future__ import division
from __future__ import print_function

__author__ = "Felix Simkovic"
__date__ = "03 Aug 2016"
__version__ = "1.0"

from enum import Enum, unique


class AminoAcidMapping(Enum):
    """Amino acid mapping to encode an alignment"""

    A = 1
    C = 2
    D = 3
    E = 4
    F = 5
    G = 6
    H = 7
    I = 8
    K = 9
    L = 10
    M = 11
    N = 12
    P = 13
    Q = 14
    R = 15
    S = 16
    T = 17
    V = 18
    W = 19
    X = 21
    Y = 20


class AminoAcidOneToThree(Enum):
    """Amino acid mapping to convert one-letter codes to three-letter codes"""

    A = "ALA"
    C = "CYS"
    B = "ASX"
    E = "GLU"
    D = "ASP"
    G = "GLY"
    F = "PHE"
    I = "ILE"
    H = "HIS"
    K = "LYS"
    J = "XLE"
    M = "MET"
    L = "LEU"
    O = "PYL"
    N = "ASN"
    Q = "GLN"
    P = "PRO"
    S = "SER"
    R = "ARG"
    U = "SEC"
    T = "THR"
    W = "TRP"
    V = "VAL"
    Y = "TYR"
    X = "XAA"
    Z = "GLX"


class AminoAcidThreeToOne(Enum):
    """Amino acid mapping to convert three-letter codes to one-letter codes"""

    ALA = "A"
    ARG = "R"
    ASN = "N"
    ASP = "D"
    CME = "C"
    CYS = "C"
    GLN = "Q"
    GLU = "E"
    GLY = "G"
    HIS = "H"
    ILE = "I"
    LEU = "L"
    LYS = "K"
    MET = "M"
    MSE = "M"
    PHE = "F"
    PRO = "P"
    PYL = "O"
    SER = "S"
    SEC = "U"
    THR = "T"
    TRP = "W"
    TYR = "Y"
    VAL = "V"
    ASX = "B"
    GLX = "Z"
    XAA = "X"
    UNK = "X"
    XLE = "J"


@unique
class ContactMatchState(Enum):
    """Enumerated class to store state constants for each contact"""

    unknown = 0
    true_positive = 1
    true_negative = 2
    false_positive = 3
    false_negative = 4


@unique
class SequenceAlignmentState(Enum):
    """Alignment states"""

    unknown = 0
    unaligned = 1
    aligned = 2
