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
    A = 'ALA'
    C = 'CYS'
    B = 'ASX'
    E = 'GLU'
    D = 'ASP'
    G = 'GLY'
    F = 'PHE'
    I = 'ILE'
    H = 'HIS'
    K = 'LYS'
    J = 'XLE'
    M = 'MET'
    L = 'LEU'
    O = 'PYL'
    N = 'ASN'
    Q = 'GLN'
    P = 'PRO'
    S = 'SER'
    R = 'ARG'
    U = 'SEC'
    T = 'THR'
    W = 'TRP'
    V = 'VAL'
    Y = 'TYR'
    X = 'XAA'
    Z = 'GLX'


class AminoAcidThreeToOne(Enum):
    """Amino acid mapping to convert three-letter codes to one-letter codes"""
    ALA = 'A'
    ARG = 'R'
    ASN = 'N'
    ASP = 'D'
    CME = 'C'
    CYS = 'C'
    GLN = 'Q'
    GLU = 'E'
    GLY = 'G'
    HIS = 'H'
    ILE = 'I'
    LEU = 'L'
    LYS = 'K'
    MET = 'M'
    MSE = 'M'
    PHE = 'F'
    PRO = 'P'
    PYL = 'O'
    SER = 'S'
    SEC = 'U'
    THR = 'T'
    TRP = 'W'
    TYR = 'Y'
    VAL = 'V'
    ASX = 'B'
    GLX = 'Z'
    XAA = 'X'
    UNK = 'X'
    XLE = 'J'


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
