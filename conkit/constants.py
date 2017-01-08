"""
Collection of constant variables used throught ConKit
"""

__author__ = "Felix Simkovic"
__date__ = "09 Sep 2016"
__version__ = 0.1

# ================================================
# Amino acid conversions
# ================================================
ONE_TO_THREE = {'A': 'ALA', 'C': 'CYS', 'B': 'ASX', 'E': 'GLU', 'D': 'ASP', 'G': 'GLY', 'F': 'PHE', 'I': 'ILE',
                'H': 'HIS', 'K': 'LYS', 'J': 'XLE', 'M': 'MET', 'L': 'LEU', 'O': 'PYL', 'N': 'ASN', 'Q': 'GLN',
                'P': 'PRO', 'S': 'SER', 'R': 'ARG', 'U': 'SEC', 'T': 'THR', 'W': 'TRP', 'V': 'VAL', 'Y': 'TYR',
                'X': 'XAA', 'Z': 'GLX'}

THREE_TO_ONE = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CME': 'C', 'CYS': 'C', 'GLN': 'Q', 'GLU': 'E', 
                'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'MSE': 'M', 'PHE': 'F', 
                'PRO': 'P', 'PYL': 'O', 'SER': 'S', 'SEC': 'U', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V', 
                'ASX': 'B', 'GLX': 'Z', 'XAA': 'X', 'XLE': 'J'}

# ================================================
# Constants for easier instance detection
# ================================================
CONTACTFILE = 1
CONTACTMAP = 2
CONTACT = 3
SEQUENCEFILE = 4
SEQUENCE = 5
UNKNOWN = 999
MATCHED = 0
UNMATCHED = -1
UNREGISTERED = -2

# ================================================
# Constants defining Contact().status color coding
# ================================================
TPCOLOR = '#2D9D00'     # color true positive
FPCOLOR = '#AB0000'     # color false positive
NTCOLOR = '#0482ff'     # color undefined
RFCOLOR = '#B5B5B5'     # color structure
