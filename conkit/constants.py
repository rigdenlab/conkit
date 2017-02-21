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

