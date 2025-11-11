"""Three to one symbol residue translation."""

PROTEIN_standard_rescodes = {
  "ALA": "A",
  "ARG": "R",
  "ASN": "N",
  "ASP": "D",
  "CYS": "C",
  "GLN": "Q",
  "GLU": "E",
  "GLY": "G",
  "HIS": "H",
  "ILE": "I",
  "LEU": "L",
  "LYS": "K",
  "MET": "M",
  "PHE": "F",
  "PRO": "P",
  "SER": "S",
  "THR": "T",
  "TRP": "W",
  "TYR": "Y",
  "VAL": "V",
}

RNA_standard_rescodes = {
  "A": "A",  
  "C": "C",  
  "G": "G",  
  "U": "U"
}

DNA_standard_rescodes = {
  "A": "A",  "DA": "A",
  "T": "T",  "DT": "T",
  "C": "C",  "DC": "C",
  "G": "G",  "DG": "G"
}

reslist = {
  "ALA": "A",
  "ARG": "R",
  "ASN": "N",
  "ASP": "D",
  "CYS": "C",
  "GLN": "Q",
  "GLU": "E",
  "GLY": "G",
  "HIS": "H",
  "ILE": "I",
  "LEU": "L",
  "LYS": "K",
  "MET": "M",
  "PHE": "F",
  "PRO": "P",
  "SER": "S",
  "THR": "T",
  "TRP": "W",
  "TYR": "Y",
  "VAL": "V",
  "SEC": "U",
  "PYL": "O"
}

mod_reslist = {
  "ABA": "A",  "AIB": "A",  "ALC": "A",  "AYA": "A",  "DAB": "A",  "MAA": "A",
  "ORN": "A",
  "2MR": "R",  "AAR": "R",  "AGM": "R",  "ARO": "R",  "CIR": "R",  "DA2": "R",
  "DAR": "R",
  "MEN": "N",  "DSG": "N",  "SNN": "N",  "AHB": "N",  "DMH": "N",
  "0TD": "D",  "BFD": "D",  "BH2": "D",  "BHD": "D",  "DAS": "D",  "IAS": "D",
  "PHD": "D",
  "CAF": "C",  "CAS": "C",  "CCC": "C",  "CME": "C",  "CMT": "C",  "CSD": "C",
  "CSO": "C",  "CSP": "C",  "CSS": "C",  "CSW": "C",  "CSX": "C",  "CSZ": "C",
  "CY3": "C",  "DCY": "C",  "OCS": "C",  "OMC": "C",  "QPA": "C",  "SCH": "C",
  "SCY": "C",  "SMC": "C",  "SNC": "C",  "YCM": "C",
  "DGN": "Q",  "MEQ": "Q",  "MGN": "Q",
  "3GL": "E",  "B3E": "E",  "CGU": "E",  "DGL": "E",  "FGA": "E",  "PCA": "E",
  "CR2": "G",  "GDP": "G",  "GHP": "G",  "GL3": "G",  "TGP": "G",  "OMG": "G",
  "SAR": "G",
  "CR8": "H",  "DDE": "H",  "DHI": "H",  "HIA": "H",  "HIC": "H",  "HIP": "H",
  "HIQ": "H",  "MHS": "H",  "NEP": "H",
  "DIL": "I",  "ILX": "I",  "IML": "I",  "TSI": "I",
  "DLE": "L",  "LED": "L",  "MK8": "L",  "MLE": "L",  "MLL": "L",  "NLE": "L",
  "ALY": "K",  "BTK": "K",  "DLY": "K",  "KCX": "K",  "KPI": "K",  "LLP": "K",
  "LYR": "K",  "LYZ": "K",  "M3L": "K",  "MCL": "K",  "MLY": "K",  "MLZ": "K",
  "PRK": "K",
  "AME": "M",  "CXM": "M",  "MSE": "M",  "FME": "M",  "MED": "M",  "MHO": "M",
  "MME": "M",  "OMT": "M",  "SME": "M",
  "DPN": "F",  "MEA": "F",  "NFA": "F",  "PHI": "F",  "PHL": "F",
  "DPR": "P",  "HY3": "P",  "HYP": "P",
  "CSH": "S",  "FGL": "S",  "GYS": "S",  "SEP": "S",  "OAS": "S",  "OSE": "S",
  "SAC": "S",
  "BMT": "T",  "CRF": "T",  "DBU": "T",  "DTH": "T",  "OLT": "T",  "TH5": "T",
  "TPO": "T",
  "0AF": "W",  "6CW": "W",  "BTR": "W",  "DTR": "W",  "HTR": "W",  "HYD": "W",
  "TOX": "W",  "TRQ": "W",  "TRY": "W",
  "DTY": "Y",  "IYR": "Y",  "NIY": "Y",  "OMY": "Y",  "PTR": "Y",  "TPQ": "Y",
  "TYC": "Y",  "TYI": "Y",  "TYS": "Y",
  "DVA": "V",  "FVA": "V",  "MVA": "V",
  "PSU": "U"
}

multiple_reslist = {
  "XAA": ["X"]+list(reslist.values()),
  "UNK": ["X"]+list(reslist.values()),
  "ASX": ["B", "D", "N"],
  "GLX": ["Z", "E", "Q"],
  "XLE": ["J", "I", "L"],
  "DAL": ["A", "S"],  "MDO": ["A", "Y"],
  "CRQ": ["Q", "Y"],
  "CH6": ["M", "Y"],  "CRK": ["M", "Y"],  "NRQ": ["M", "Y"],
  "DAH": ["F", "Y"],
  "DHA": ["S", "Y"],  "DSN": ["S", "C"],
  "CRO": ["T", "Y"],
  "2KT": ["X", "T"],  "4MM": ["X", "M"],  "ACE": ["X"]+list(reslist.values()),
  "AMP": ["X", "G", "Y"],  "BTN": ["X", "K"],  "DMG": ["X", "G"],
  "FMT": ["X", "K"],  "LAC": ["X", "S"],  "NH2": ["X"]+list(reslist.values()),
  "PBE": ["X", "P"],  "PXU": ["X", "P"],  "PYR": ["X", "Y"],
  "SC2": ["X", "C"],  "SIN": ["X", "C"]
}

nuclist = {
  "A": "A",  "DA": "A",
  "T": "T",  "DT": "T",
  "C": "C",  "DC": "C",
  "G": "G",  "DG": "G",
  "U": "U",  "DU": "U",
  "I": "I",  "DI": "I"
}

mod_nuclist = {
  "1MA": "A",  "2MA": "A",  "6MZ": "A",  "MA6": "A",  "MIA": "A",
  "4OC": "C",  "5CM": "C",  "CBR": "C",  "DOC": "C",
  "YG":  "G",  "1MG": "G",  "2MG": "G",  "7MG": "G",  "8OG": "G",  "M2G": "G",
  "4SU": "U",  "5BU": "U",  "5MU": "U",  "BRU": "U",  "H2U": "U",  "OMU": "U",
  "UR3": "U",
  "C5P": "X",  "U5P": "X"
}

multiple_nuclist = {
  "R": ["R", "A", "G"],
  "Y": ["Y", "C", "T", "U"],
  "K": ["K", "G", "T", "U"],
  "M": ["M", "A", "C"],
  "S": ["S", "C", "G"],
  "W": ["W", "A", "T", "U"],
  "B": ["B", "C", "G", "T", "U"],
  "D": ["D", "A", "G", "T", "U"],
  "H": ["H", "A", "C", "T", "U"],
  "V": ["V", "A", "C", "G"],
  "N": ["N", "A", "C", "G", "T", "U"]
}
