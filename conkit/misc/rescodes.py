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

#list obtained form nakb.org (14/04/26), likely needs updating somewhat frequently
mod_nuclist = {
 '02I': 'A', '05A': 'T', '05H': 'T', '05K': 'T', '08Q': 'C', '0A': 'A', '0AD': 'G', '0C': 'C', '0DA': 'A', '0DC': 'C', '0DG': 'G', '0DT': 'T', '0G': 'G', '0KZ': 'T', '0R8': 'C', '0U': 'U', '0U1': 'T',
 '10C': 'C', '125': 'U', '126': 'U', '127': 'U', '12A': 'A', '16B': 'C', '18M': 'G', '18Q': 'T', '1AP': 'A', '1CC': 'C', '1DP': 'A', '1FC': 'C', '1MA': 'A', '1MG': 'G', '1RN': 'U', '1SC': 'C', '1TL': 'T', '1TW': 'G', '1W5': 'C', '1WA': 'G',
 '23G': 'G', '2AD': 'A', '2AR': 'A', '2AT': 'T', '2AU': 'U', '2BD': 'G', '2BT': 'T', '2BU': 'A', '2DA': 'A', '2DF': 'C', '2DM': '?', '2DT': 'T', '2EG': 'G', '2FE': 'A', '2FI': 'G', '2GF': 'C', '2GT': 'T', '2IA': 'A', '2JU': 'T', '2JV': 'G', '2L8': 'T', '2LA': 'G', '2LF': 'G', '2MA': 'A', '2MG': 'G', '2MU': 'U', '2NT': 'T', '2OP': '?', '2OT': 'T', '2PR': 'G', '2SG': 'G', '2ST': 'T', '2YR': 'C',
 '30U': '?', '31H': 'A', '31M': 'A', '365': 'A', '3AU': 'U', '3D1': 'A', '3DA': 'A', '3DR': 'A', '3KA': 'T', '3ME': 'U', '3PO': '?', '3TD': 'U', '3ZO': 'G',
 '40A': 'A', '40C': 'C', '40G': 'G', '40T': 'T', '45A': 'A', '47C': 'C', '48Z': 'A', '4AC': 'C', '4DG': 'G', '4DU': 'A', '4E9': 'G', '4EN': 'A', '4JA': '?', '4MF': 'A', '4OC': 'C', '4PC': 'C', '4SC': 'C', '4SU': 'U', '4U3': 'C',
 '50L': 'G', '50N': 'C', '56B': 'G', '574': 'A', '5AA': 'A', '5AT': 'T', '5BT': 'C', '5BU': 'U', '5CF': 'C', '5CG': 'G', '5CM': 'C', '5CY': '?', '5DB': 'T', '5EJ': 'T', '5FC': 'C', '5FU': 'U', '5GP': 'G', '5HC': 'C', '5HM': 'C', '5HT': 'T', '5HU': 'T', '5IC': 'C', '5IU': 'T', '5JO': 'A', '5MC': 'C', '5MU': 'U', '5NC': 'C', '5OC': 'C', '5PC': 'C', '5PY': 'T', '5SE': 'T', '5UA': 'A', '5UD': 'U',
 '61H': 'C', '63G': 'G', '63H': 'G', '63T': 'A', '64P': 'T', '64T': 'T', '68Z': 'G', '6F7': 'U', '6FC': 'C', '6FK': 'G', '6FM': 'A', '6FU': 'T', '6HA': 'A', '6HB': 'A', '6HC': 'C', '6HG': 'G', '6HT': 'T', '6IA': 'A', '6MA': 'A', '6MD': 'A', '6MI': 'G', '6MZ': 'A', '6NW': 'A', '6OG': 'G', '6OO': 'C', '6OP': 'U', '6PO': 'G', '6TW': 'G', '6U0': 'A',
 '70U': 'U', '73W': 'C', '75B': 'U', '77Y': 'T', '7AT': 'A', '7BG': 'G', '7DA': 'A', '7GU': 'G', '7MG': 'G', '7OK': 'C', '7S3': 'G', '7SN': 'G', '7TE': 'A',
 '84E': 'T', '85Y': 'U', '8AA': 'G', '8AF': 'G', '8AG': 'G', '8AH': 'A', '8AN': 'A', '8AZ': 'G', '8BA': 'A', '8DT': 'T', '8FG': 'G', '8GM': 'G', '8MG': 'G', '8NI': 'G', '8OG': 'G', '8OS': 'G', '8PI': 'G', '8PY': 'G', '8RJ': 'U', '8RO': 'C', '8XA': 'A', '8XC': 'C', '8XG': 'G', '8XU': 'U', '8Y9': 'T', '8YN': 'C',
 '91N': 'A', '92F': 'T', '93D': 'C', '94O': 'T', '96T': '?', '9O4': 'A', '9QV': 'U', '9SI': 'A', '9SY': 'A', '9V9': 'U',
 'A1A0L': 'U', 'A1AAZ': 'T', 'A1B8A': 'A', 'A1BBA': 'C', 'A1EFN': 'A', 'A1ELZ': 'A', 'A1H3G': 'A', 'A1I9V': 'U', 'A1IC0': 'A', 'A1IC1': 'U', 'A1ID5': 'A', 'A1IEA': 'A', 'A1L3P': 'U', 'A1L82': 'C', 'A1L89': 'A', 'A1LXS': 'A', 'A1LZ3': 'U', 'A1MA9': 'G',
 'A23': 'A', 'A2L': 'A', 'A2M': 'A', 'A2P': 'A', 'A38': 'A', 'A3A': 'A', 'A3P': 'A', 'A40': 'A', 'A43': 'A', 'A44': 'A', 'A47': 'A', 'A5L': 'A', 'A5M': 'C', 'A5O': 'A', 'A66': 'A', 'A7C': 'A', 'A7E': 'A', 'A9Z': 'A', 'AAB': 'A', 'ABR': 'A', 'ABS': 'A', 'ACA': '?', 'ACE': '?', 'AD2': 'A', 'ADN': 'A', 'ADP': 'A', 'ADS': 'A', 'AET': 'A', 'AF2': 'A', 'AFG': 'G', 'AG9': 'C', 'AGD': 'G', 'AI5': 'C', 'AMP': 'A', 'ANZ': 'A', 'AP7': 'A', 'APC': 'A', 'APN': 'A', 'AS': 'A', 'ASU': 'A', 'AT7': 'A', 'ATD': 'T', 'ATL': 'T', 'ATM': 'T', 'ATP': 'A', 'AVC': 'A', 'AWC': 'T', 'AZW': '?',
 'B4P': 'A', 'B7C': 'C', 'B86': 'C', 'B8H': 'U', 'B8K': 'G', 'B8N': 'U', 'B8Q': 'C', 'B8T': 'C', 'B8W': 'G', 'B9B': 'G', 'B9H': 'C', 'BGH': 'G', 'BGM': 'G', 'BGR': 'G', 'BMN': 'C', 'BOE': 'T', 'BRU': 'T', 'BTN': '?', 'BZG': 'G',
 'C2L': 'C', 'C2S': 'C', 'C31': 'C', 'C34': 'C', 'C36': 'C', 'C37': 'C', 'C38': 'C', 'C42': 'C', 'C43': 'C', 'C45': 'C', 'C46': 'C', 'C49': 'C', 'C4J': 'U', 'C4S': 'C', 'C5L': 'C', 'C5P': 'C', 'C66': 'C', 'C6G': 'G', 'C7R': 'C', 'C7S': 'C', 'CAR': 'C', 'CBR': 'C', 'CBV': 'C', 'CCC': 'C', 'CDW': 'C', 'CFL': 'C', 'CFV': 'C', 'CFZ': 'C', 'CG1': 'G', 'CGY': 'C', 'CH': 'C', 'CJ1': 'G', 'CM0': 'U', 'CMR': 'C', 'CNV': '?', 'CP1': 'C', 'CPN': 'C', 'CSL': 'C', 'CSM': 'T', 'CTG': 'T', 'CTP': 'C', 'CUD': 'C', 'CVC': 'C', 'CX2': 'C',
 'D00': 'C', 'D3': 'A', 'D33': 'A', 'D3N': 'C', 'D4B': 'C', 'D5M': 'A', 'DCM': 'C', 'DCZ': 'C', 'DDG': 'G', 'DDN': 'T', 'DDX': 'A', 'DFT': 'T', 'DG8': 'G', 'DI': 'G', 'DJF': 'A', 'DLY': '?', 'DN': 'A', 'DNR': 'C', 'DOC': 'C', 'DOP': '?', 'DP': 'G', 'DPY': 'C', 'DRP': 'C', 'DRZ': '?', 'DU': 'T', 'DUR': 'U', 'DUZ': 'T', 'DV3': 'A', 'DX': 'G', 'DXD': 'A', 'DZ': 'C', 'DZM': 'A',
 'E': 'A', 'E1X': 'A', 'E3C': 'C', 'E6G': 'G', 'E7G': 'G', 'EAN': 'T', 'EDA': 'A', 'EDC': 'C', 'EDI': 'G', 'EFG': 'G', 'EHG': 'G', 'EIT': 'T', 'EIX': 'C', 'EQ0': 'G', 'EW3': 'T', 'EWC': 'G', 'EXC': 'C',
 'F2T': 'U', 'F3H': 'T', 'F3N': 'A', 'F4H': 'T', 'F4Q': 'G', 'F5H': 'T', 'F6H': 'T', 'F6U': 'A', 'F6X': 'C', 'F73': 'G', 'F74': 'G', 'F7H': 'C', 'F7K': 'G', 'F7O': 'A', 'F7R': 'C', 'F7U': 'G', 'F7X': 'G', 'F86': 'A', 'FA2': 'A', 'FAG': 'G', 'FAX': 'A', 'FDG': 'G', 'FFD': 'C', 'FHU': 'U', 'FME': '?', 'FMG': 'G', 'FMU': 'U', 'FOX': 'G',
 'G2L': 'G', 'G2M': 'G', 'G2S': 'G', 'G31': 'G', 'G35': 'G', 'G36': 'G', 'G38': 'G', 'G3A': 'G', 'G46': 'G', 'G47': 'G', 'G48': 'G', 'G49': 'G', 'G5J': 'G', 'G7M': 'G', 'GAO': 'G', 'GCK': 'C', 'GCP': 'G', 'GDO': 'G', 'GDP': 'G', 'GF0': 'G', 'GF2': 'G', 'GFL': 'G', 'GMP': 'G', 'GMS': 'G', 'GMU': 'U', 'GMX': 'G', 'GN7': 'G', 'GNE': 'G', 'GNG': 'G', 'GP3': 'G', 'GPN': 'G', 'GRB': 'G', 'GS': 'G', 'GSR': 'G', 'GSS': 'G', 'GT3': 'G', 'GTA': 'G', 'GTG': 'G', 'GTP': 'G', 'GX1': 'G',
 'H2U': 'U', 'HCX': '?', 'HEU': 'T', 'HFA': '?', 'HGL': 'G', 'HHU': 'T', 'HHX': 'T', 'HN0': 'G', 'HN1': 'G', 'HOB': 'A', 'HOL': 'A', 'HPD': '?', 'HYJ': 'G',
 'I': 'G', 'I2T': 'U', 'I4U': 'U', 'IC': 'C', 'IG': 'G', 'IGU': 'G', 'IKS': 'C', 'ILK': 'U', 'IMC': 'C', 'IMP': 'G', 'IPN': 'U', 'IQG': 'G', 'IU': 'U',
 'J0X': 'C', 'J4T': 'T', 'JDT': 'T', 'JMC': 'C', 'JMH': 'C', 'JSP': 'T',
 'K12': '?', 'K1F': 'C', 'K2F': 'A', 'K39': 'G', 'KAG': 'G', 'KAK': 'G', 'KBC': 'T', 'KGV': 'A', 'KPN': '?',
 'L1J': 'G', 'L2B': 'U', 'L3X': 'A', 'L5R': 'T', 'L8P': 'C', 'LCA': 'A', 'LCC': 'C', 'LCG': 'G', 'LDG': 'G', 'LGP': 'G', 'LHC': 'C', 'LHH': 'C', 'LHO': 'C', 'LKC': 'C', 'LR6': 'T', 'LSH': 'T', 'LST': 'T', 'LV2': 'C', 'LVR': '?', 'LWM': 'G',
 'M1G': 'G', 'M1Y': 'U', 'M2G': 'G', 'M3X': 'C', 'M5M': 'C', 'M7A': 'A', 'M7G': 'G', 'M7M': 'G', 'MA6': 'A', 'MA7': 'A', 'MBZ': 'A', 'MCY': 'C', 'MDJ': 'C', 'MDK': 'C', 'MDQ': 'C', 'MDU': 'T', 'MDV': 'A', 'ME6': 'C', 'MF7': 'G', 'MFO': 'G', 'MFT': 'T', 'MG1': 'G', 'MGT': 'G', 'MHG': 'G', 'MIA': 'A', 'MKX': 'A', 'MM7': 'C', 'MMT': 'T', 'MMX': 'C', 'MNU': 'U', 'MRG': 'G', 'MTR': 'T', 'MTU': 'A', 'MUM': 'U',
 'N2G': 'G', 'N4S': 'C', 'N5I': 'A', 'N5M': 'C', 'N68': '?', 'N6G': 'G', 'N79': 'A', 'N7X': 'C', 'NCU': 'C', 'NCX': 'A', 'NF2': 'U', 'NH2': '?', 'NME': '?', 'NMS': 'T', 'NMT': 'T', 'NP3': 'G', 'NR0': 'A', 'NR1': 'A', 'NRI': 'A', 'NRL': 'C', 'NSF': 'U', 'NSU': 'T', 'NTT': 'T', 'NYM': 'T',
 'O2C': 'C', 'O2G': 'G', 'O2Z': 'A', 'OBX': 'C', 'OFC': 'C', 'OGX': 'G', 'OHU': 'T', 'OIQ': 'T', 'OKN': 'C', 'OKQ': 'C', 'OKT': 'T', 'OMC': 'C', 'OMG': 'G', 'OMU': 'U', 'ONE': 'U', 'OPN': '?', 'ORP': 'A', 'OWR': 'A',
 'P': 'G', 'P2T': 'T', 'P2U': 'T', 'P4U': 'U', 'P5P': 'A', 'P7G': 'G', 'P9G': 'G', 'PA9': '?', 'PAE': '?', 'PBT': 'T', 'PDI': '?', 'PDU': 'T', 'PE6': '?', 'PED': '?', 'PGE': '?', 'PGN': 'G', 'PGP': 'G', 'PHA': '?', 'PO2': '?', 'PO4': '?', 'POP': '?', 'PPS': 'A', 'PPU': 'A', 'PPW': 'G', 'PQ1': 'G', 'PRN': 'A', 'PST': 'T', 'PSU': 'U', 'PU': 'A', 'PVX': 'C', 'PYO': 'U', 'PYP': 'C', 'PYY': 'C',
 'Q61': 'G', 'QBT': 'T', 'QCE': 'T', 'QCK': 'T', 'QGJ': 'T', 'QRV': 'G', 'QSK': 'C', 'QSQ': 'A', 'QUO': 'G',
 'R': 'A', 'RBD': 'A', 'RCE': 'T', 'RDG': 'G', 'RF5': 'A', 'RFJ': 'G', 'RIA': 'A', 'RMP': 'A', 'RP5': '?', 'RPC': 'C', 'RSP': 'C', 'RSQ': 'C', 'RUS': 'U', 'RY': 'C',
 'S02': '?', 'S2M': 'T', 'S4A': 'A', 'S4C': 'C', 'S4G': 'G', 'S6G': 'G', 'S6M': 'T', 'S8U': 'G', 'S9L': '?', 'SAY': 'T', 'SC': 'C', 'SDE': 'A', 'SDG': 'G', 'SDH': 'G', 'SJO': 'G', 'SMP': 'A', 'SMT': 'T', 'SOS': 'G', 'SPT': 'T', 'SRA': 'A', 'SSU': 'U', 'SUR': 'U',
 'T0N': 'T', 'T0P': 'G', 'T0T': 'C', 'T23': 'U', 'T2S': 'T', 'T2T': 'U', 'T32': 'T', 'T38': 'T', 'T39': 'T', 'T3P': 'T', 'T41': 'T', 'T48': 'T', 'T49': 'T', 'T4S': 'T', 'T5O': 'T', 'T5S': 'T', 'T64': 'T', 'T66': 'T', 'T6A': 'A', 'TA3': 'T', 'TAF': 'T', 'TC': 'C', 'TC1': 'C', 'TCJ': 'C', 'TCP': 'T', 'TCY': 'A', 'TDY': 'T', 'TED': 'T', 'TFE': 'T', 'TFO': 'A', 'TFT': 'T', 'TG': 'G', 'TGP': 'G', 'THM': 'T', 'THP': 'T', 'THX': 'T', 'TJU': 'C', 'TKW': 'C', 'TLB': 'U', 'TLC': 'T', 'TLN': 'U', 'TM2': 'U', 'TMS': '?', 'TP1': 'T', 'TPC': 'C', 'TPN': 'T', 'TS6': '?', 'TSE': '?', 'TSP': 'T', 'TT': 'T', 'TTD': 'T', 'TTI': 'T', 'TTM': 'T', 'TTP': 'T', 'TX2': 'C', 'TYD': 'T',
 'U23': 'U', 'U2L': 'U', 'U2M': '?', 'U31': 'U', 'U33': 'T', 'U34': 'U', 'U36': 'U', 'U37': 'U', 'U48': 'C', 'U4M': 'T', 'U5M': 'T', 'U5P': 'U', 'U5R': 'T', 'U7B': 'C', 'U7T': 'A', 'U8U': 'U', 'UAR': 'U', 'UBB': 'U', 'UBD': 'U', 'UBI': 'T', 'UCL': 'T', 'UD5': 'U', 'UEL': 'G', 'UF2': 'T', 'UFB': 'T', 'UFP': 'T', 'UFR': 'T', 'UFT': 'T', 'ULF': 'T', 'UMC': 'T', 'UMO': 'U', 'UMP': 'T', 'UMS': 'T', 'UMX': 'T', 'UNK': '?', 'UOA': 'U', 'UOB': 'U', 'UPE': 'T', 'UPS': 'T', 'UPV': 'U', 'UR3': 'U', 'URT': 'A', 'URU': 'U', 'URX': 'T', 'US1': 'T', 'US2': 'T', 'US3': 'T', 'US4': 'T', 'US5': 'U', 'USM': 'T', 'UTB': 'T', 'UTP': 'U', 'UVP': 'U', 'UVX': 'T', 'UWJ': 'T', 'UY1': 'U', 'UY4': 'A', 'UZL': 'C', 'UZR': 'U',
 'VC7': 'G', 'VET': 'T', 'VKJ': 'G', 'VM6': '?',
 'WC7': 'C', 'WUH': 'T', 'WVQ': 'G',
 'X0K': '?', 'X4A': 'A', 'X6H': 'G', 'XAD': 'A', 'XAE': 'A', 'XAL': 'A', 'XAN': 'G', 'XAR': 'A', 'XB9': 'T', 'XC': 'C', 'XCI': 'T', 'XCL': 'C', 'XCR': 'C', 'XCS': 'C', 'XCT': 'C', 'XCY': 'C', 'XDD': 'A', 'XDJ': 'G', 'XDV': '?', 'XDY': 'T', 'XE6': 'C', 'XEC': 'U', 'XFC': 'C', 'XGA': 'G', 'XGL': 'G', 'XGR': 'G', 'XGU': 'G', 'XMP': 'G', 'XNY': 'A', 'XPB': 'G', 'XSX': 'U', 'XTF': 'T', 'XTH': 'T', 'XTL': 'T', 'XTR': 'T', 'XTY': 'T', 'XUA': 'A', 'XUG': 'G', 'XY7': 'T',
 'Y': 'A', 'Y5P': 'U', 'YA4': 'A', 'YB9': '?', 'YCO': 'C', 'YG': 'G', 'YPE': '?', 'YPF': '?', 'YQS': 'T', 'YRR': 'A', 'YTI': 'T', 'YYG': 'G',
 'Z': 'C', 'ZAD': 'A', 'ZBC': 'C', 'ZBU': 'U', 'ZCY': 'C', 'ZDU': 'T', 'ZGU': 'G', 'ZHP': 'U', 'ZIV': 'G', 'ZJS': 'A', 'ZTH': 'T'
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

IRRELEVANT_res_codes = { 
  'HOH', 
  'MG',
  ''}
