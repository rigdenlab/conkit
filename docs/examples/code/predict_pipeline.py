"""
Simple contact prediction pipeline
==================================

This script contains a simple example of how we can create pipelines
using ConKit.

.. warning::
   You need to exchange the paths to the executables

"""

import conkit

# Define the input variables
sequence_file = "toxd/toxd.fasta"
sequence_format = "fasta"

# Define the paths to the software we use
hhblits_exe = "path/to/hhblits"                 # <-- MODIFY THIS
hhblits_database = "path/to/hhblits_database"   # <-- MODIFY THIS
ccmpred_exe = "path/to/ccmpred"                 # <-- MODIFY THIS

# Generate a Multiple Sequence Alignment
print("Generating the Multiple Sequence Alignment")
a3m_file = "toxd/toxd.a3m"
hhblits_cline = conkit.applications.HHblitsCommandLine(
    cmd=hhblits_exe, database=hhblits_database, input=sequence_file, oa3m=a3m_file
)
hhblits_cline()  # Execute HHblits

# Analyse the alignment
msa = conkit.io.read(a3m_file, "a3m")
print("Length of the Target Sequence:    %d" % msa.top_sequence.seq_len)
print("Total Number of Sequences:        %d" % msa.nseqs)
print("Number of Effective Sequences:    %d" % msa.calculate_meff())

# Plot the amino acid coverage per position in the alignment
seq_cov_file = "toxd/toxd.freq.png"
conkit.plot.SequenceCoverageFigure(msa, file_name=seq_cov_file)
print("Sequence Coverage Plot:           %s" % seq_cov_file)

# Convert the alignment into a CCMpred-readable format
jones_file = "toxd/toxd.jones"
conkit.io.write(jones_file, "jones", msa)

# Predict the contacts
print("Predicting contacts")
mat_file = "toxd/toxd.mat"
ccmpred_cline = conkit.applications.CCMpredCommandLine(
    alnfile=jones_file, matfile=mat_file
)
ccmpred_cline()  # Execute CCMpred

# Plot the top-30 contacts
conpred = conkit.io.read(mat_file,  "ccmpred").top_map
conpred.remove_neighbors(inplace=True)                  # Remove contacts of neigbouring residues
conpred.sort("raw_score", reverse=True, inplace=True)   # Sort the list of contacts by their score
conpred = conpred[:30]                                  # Slice the contact map
contact_map_file = "toxd/toxd.map.png"
conkit.plot.ContactMapFigure(conpred, file_name=contact_map_file)
print("Contact Map Plot:                 %s" % contact_map_file)

# Convert the contact prediction to a standardised format
casp_file = "toxd/toxd.rr"
conkit.io.convert(mat_file, "ccmpred", casp_file, "casprr")
print("Final Contact Prediction File:    %s" % casp_file)
