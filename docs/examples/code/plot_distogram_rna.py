"""
RNA distogram plotting
======================

This script shows how to plot an RNA C1' inter-nucleotide distogram
using ConKit with an AlphaFold 3 distogram prediction.

AlphaFold 3 outputs a pre-computed distance probability matrix (NPZ
format), so no atom-type or distance-cutoff parameters are needed at
read time — those are only relevant when computing distances from
structure coordinates.

"""

import conkit.io
import conkit.plot

# Define the input variables
sequence_file = "rna/target.fasta"
sequence_format = "fasta"
distance_file = "rna/target_distogram.npz"
distance_format = "af3npz"

# Read the sequence
seq = conkit.io.read(sequence_file, sequence_format).top

# Read the AF3 distogram
distpred = conkit.io.read(distance_file, distance_format).top

# Assign the sequence register
distpred.sequence = seq
distpred.set_sequence_register()

# Plot the distogram
fig = conkit.plot.DistogramHeatmapFigure(distpred)
fig.savefig("rna/target_distogram.png")
