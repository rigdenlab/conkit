"""
Simple contact cmap plotting 1
=============================

This script contains a simple example of how you can plot
contact cmaps using ConKit

"""

import conkit.io
import conkit.plot

# Define the input variables
sequence_file = "toxd/toxd.fasta"
sequence_format = "fasta"
contact_file = "toxd/toxd.mat"
contact_format = "ccmpred"

# Create ConKit hierarchies
#       Note, we only need the first Sequence/ContactMap
#       from each file
seq = conkit.io.read(sequence_file, sequence_format).top
conpred = conkit.io.read(contact_file, contact_format).top

# Assign the sequence register to your contact prediction
conpred.sequence = seq
conpred.assign_sequence_register()

# Then we can plot the cmap
fig = conkit.plot.ContactMapMatrixFigure(cmap)
fig.savefig("toxd/toxd.png")
