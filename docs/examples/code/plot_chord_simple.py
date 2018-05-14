"""
Simple contact cmap Chord plotting
=================================

This script contains a simple example of how you can plot
contact cmaps in Chord Diagram style using ConKit

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
conpred.set_sequence_register()

# We need to tidy our contact prediction before plotting
conpred.remove_neighbors(inplace=True)
conpred.sort('raw_score', reverse=True, inplace=True)

# Finally, we don't want to plot all contacts but only the top-L,
# so we need to slice the contact cmap
cmap = conpred[:conpred.sequence.seq_len]

# Then we can plot the cmap
fig = conkit.plot.ContactMapChordFigure(cmap, legend=True)
fig.savefig("toxd/toxd.png")
