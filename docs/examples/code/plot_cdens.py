"""
Contact Density Plotting
========================

This script contains a simple example of how you can evaluate
the contact density of your contact map using ConKit

"""

import conkit.io
import conkit.plot

# Define the input variables
sequence_file = "4p9g/4p9g.fasta"
sequence_format = "fasta"
contact_file = "4p9g/4p9g.mat"
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

# Truncate your contact list to 10L contacts
conpred = conpred[:int(seq.seq_len * 10.)]

# Then we can plot the density plot
fig = conkit.plot.ContactDensityFigure(conpred, legend=True)
fig.savefig("4p9g/4p9g_cdens.png")
