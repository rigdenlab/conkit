
Contact Map Chord Diagram Plotting
----------------------------------

If you would like to plot a contact map in Chord diagram style using ConKit without the overhead of using Python, you can simply use the ``conkit.plot`` script.

.. code-block:: bash

   $> conkit.plot chord toxd/toxd.fasta fasta toxd/toxd.mat ccmpred

The call above uses the contact prediction file ``toxd.mat`` file, which is in ``ccmpred`` format, and plots the following contact map stored in the file ``toxd/toxd.png``

.. image:: ../images/toxd_chord_simple.png
   :alt: Toxd Chord Simple
   :width: 500px
