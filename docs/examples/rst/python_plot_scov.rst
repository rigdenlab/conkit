
Sequence Coverage Plotting
--------------------------

.. code-block:: python

   >>> msa = conkit.io.read('toxd/toxd.a3m', 'a3m')
   >>> conkit.plot.SequenceCoverageFigure(msa)

The following plot will be produced. Your alignment coverage is shown with the black line with each point corresponding to a residue in the alignment. The red and green lines give you indicators of how good your alignment is.

If parts or all of your coverage fall below the red line, the suitability for covariance-based contact prediction is very little. If you parts of the alignment are well above the red line, possibly even the green, and bigger chunks are below, then you might want to consider re-defining your sequence boundaries to predict contacts only for the well-covered area.

If most residues in your alignment have a better coverage than the green line, the alignment is well-suited for covariance-based predictions.

.. image:: ../images/toxd_scov_plot.png
   :alt: Toxd Sequence Coverage Plot
