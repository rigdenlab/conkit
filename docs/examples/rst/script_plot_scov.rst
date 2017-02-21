
Sequence Coverage Plotting
--------------------------

The ``scov`` subcommand of the ``conkit-plot`` script is used to plot the coverage plot of the multiple sequence alignment.

.. code-block:: bash

   $> conkit-plot scov toxd/toxd.a3m a3m

The following plot will be produced. Your alignment coverage is shown with the black line with each point corresponding to a residue in the alignment. The red and green lines give you indicators of how good your alignment is.

.. image:: ../images/toxd_scov_plot.png
   :alt: Toxd Sequence Coverage Plot

If parts or all of your coverage fall below the red line, the suitability for covariance-based contact prediction is very little. If you parts of the alignment are well above the red line, possibly even the green, and bigger chunks are below, then you might want to consider re-defining your sequence boundaries to predict contacts only for the well-covered area.

If most residues in your alignment have a better coverage than the green line, the alignment is well-suited for covariance-based predictions.

