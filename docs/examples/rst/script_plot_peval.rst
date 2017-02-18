
Precision Evaluation Plotting
-----------------------------

.. code-block:: bash

   $> conkit.plot peval -j 0.1 -min 0 -max 5 toxd/toxd.pdb pdb toxd/toxd.fasta fasta toxd/toxd.mat ccmpred

Three command line flags above are important to note. The ``-min`` flag is the minimum factor for contact selection, i.e. ``L * min`` number of contacts. The ``-max`` flag is the maximum factor for contact selection, i.e. ``L * max`` number of contacts. The ``-j`` flag defines the stepwise increase inbetween ``-min`` and ``-max``.

.. image:: ../images/toxd_peval_plot.png
   :alt: Toxd Precision Evaluation Plot
