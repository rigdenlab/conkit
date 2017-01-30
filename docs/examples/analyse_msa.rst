.. _example_analyse_msa:

Multiple Sequence Alignment Analysis
====================================

.. warning::
   You require the `SciPy <http://scipy.org/>`_ package to use this functionality. If you are unsure if it is installed on your system, refer to the :ref:`installation` documentation


.. note::
   ConKit provides a script that allows you to plot contact maps from the command line. The script is called ``conkit.msatool`` and installed automatically.


To start, you need to read the file into its hierarchy.

.. code-block:: python

   >>> import conkit
   >>> msa = conkit.io.read('pathto/msa.file', 'msa_format')

With the code above, we created a :obj:`conkit.core.SequenceFile` hierarchy. This allows you to store a multiple sequence file in a single hierarchy.

Your alignment stored in this hierarchy allows you to extract some information immediately without any computation. For example,

1. The length of your input sequence:

.. code-block:: python

   >>> print(msa.top_sequence.seq_len)
   100

2. The total number of sequences in the alignment:

.. code-block:: python

   >>> print(msa.nseqs)
   12345

If you would like to analyse this further and determine the usefulness of the alignment for covariance-based contact prediction, you might want to compute the alignment depth, sometimes also referred to as the number of effective sequences. In ConKit, you can do this by running the following command:

.. code-block:: python

   >>> print(msa.calculate_meff())
   255

This number indicates that out of the 12345 total sequences 255 are effective and will contribute significantly to your covariance-based contact prediction.