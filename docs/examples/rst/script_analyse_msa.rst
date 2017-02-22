
Multiple Sequence Alignment Analysis
------------------------------------

.. warning::
   You require the optional dependency `SciPy <http://scipy.org/>`_ package to use this script. If you are unsure if it is installed on your system, refer to the :ref:`Installation` documentation


If you would like to analyse a Multiple Sequence Alignment (MSA) file, you can do so using ConKit's provided script, which is called ``conkit-msatool``.

.. code-block:: bash

   $> conkit-msatool toxd/toxd.a3m a3m

The call above analyses the ``toxd.a3m`` MSA file, which is in ``a3m`` format. This call with will procude the following output:


.. code-block:: none

   Input MSA File:                   toxd/toxd.a3m
   Input MSA Format:                 a3m
   Sequence Identity Threshold:      0.7
   Length of the Target Sequence:    59
   Total Number of Sequences:        13448
   Number of Effective Sequences:    3318
   Sequence Coverage Plot:           toxd/toxd.png

The output tells you the file and format you have provided. It also prints which identity threshold was used to compare sequences during the analysis. Furthermore, it tells you the total number of sequences in your alignment and the number of effective sequences, or depth, of your alignment. Finally, this script will produce a plot that illustrates the coverage of your alignment in individual positions. The plot is shown below:

.. _Toxd Frequency Plot:

.. image:: ../images/toxd_scov_plot.png
   :alt: Toxd Sequence Coverage Plot
   :scale: 30
   :align: center

