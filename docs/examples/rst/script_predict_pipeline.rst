
Contact Prediction
------------------

.. warning::

   External software is required to execute this script. For further information, refer to the :ref:`installation` page.

The script to run contact prediction using `HHblits <https://github.com/soedinglab/hh-suite>`_ to generate your sequence alignment and `CCMpred <https://github.com/soedinglab/CCMpred>`_ to predict the contacts is called ``conkit-predict``. You can run this script using two different modes.

**1. Starting with a sequence**

.. code-block:: bash

   $> conkit-predict seq <path/to/hhblits> <path/to/hhblits_database> <path/to/ccmpred> toxd/toxd.fasta fasta

The call above uses your sequence file ``toxd/toxd.fasta`` in ``fasta`` format to first generate a Multiple Sequence Alignment. It will then analyse your alignment identical to the ``conkit-msatool`` script. It will also sort out all the required conversions before executing CCMpred to run the contact prediction. Finally, it will analyse your contact prediciton and plot a contact map, just like the ``conkit-plot`` script does.

**2. Starting with a Multiple Sequence Alignment**

.. code-block:: bash

   $> conkit-predict aln <path/to/ccmpred> toxd/toxd.a3m a3m

This call performs identical operations to the full call under point 1, except that it skips the generation of the alignment. This might be particularly useful if you have limited disk space and cannot store the rather large sequence database that HHblits requires. You can generate your alignment using online servers, such as the `HHblits Server <https://toolkit.tuebingen.mpg.de/hhblits>`_ or the `Jackhmmer Server <https://www.ebi.ac.uk/Tools/hmmer/search/jackhmmer>`_. Both formats are accepted by ConKit, the keywords can be found in the :ref:`file_formats`.
