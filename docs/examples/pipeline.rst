.. _example_pipeline_creation:

Pipeline Design
===============

ConKit can be a very useful tool to create custom prediction pipelines without having to code the tedious file conversions, etc.

To illustrate the simplicity, let's create one right now.

0. Required dependencies
^^^^^^^^^^^^^^^^^^^^^^^^

Before we are getting started, please ensure that you have the following dependencies installed:

* `HH-suite Repo`_
* `HHblits database`_
* `CCMpred Repo`_

1. Importing ConKit modules
^^^^^^^^^^^^^^^^^^^^^^^^^^^

At the top of our pipeline script, we need to import ConKit Python modules to be able to run a contact prediction

To handle file conversions, we need the ConKit I/O package.

.. code-block:: python

   >>> from conkit import io

We also need to import the command line application wrappers so we can execute our applications from within Python.

.. code-block:: python

   >>> from conkit import applications

2. Creating I/O filenames
^^^^^^^^^^^^^^^^^^^^^^^^^

The next step is to define some filenames and variables so we do not get confused as we go along.

.. code-block:: python

   >>> name = 'conkit_example'
   >>> a3m_fname = name + '.a3m'
   >>> a3m_filtered_fname = name + '.filtered.a3m'
   >>> jones_filtered_fname = name + '.filtered.jones'
   >>> matrix_fname = name + '.mat'
   >>> casprr_fname = name + '.rr'
   >>> hhblits_db = 'path_to_hhblits_db'

The files defined are all the files we need in a very simplified contact prediction pipeline. I opted for starting all my filenames with ``conkit_example`` so we could identify them in a single directory.

Now, we also need to provide a FASTA file. Here, my arbitrary FASTA file is called ``conkit.fasta``. Let's define it

.. code-block:: python

   >>> fasta = 'conkit.fasta'

3. Generating a Multiple Sequence Alignment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Above, we defined all filenames required in this pipeline, so let us continue and generate our first Multiple Sequence Alignment.

.. code-block:: python

   >>> hhblits_cline = HHblitsCommandLine(
   ...     input=seq_fname, database=hhblits_db, niterations=2, evalue=0.001, oa3m=a3m_fname
   ... )

The above code looks more complex than the previous ones. However, this design gives us great flexibility when it comes to trialing different options.

First, we create a new :obj:`conkit.applications.HHblitsCommandLine` instance and provide certain options, which instruct HHblits on which parameters to use during execution. Most of them are pretty self-explanatory, but as you can see the filenames we created earlier are now provided to this instance.

.. note::

   Almost all HHblits command line flags are accessible through this interface. For a full list, consult the HHblits documentation.

At this point we have not actually generated a Multiple Sequence Alignment. All we have done is instructed HHblits on how to run our query. To validate the command that will be executed, you can print the full statement.

.. code-block:: python

   >>> print(hhblits_cline)
   hhblits -n 2 -evalue 0.001 -oa3m conkit_example.a3m -i conkit.fasta -d path_to_hhblits_db

Finally, to invoke HHblits, run the following command.

.. code-block:: python

   >>> hhblits_cline()

Great, now we have generated a sequence alignment, which is stored in the filename defined in the ``a3m_fname`` variable.

4. Filtering the MSA to remove redundant sequences
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Next, we would like to remove some of the redundancy in the sequence alignment file to reduce the prediction bias. Luckily, the HH-suite provides us with a filtering algorithm that does exactly that. Let us create a command line instance just like before, only using the :obj:`conkit.applications.HHfilterCommandLine` class this time.

.. code-block:: python

   >>> hhfilter_cline = HHfilterCommandLine(
   ...     input=a3m_fname, output=a3m_filtered_fname, pairwise_identity=90
   ... )

The above command follows the identical style as the :obj:`conkit.applications.HHblitsCommandLine` wrapper. Note, we have defined our pairwise sequence identity to be 90%, i.e. all sequences with a higher sequence identity will be removed. Finally, let us invoke this command and filter our alignment.

.. code-block:: python

   >>> hhfilter_cline()

5. Converting the sequence alignment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This command results in a final, filtered alignment which is nearly ready to be subjected to CCMpred for contact prediction. However, CCMpred does not like the A3M format the HHblits and HHfilter produced. Thus, we need to convert it into a format that CCMpred recognises. This is the point where ConKit main functionality comes in, the conversion of files.

First, we need to create file handlers for the input and output files:

.. code-block:: python

   >>> f_in = open(a3m_filtered_fname, 'r')
   >>> f_out = open(jones_filtered_fname, 'w')

Once these files are open, we can parse them to the ConKit I/O package for conversion, whereby we need to specify the input format, here ``a3m`` and output format ``jones``. For a full list of file formats available, head over to the :ref:`file_formats`.

.. code-block:: python

   >>> io.convert(f_in, 'a3m', f_out, 'jones')

6. Predicting contacts
^^^^^^^^^^^^^^^^^^^^^^

Finally, we can predict contacts using our generated alignment file. To do this, we use CCMpred [note: the syntax is always the same for command line applications].

.. code-block:: python

   >>> ccmpred_cline = CCMpredCommandLine(
   ...     alnfile=jones_filtered_fname, matfile=matrix_fname, renormalize=True
   ... )
   >>> ccmpred_cline()

Our final contact prediction matrix is now stored in the file with the name ``conkit_example.mat``. Again, this format is not really human-readable and you might want to convert it to a more standardised format, e.g.

.. code-block:: python

   >>> f_in = open(matrix_fname, 'r')
   >>> f_out = open(casprr_fname, 'w')
   >>> io.convert(f_in, 'ccmpred', f_out, 'casprr')

This will produce your final contact prediction in Casp RR format in the file ``conkit_example.rr``

.. note::

   Did you notice that the function call for converting files is identical for sequence and contact files?


.. _CCMpred Repo: https://github.com/soedinglab/ccmpred
.. _HH-suite Repo: https://github.com/soedinglab/hh-suite
.. _HHblits database: http://wwwuser.gwdg.de/%7Ecompbiol/data/hhsuite/databases/hhsuite_dbs/
