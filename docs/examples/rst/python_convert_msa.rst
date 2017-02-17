
File Format Conversion
----------------------

In order to convert files in ConKit, we need to use the ConKit I/O framework.

.. note::
   ConKit I/O framework consists of three main functions that handle the relevant parsers: :func:`read() <conkit.io.read>`, :func:`write() <conkit.io.write>` and :func:`convert() <conkit.io.convert>`. The latter effectively uses the former two but handles everything in one step.

**1. Files can be read in ConKit's internal hierarchies using simple Python code.**

.. code-block:: python
  
   >>> msa = conkit.io.read('toxd/toxd.a3m', 'a3m')

**2. Sequence hierarchies can also be written in a similarly easy format. Using the ``msa`` hierarchy we have created above:**

.. code-block:: python

   >>> conkit.io.write('toxd/toxd.fasta', 'fasta', msa)

**3. To convert file formats in single call, you can use the :func:`conkit.io.convert` function.**

.. code-block:: python

   >>> conkit.io.convert('toxd/toxd.a3m', 'a3m', 'toxd/toxd.fasta', 'fasta')

You can convert these files to many different other formats, for a full list check out the :ref:`file_formats`.
