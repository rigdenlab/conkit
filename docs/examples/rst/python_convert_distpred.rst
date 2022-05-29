.. _python_convert_distpred:

File Format Conversion
----------------------

In order to convert files in ConKit, we need to use the ConKit I/O framework.

.. note::
   ConKit I/O framework consists of three main functions that handle the relevant parsers: :func:`~conkit.io.read`, :func:`~conkit.io.write` and :func:`~conkit.io.convert`. The latter effectively uses the former two but handles everything in one step.

**1. Files can be read in ConKit's internal hierarchies using simple Python code.**

.. code-block:: python

   >>> import conkit.io
   >>> distpred = conkit.io.read('7l6q/7l6q.af2', 'alphafold2')

**2. Residue distance prediction hierarchies can also be written in a similarly easy format. Using the ``distpred`` hierarchy we have created above:**

.. code-block:: python

   >>> import conkit.io
   >>> conkit.io.write('7l6q/7l6q.casprr2', 'caspmode2', distpred)

**3. To convert file formats in single call, you can use the :func:`~conkit.io.convert` function.**

.. code-block:: python

   >>> import conkit.io
   >>> conkit.io.convert('7l6q/7l6q.af2', 'alphafold2', '7l6q/7l6q.casprr2', 'caspmode2')

You can convert these files to many different other formats, for a full list check out the :ref:`file_formats`. It is possible to convert distance prediction files into contact prediction files, but not the other way around.
