.. _example_file_conversion:

Converting file formats
-----------------------

Using a script
^^^^^^^^^^^^^^

If you would like to convert a file from one format to another, you can do so using ConKit's provided script, which is called ``conkit.convert``.

.. code-block:: bash

   $> conkit.convert toxd/toxd.mat ccmpred toxd/toxd.rr casprr

The call above converts the ``toxd.mat`` file, which is in ``ccmpred`` format, to the ``toxd.rr`` file in ``casprr`` format.

You can convert these files to many different other formats, for a full list check out the :ref:`file_formats`.

--------------------------------------------------------

Using Python
^^^^^^^^^^^^

In order to convert files in ConKit, we need to use the ConKit I/O framework.

.. note::
   ConKit I/O framework consists of three main functions that handle the relevant parsers: :func:`conkit.io.read`, :func:`conkit.io.write` and :func:`conkit.io.convert`. The latter effectively uses the former two but handles everything in one step.

1. Files can be read in ConKit's internal hierarchies using simple Python code.

.. code-block:: python
  
   >>> conpred = conkit.io.read('toxd/toxd.mat', 'ccmpred')

2. Contact prediction hierarchies can also be written in a similarly easy format. Using the ``conpred`` hierarchy we have created above:

.. code-block:: python

   >>> conkit.io.write('toxd/toxd.rr', 'casprr', conpred)

3. To convert file formats in single call, you can use the :func:`conkit.io.convert` function.

.. code-block:: python

   >>> conkit.io.convert('toxd/toxd.mat', 'ccmpred', 'toxd/toxd.rr', 'casprr')

You can convert these files to many different other formats, for a full list check out the :ref:`file_formats`.