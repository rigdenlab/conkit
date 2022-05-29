.. _script_convert_distpred:

File Format Conversion
----------------------

If you would like to convert a file from one format to another, you can do so using ConKit's provided script, which is called ``conkit-convert``.

.. code-block:: bash

   $> conkit-convert 7l6q/7l6q.af2 alphafold2 7l6q/7l6q.rr2 caspmode2

The call above converts the ``7l6q.af2`` file, which is in ``alphafold2`` format, to the ``7l6q.rr2`` file in ``caspmode2`` format. These two file formats contain inter-residue distance predictions. You can also choose to convert residue distance predictions into a residue contact prediction file, for example:

.. code-block:: bash

   $> conkit-convert 7l6q/7l6q.af2 alphafold2 7l6q/7l6q.psicov psicov

The call above converts the ``7l6q.af2`` file, which is in ``alphafold2`` distance prediction format, to the ``7l6q.psicov`` file in ``psicov`` contact prediction format. It is not possible to perform the inverse conversion (residue contact predictions into inter-residue distance predictions).

Conkit supports many different other formats, for a full list check out the :ref:`file_formats`.
