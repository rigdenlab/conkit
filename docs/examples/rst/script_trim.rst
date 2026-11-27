.. _script_trim:

Structure Trimming
------------------

``conkit-trim`` trims a set of structure files to keep only residues involved in highly conserved contacts. It takes a contact probability map (typically produced by ``conkit-summarise``) and retains only residues that appear in contacts whose conservation score meets a threshold. This is useful for reducing structures to their most informative, well-conserved core before further analysis.

Basic usage
^^^^^^^^^^^

.. code-block:: bash

   $> conkit-trim summary.mat mapalign "models/*.pdb" pdb \
          --output_prefix trimmed

The first two positional arguments are the contact map and its format; the third and fourth are the glob pattern and format for the structure files to trim. Each trimmed structure is written as ``trimmed_<wildcard>.<format>``, where ``<wildcard>`` is the part of the original filename that matched the glob pattern.

Adjusting the conservation threshold
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

By default only contacts with a score ≥ 0.9 (present in at least 90 % of the ensemble) are considered. Use ``--score_threshold`` to change this:

.. code-block:: bash

   $> conkit-trim summary.mat mapalign "models/*.pdb" pdb \
          --score_threshold 0.7 \
          --output_prefix trimmed_70pct

Lowering the threshold keeps more residues; raising it produces a more compact, highly-conserved core.

Using a structure as the contact reference
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The contact map can also be derived from a single structure rather than a probability map:

.. code-block:: bash

   $> conkit-trim reference.pdb pdb "models/*.pdb" pdb \
          --output_prefix trimmed

RNA and custom contact definitions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For RNA, pass ``--moltype RNA`` to use C1′–C1′ distances with a 12.5 Å cutoff when reading a structure-format contact map:

.. code-block:: bash

   $> conkit-trim rna_summary.mat mapalign "rna_models/*.cif" mmcif \
          --moltype RNA \
          --output_prefix rna_trimmed

.. note::

   When the contact map is provided in a non-structure format (e.g. ``mapalign``), the contact definition stored in the file is used and ``--moltype``, ``--rep_atom``, and ``--contact_dist`` have no effect on the contact map itself. They only affect how structure-format contact maps are read.

Optional flags
^^^^^^^^^^^^^^^

``--score_threshold FLOAT``
  Minimum contact score (conservation fraction) to keep. Residues in contacts below this threshold are removed. Default: ``0.9``.

``--output_prefix PREFIX``
  Prefix for output filenames. Each trimmed structure is saved as ``<prefix>_<wildcard>.<format>``. Default: ``trimmed``.

``--moltype Protein|RNA|Custom``
  Molecule type controlling contact definition when reading a structure-format contact map. Default: ``Protein``.

``--contact_dist FLOAT``
  Distance cutoff in Ångstrom when using ``--moltype Custom``.

``--rep_atom ATOM``
  Representative atom name when using ``--moltype Custom`` (e.g. ``CA``, ``C1'``).

``--overwrite``
  Overwrite output files if they already exist.
