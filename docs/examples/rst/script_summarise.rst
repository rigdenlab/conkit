.. _script_summarise:

Contact Map Summarisation
--------------------------

``conkit-summarise`` builds a consensus contact probability map from a set of structure files. Each contact's score in the output map is the fraction of input structures in which that contact is present. This is useful for identifying contacts that are conserved across an ensemble of models (e.g. from a molecular dynamics trajectory, a set of homology models, or multiple AlphaFold predictions).

Basic usage
^^^^^^^^^^^

.. code-block:: bash

   $> conkit-summarise "models/*.pdb" pdb \
          --output summary.png \
          --output_contact_map summary.mat

The first positional argument is a glob pattern (quoted to prevent shell expansion) that matches all structure files to summarise; the second is their format. The tool writes a contact map probability figure (``summary.png``) and the summary contact matrix in map_align format (``summary.mat``).

Each contact's score in ``summary.mat`` is the fraction of structures (0.0–1.0) in which that contact was present. Contacts absent from all structures are not written.

Overlaying a reference structure
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A single reference structure can be overlaid on the summary map for comparison:

.. code-block:: bash

   $> conkit-summarise "models/*.pdb" pdb \
          --overlay_structure reference.pdb \
          --overlay_struct_type pdb \
          --output summary_overlay.png

The reference contacts appear as a second layer on the same matrix plot, allowing visual comparison between the consensus and a specific structure.

Comparing two ensembles (split diagonal)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To place a second ensemble below the diagonal of the matrix plot, use ``--other_files``:

.. code-block:: bash

   $> conkit-summarise "ensemble_A/*.pdb" pdb \
          --other_files "ensemble_B/*.pdb" \
          --other_struct_type pdb \
          --output comparison.png

The upper triangle shows the first ensemble; the lower triangle shows the second.

RNA and custom contact definitions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

By default the tool uses Cα–Cα distances with a cutoff of 8 Å (protein convention). For RNA, pass ``--moltype RNA`` to switch to C1′–C1′ distances with a 12.5 Å cutoff:

.. code-block:: bash

   $> conkit-summarise "rna_models/*.cif" mmcif \
          --moltype RNA \
          --output rna_summary.png

Custom atom and distance definitions can be set with ``--rep_atom`` and ``--contact_dist`` when using ``--moltype Custom``.

Optional flags
^^^^^^^^^^^^^^^

``--output PATH``
  Output figure file. Default: ``conkit_summary.png``.

``--output_contact_map PATH``
  Output contact matrix file. Default: ``conkit_summary.mat``.

``--output_contact_map_format FORMAT``
  Format for the output contact matrix. Default: ``mapalign``.

``--color_map CMAP``
  Matplotlib colormap for the figure. Default: ``Greys``.

``--overlay_structure PATH`` / ``--overlay_struct_type FORMAT``
  A single additional structure to overlay on the figure.

``--other_files PATTERN`` / ``--other_struct_type FORMAT``
  A second set of structures to display below the diagonal.

``--moltype Protein|RNA|Custom``
  Molecule type controlling contact definition. Default: ``Protein``.

``--contact_dist FLOAT``
  Distance cutoff in Ångstrom when using ``--moltype Custom``.

``--rep_atom ATOM``
  Representative atom name when using ``--moltype Custom`` (e.g. ``CA``, ``C1'``).

``--overwrite``
  Overwrite the output figure if it already exists.
