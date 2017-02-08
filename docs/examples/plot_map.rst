.. _example_plotting_a_map:

Plotting a contact map
----------------------

Using a script
^^^^^^^^^^^^^^

If you would like to plot a contact map using ConKit without the overhead of using Python, you can simply use the ``conkit.plot_map`` script.

.. code-block:: bash

   $> conkit.plot_map toxd/toxd.fasta fasta toxd/toxd.mat ccmpred

The call above uses the contact prediction file ``toxd.mat`` file, which is in ``ccmpred`` format, and plots the following 2D contact map stored in the file ``toxd/toxd.png``

.. _Toxd Simple Contact Map:

.. image:: images/toxd_cmap_simple.png
   :alt: Toxd CMap Simple

--------------------------------------------------------------

You can also add a reference structure to determine which contacts are true and false positive contacts. By default, all contacts are identified in the reference structure by measuring the distance between Cβ atoms, whereby all atoms closer than 8Å are considered to be in contact.

.. code-block:: bash

   $> conkit.plot_map -p toxd/toxd.pdb toxd/toxd.fasta fasta toxd/toxd.mat ccmpred

The call above produces a contact map plot looking like this. The gray points are the reference contacts, green show true positive contacts in your prediction and red false positive ones.

.. _Toxd Reference Contact Map:

.. image:: images/toxd_cmap_reference.png
   :alt: Toxd CMap Reference

--------------------------------------------------------------

You could also add a second contact prediction file to the call to compare two maps against each other.

.. code-block:: bash

   $> conkit.plot_map -e toxd/toxd.psicov -ef psicov -p toxd/toxd.pdb toxd/toxd.fasta fasta toxd/toxd.mat ccmpred

The call above produces a contact map plot looking like this. The gray points are the reference contacts, green show true positive contacts in your prediction and red false positive ones. The top triangle is the second contact map from file ``toxd/toxd.psicov`` whereas the bottom one is from ``toxd/toxd.mat``.

.. image:: images/toxd_cmap_advanced.png
   :alt: Toxd CMap Advanced

--------------------------------------------------------------

Finally, you could also illustrate the confidence with which each contact was predicted.

.. code-block:: bash

   $> conkit.plot_map --confidence -e toxd/toxd.psicov -ef psicov -p toxd/toxd.pdb toxd/toxd.fasta fasta toxd/toxd.mat ccmpred

The call above produces a contact map plot looking like this. All parameters and settings are identical to the previous map except the ``--confidence`` flag, which will show more confidently predicted contacts as larger markers.

.. image:: images/toxd_cmap_confidence.png
   :alt: Toxd CMap Confidence

.. note::

   You can use the last two examples also **without** a reference structure!

--------------------------------------------------------

Using Python
^^^^^^^^^^^^

Two simplified versions of the ``conkit.plot_map`` script is shown below to illustrate how you can plot your own contact map using Python.

.. only:: html

   .. note::

      You can download this script :download:`here <code/map_plotting_simple.py>`.

.. literalinclude:: /../docs/examples/code/map_plotting_simple.py
   :language: python
   :linenos:

This will produce the `Toxd Simple Contact Map`_ plot.

--------------------------------------------------------------

.. only:: html

   .. note::

      You can download this script :download:`here <code/map_plotting_reference.py>`.

.. literalinclude:: /../docs/examples/code/map_plotting_reference.py
   :language: python
   :linenos:

This will produce the `Toxd Reference Contact Map`_ plot.
