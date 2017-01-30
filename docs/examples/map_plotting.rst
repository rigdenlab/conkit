.. _example_plotting_a_map:

Contact Map Plotting
====================

To plot contacts in form of a contact map using ConKit, you need a contact prediction file.

.. warning::
   You require the `Matplotlib <http://matplotlib.org/>`_ package to use this functionality. If you are unsure if it is installed on your system, refer to the :ref:`installation` documentation


.. note::
   ConKit provides a script that allows you to plot contact maps from the command line. The script is called ``conkit.plot_map`` and installed automatically.


To start, you need to read the file into its hierarchy.

.. code-block:: python

   >>> import conkit
   >>> contact_h = conkit.io.read('pathto/contact.file', 'contact_format')

With the code above, we created a :obj:`conkit.core.ContactFile` hierarchy. This allows you to store contact maps in a single hierarchy. Here, we are only interested in the first, and thus we can remove the rest.

.. code-block:: python

   >>> contact_map = contact_h.top_map

Finally, to plot the contact map, all we need to do is invoke the relevant function.

.. code-block:: python

   >>> contact_map.plot_map()

This will create a contact map plot that will be saved to your disk.
