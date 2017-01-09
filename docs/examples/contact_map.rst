.. _example_plotting_a_map:

Contact Map plotting
====================

To plot contacts in form of a contact map using ConKit, you need a contact prediction file.

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
