
ConKit ContactFile Hierarchy Construction
-----------------------------------------

If you wish to construct it as part of a new development to store your contact information, you might find the following helpful.

Entities
++++++++

**1. How to create a Contact?**

.. code-block:: python

   >>> import conkit
   >>> contact = conkit.core.Contact(1, 10, 1.0)

The example above creates a contact between residues ``1`` and ``10`` and assigns a ``raw_score`` of ``1.0`` to it. By default, this contact has many more default attributes assigned, such as the distance value often seen in columns 3 and 4 in the Casp RR format.

**2. How to create a ContactMap?**

.. code-block:: python

   >>> import conkit
   >>> cmap = conkit.core.ContactMap('example')

This example shows you how to create a :obj:`ContactMap <conkit.core.ContactMap>` which can store a :obj:`Contact <conkit.core.Contact>`.

**3. How to create a ContactFile?**

.. code-block:: python

   >>> import conkit
   >>> cmap = conkit.core.ContactFile('example')

This example shows you how to create a :obj:`ContactFile <conkit.core.ContactFile>` which can store a :obj:`ContactMap <conkit.core.ContactMap>`.

Hierarchy
+++++++++

Above is an outline for the different contact-related entities. Each higher entity allows you to store one or more lower-level ones, i.e. you can store one or more :obj:`ContactMap <conkit.core.ContactMap>` entities in a single :obj:`ContactFile <conkit.core.ContactFile>`. Similarly, you could many :obj:`Contact <conkit.core.Contact>` entities in a :obj:`ContactMap <conkit.core.ContactMap>`; however, be aware that all **must** have unique IDs.

To illustrate how you can combine the entities, look at the following:

.. code-block:: python

   >>> import conkit
   >>> cfile = conkit.core.ContactFile('example_file')
   >>> cmap = conkit.core.ContactMap('example_map')
   >>> contact = conkit.core.Contact(1, 10, 1.0)
   >>> # Add the contact to the contact map
   >>> cmap.add(contact)
   >>> # Add the contact map to the contact file
   >>> cfile.add(cmap)

Note, the order in which you add entities does not matter. We could also add the ``cmap`` to the ``cfile`` before adding the ``contact`` to the ``cmap``.

Once you have constructed your hierarchy, all related functions are available to you.
