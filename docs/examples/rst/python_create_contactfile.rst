.. _python_create_contactfile:

ConKit ContactFile Hierarchy Construction
-----------------------------------------

If you wish to construct it as part of a new development to store your contact information, you might find the following helpful.

Entities
++++++++

**1. How to create a :obj:`~conkit.core.contact.Contact`?**

.. code-block:: python

   >>> from conkit.core import Contact
   >>> contact = Contact(1, 10, 1.0)

The example above creates a contact between residues ``1`` and ``10`` and assigns a :attr:`~conkit.core.contact.Contact.raw_score` of ``1.0`` to it. By default, this contact has many more default attributes assigned, such as the distance value often seen in columns 3 and 4 in the Casp RR format.

**2. How to create a :obj:`~conkit.core.contactmap.ContactMap`?**

.. code-block:: python

   >>> from conkit.core import ContactMap
   >>> cmap = ContactMap('example')

This example shows you how to create a :obj:`~conkit.core.contactmap.ContactMap` which can store a :obj:`~conkit.core.contact.Contact`.

**3. How to create a :obj:`~conkit.core.contactfile.ContactFile`?**

.. code-block:: python

   >>> from conkit.core import ContactFile
   >>> cmap = ContactFile('example')

This example shows you how to create a :obj:`~conkit.core.contactfile.ContactFile` which can store a :obj:`~conkit.core.contactmap.ContactMap`.

Hierarchy
+++++++++

Above is an outline for the different contact-related :obj:`~conkit.core.entity.Entity` classes. Each higher entity allows you to store one or more lower-level ones, i.e. you can store one or more :obj:`~conkit.core.contactmap.ContactMap` entities in a single :obj:`~conkit.core.contactfile.ContactFile`. Similarly, you could many :obj:`~conkit.core.contact.Contact` entities in a :obj:`~conkit.core.contactmap.ContactMap`; however, be aware that all **must** have unique IDs.

To illustrate how you can combine the :obj:`~conkit.core.entity.Entity` classes, look at the following:

.. code-block:: python

   >>> from conkit.core import Contact, ContactMap, ContactFile
   >>> cfile = ContactFile('example_file')
   >>> cmap = ContactMap('example_map')
   >>> contact = Contact(1, 10, 1.0)
   >>> # Add the contact to the contact map
   >>> cmap.add(contact)
   >>> # Add the contact map to the contact file
   >>> cfile.add(cmap)

Note, the order in which you add :obj:`~conkit.core.entity.Entity` instances does not matter. We could also add the ``cmap`` to the ``cfile`` before adding the ``contact`` to the ``cmap``.

Once you have constructed your hierarchy, all related functions are available to you.
