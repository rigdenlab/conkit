.. _example_constructing_a_hierarchy:

Hierarchy Construction
======================

.. note::
   The easiest way to construct a contact file hierarchy is to choose the corresponding parser.

.. warning::
   This example is kept brief, if you are unable to follow the process or want to know more, check out the source code.

If you wish to construct it as part of a new development to store your contact information, you might find the following helpful.

1. Constructing instances
^^^^^^^^^^^^^^^^^^^^^^^^^

The following blocks of code a purely an outline on how each instance can be created.

1.1 ContactFile
...............

.. code-block:: python

   >>> from conkit.core import ContactFile
   >>> contact_file = ContactFile('<YOUR_ID_HERE>')

1.2 ContactMap
..............

.. code-block:: python

   >>> from conkit.core import ContactMap
   >>> contact_map = ContactMap('<YOUR_ID_HERE')

1.3 Contact
...........

.. code-block:: python

   >>> from conkit.core import Contact
   >>> contact = Contact('res1_seq', 'res2_seq', 'raw_score')

.. note::
   The :obj:`conkit.core.Contact` ID is automatically assigned based on ``res1_seq`` and ``res2_seq`` attributes.

2. Assembling the hierarchy
^^^^^^^^^^^^^^^^^^^^^^^^^^^

To assemble a hierarchy, you want to add :obj:`conkit.core.Contact` instances to a :obj:`conkit.core.ContactMap` instance. Then, you can add one or more :obj:`conkit.core.ContactMap` instances to a single :obj:`conkit.core.ContactFile`.

Be aware, the IDs for all instances need to be unique at their level, but can be repeated across multiple instances in higher levels. I.e. a ``Contact(1, 2, 0.1)`` can be added only once to ``ContactMap('map_1')`` but added to ``ContactMap('map_2')`` if desired.

.. code-block:: python

   >>> from conkit.core import Contact, ContactMap, ContactFile
   >>> contact_file = ContactFile('example')
   >>> for i in range(1, 4):
   ...     contact_map = ContactMap('map_{0}'.format(i))
   ...
   ...     for j in range(1, 11):
   ...         contact = Contact(i, j, 0.1)
   ...         contact_map.add(contact)
   ...
   ...     contact_file.add(contact_map)

In the example above we create a single contact file hierarchy that will contain three contact maps with each 10 contacts. Can you see it?

