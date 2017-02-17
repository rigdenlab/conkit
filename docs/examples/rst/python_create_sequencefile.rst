
ConKit SequenceFile Hierarchy Construction
------------------------------------------

If you wish to construct it as part of a new development to store your sequence information, you might find the following helpful.

Entities
++++++++

1. **How to create a Sequence?**

.. code-block:: python

   >>> import conkit
   >>> seq = conkit.core.Sequence("example", "ABCDEF")

The example above creates a :obj:`Sequence <conkit.core.Sequence>` object with id "example" and sequence "ABCDEF". This sequence contains numerous functions, such as :func:`align_global() <conkit.core.Sequence.align_global>` for a global pairwise alignment with a second sequence

2. **How to create a SequenceFile?**

.. code-block:: python

   >>> import conkit
   >>> sfile = conkit.core.SequenceFile("example)

This example shows you how to create a :obj:`SequenceFile <conkit.core.SequenceFile>` which can store one or more :obj:`Sequence <conkit.core.Sequence>` objects.

Hierarchy
+++++++++

Above is an outline for the different sequence-related entities. Each higher entity allows you to store one or more lower-level ones, i.e. you can store one or more :obj:`Sequence <conkit.core.Sequence>` entities in a single :obj:`SequenceFile <conkit.core.SequenceFile>`; however, be aware that all **must** have unique IDs.

To illustrate how you can combine the entities, look at the following:

.. code-block:: python

   >>> import conkit
   >>> sfile = conkit.core.SequenceFile("example)
   >>> seq1 = conkit.core.Sequence("example", "ABCDEF")
   >>> seq2 = conkit.core.Sequence("elpmaxe", "FEDCBA")
   >>> # Add the sequences to the sequence file
   >>> sfile.add(seq1)
   >>> sfile.add(seq2)

Note, the order in which you add entities does not matter.

Once you have constructed your hierarchy, all related functions are available to you.
