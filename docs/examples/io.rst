.. _example_file_conversion:

File Conversions
================

In order to convert files in ConKit, we need to import the ConKit I/O framework.

.. code-block:: python

   >>> from conkit import io

.. note::
   ConKit I/O framework consists of three main functions that handle the relevant parsers: :func:`conkit.io.read`, :func:`conkit.io.write` and :func:`conkit.io.convert`. The latter effectively uses the former but handles everything in one step.

--------------------------------------------------------

Contact Prediction Files
------------------------

File reading
^^^^^^^^^^^^

Let's say you are designing a new pipeline and you want to manipulate your contact file. Thus, you need to have a way of reading the file into memory before you can manipulate it.

ConKit provides the following simple syntax to do exactly that.

.. code-block:: python
  
   >>> with open('query.pconsc3.txt', 'r') as f_in:
   ...     contact_hierarchy = io.read(f_in, 'pconsc3')

File writing
^^^^^^^^^^^^

Now, you've manipulated the file from the previous example and you now want to write it back out. Look no further, this is what you need to do.

.. code-block:: python
   
   >>> with open('query.pconsc3.new.txt', 'w') as f_out:
   ...     io.write(f_out, contact_hierarchy, 'pconsc3')

File conversion
^^^^^^^^^^^^^^^

Imagine the following scenario: you have just used an online server, let's say `PconsC3`_, and you have received a contact prediction file. Now, the file in text format is little use to you and although you have a contact map, there are many cool tools that use predictions to visualise them in a protein structure. These tools rely on the official Casp format and you really don't want to convert every line manually or write a script to do it. Look no further, ConKit does this in a heartbeat.

Assuming our PconsC3 contact file is called ``query.pconsc3.txt`` and you want to write it to a file called ``query.rr.txt``, you can convert it using the following code:

.. code-block:: python
   
   >>> with open('query.pconsc3.txt', 'r') as f_in, open('query.casp.rr', 'w') as f_out:
   ...     io.convert(f_in, 'pconsc3', f_out, 'casprr')

This is it, you now have converted your contact file from PconsC3 format to Casp format.

--------------------------------------------------------

(Multiple) Sequence Files
-------------------------

The functions used to convert sequence files are identical to those used for contact prediction files. The only difference is the format parsed to those functions.


.. _PconsC3: http://pconsc3.bioinfo.se/

 
