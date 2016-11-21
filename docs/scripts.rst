.. _scripts:


Scripts
=======

ConKit is shipped with a small number of convenience scripts that are intended to give examples as well as provide convenience wrappers for common tasks.

The following scripts are currently shipped with ConKit:

- :obj:`conkit.convert` - script to convert file formats
- :obj:`conkit.msatool` - script to analyse a Multiple Sequence Alignment 
- :obj:`conkit.predict` - simple contact prediction pipeline using HHblits, HHfilter and CCMpred

All scripts can be invoked from the command line using the name listed above, i.e to convert a prediction file to the Casp RR format you can call

.. code-block:: bash

   $> conkit.convert input.mat ccmpred output.rr casprr
