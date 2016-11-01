..

**************************
Contact Prediction ToolKit
**************************

**A Python Interface for contact predictions**

.. warning::

   ConKit is still undergoing major development and changes prior to a 1.0 release will happen. It has also not yet been extensively validated.

This project is a result of the continuous struggle using residue-residue contact prediction pipelines, visualisation tools and related software. The aim was to reduce this complexity and provide one unified interface as a basic platform. The resulting platform is the Contact Prediction ToolKit, or ConKit in short.

Key Features
~~~~~~~~~~~~

- Parsers for Multiple Sequence Alignment and contact prediction files
- Easy-to-use wrappers for the most-commonly used software, e.g. HHblits, Jackhmmer, CCMpred, PSICOV, ...
- Many convenience functions for prediction-structure matching, determination of True Positive contacts, calculation of MSA and contact prediction statistics, etc.

Contributing
~~~~~~~~~~~~
There are two ways by which you can contribute to ConKit:

1. Submit any suggestions to the `GitHub Issue Tracker`_, or
2. Fork this repository, commit your changes and submit a pull request.

Found a Bug?
~~~~~~~~~~~~
Please use the `GitHub Issue Tracker`_.

--------------------------------------------------------------------------------

Documentation
~~~~~~~~~~~~~

.. toctree::
   :caption: Table of Contents
   :maxdepth: 1

   installation
   examples
   scripts
   pythonapi

--------------------------------------------------------------------------------

Links
~~~~~

- `ConKit GitHub Repository`_

--------------------------------------------------------------------------------

Ackowledgements
~~~~~~~~~~~~~~~
- `Jens Thomas`_
- `Stefan Seemayer`_
- `BioPython`_



.. _BioPython: https://www.biopython.org
.. _ConKit GitHub Repository: https://github.com/fsimkovic/conkit
.. _GitHub Issue Tracker: https://github.com/fsimkovic/conkit/issues
.. _Jens Thomas: http://www.farmurban.co.uk/about
.. _Stefan Seemayer: https://github.com/sseemayer
