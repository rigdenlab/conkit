
Changelog
=========

[Unreleased]
------------
Added
~~~~~
- ``Entity.top`` property to always return the first child in the list
- ``ContactMap.find`` function accepts ``strict`` keyword argument to find contact pairs with both residues in ``register``
- ``PdbParser`` takes a distance cutoff of ``0`` to include all Cb-Cb contacts in the protein structure
- ``ContactMatchState`` enumerated type for definitions of state constants for contact
- ``SequenceAlignmentState`` enumerated type for definitions of state constants for each sequence file 
Changed
~~~~~~~
- Optimized some functions and comparisons according to the recommended Python optimization instructions 
- ``ContactMap.match`` does __not__ modifiy ``other`` by default anymore. Specify ``match_other=True`` as kwarg!
- ``ContactMap.calculate_kernel_density`` renamed to ``ContactMap.calculate_contact_density`` 

[0.8.3]
-------
Added
~~~~~
- ``requirements.txt`` file re-added for easier dependency installation
- ``LinearBW`` calculator added for linear bandwidth calculation in analysis
- ``seq_ascii`` property to ``Sequence`` for encoded sequence
- ``ascii_matrix`` property to ``SequenceFile`` for encoded alignment 
- ``SequenceFile`` and ``ContactFile`` classes have new ``empty`` properties
- ``flib`` format for ``ContactFile`` classes to allow easier conversions for the Flib-Coevo fragment picking library
Changed
~~~~~~~
- Distance definitions accept floating point values
- ``_BandwidthCalc`` class renamed to ``BandwidthBase``
- Abstractified ``BandwidthBase``, and ``Parser`` with all subparser classes 
- Refactored ``conkit/io/__init__.py`` to avoid duplication of code
Fixed
~~~~~
- ``PconsParser`` class accepts negative ``raw_score`` values
- ``SequenceFile.neff`` returns ``float`` instead of ``int``
- ``CCMpredParser.read()`` returns empty ``ContactFile`` when matrix file empty

[0.8.2]
-------
Added
~~~~~
- Test function skipping added for ``SequenceFile.filter()`` when SciPy not installed
Changed
~~~~~~~
- Renamed conkit/io/tests files for filenames to agree with modules in conkit/io
- Performance of ``write()`` in parsers improved by construction of string and single call to ``write()`` of filehandle
Fixed
~~~~~
- Critical bug fix for automated opening of filehandle in Python2.7 

[0.8.1]
-------
Changed
~~~~~~~
- Revoked catching of ``SystemExit(0)`` exception in scripts when invoked with ``--help`` flag 
Fixed
~~~~~
- Bug fix relating to Python3 automatic opening of file handles - Thanks to Miguel Correa for reporting this bug

[0.8]
-----
Added
~~~~~
- Logging message coloring according to message level
- ``filter()`` function added for redundancy/distant homolog removal from ``SequenceFile``
- License text added to each module
- ``io`` sub-package caches modules and imports upon request
Changed
~~~~~~~
- Default value in ``calculate_meff()`` and ``calculate_weights()`` changed from 0.7 to 0.8 [more commonly used in literature]
- ``core`` classes extracted to individual module files
Fixed
~~~~~
- Bug fix with PyPi installation where ``requirements.txt`` not found; fix includes removal of ``requirements.txt`` and addition of ``install_requires`` to ``setup.py`` instead. - Thanks to Miguel Correa for reporting this bug
