
Changelog
=========

v0.8.3dev1
----------
- ``requirements.txt`` file re-added for easier dependency installation
- Distance definitions accept floating point values
- ``LinearBW`` calculator added for linear bandwidth calculation in analysis
- ``_BandwidthCalc`` class renamed to ``BandwidthBase``
- Abstractified ``BandwidthBase``, and ``Parser`` with all subparser classes 
- Bug fix
    - ``PconsParser`` class accepts negative ``raw_score`` values
    - ``SequenceFile.neff`` returns ``float`` instead of ``int``
    - ``CCMpredParser.read()`` returns empty ``ContactFile`` when matrix file empty

v0.8.2
------

- Critical bug fix for automated opening of filehandle in Python2.7 
- Test function skipping added for ``SequenceFile.filter()`` when SciPy not installed
- Renamed conkit/io/tests files for filenames to agree with modules in conkit/io
- Performance of ``write()`` in parsers improved by construction of string and single call to ``write()`` of filehandle

v0.8.1
------

- Revoked catching of ``SystemExit(0)`` exception in scripts when invoked with ``--help`` flag 
- Bug fix relating to Python3 automatic opening of file handles - Thanks to Miguel Correa for reporting this bug

v0.8
----

- Default value in ``calculate_meff()`` and ``calculate_weights()`` changed from 0.7 to 0.8 [more commonly used in literature]
- Bug fix with PyPi installation where ``requirements.txt`` not found; fix includes removal of ``requirements.txt`` and addition of ``install_requires`` to ``setup.py`` instead. - Thanks to Miguel Correa for reporting this bug
- Logging message coloring according to message level
- ``filter()`` function added for redundancy/distant homolog removal from ``SequenceFile``
- ``io`` sub-package caches modules and imports upon request
- ``core`` classes extracted to individual module files
- License text added to each module
