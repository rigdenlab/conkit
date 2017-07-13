
Changelog
=========

v0.8
----

- Default value in ``calculate_meff()`` and ``calculate_weights()`` changed from 0.7 to 0.8 [more commonly used in literature]
- Bug fix with PyPi installation where ``requirements.txt`` not found; fix includes removal of ``requirements.txt`` and addition of ``install_requires`` to ``setup.py`` instead. - Thanks to Miguel Correa for reporting this bug
- Logging message coloring according to message level
- ``filter()`` function added for redundancy/distant homolog removal from ``SequenceFile``
- ``io`` sub-package caches modules and imports upon request
- ``core`` classes extracted to individual module files
- License text added to each module


----
