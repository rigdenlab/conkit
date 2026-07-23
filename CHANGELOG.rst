**[0.14.1]**

*Added*

- Residue numbering anomaly detection in ``conkit-validate``: MISUSED_ICODE and EXTRA_CANONICAL events are flagged in the terminal table, the JSON output, and as annotated vertical lines in the validation figure
- Back-mapping of original PDB residue numbers in the terminal table (``orig→new`` format) and in JSON output (``orig_resnum`` field per residue)

*Changed*

- ``ModelValidationFigure`` x-axis label updated to "Residue Number (FASTA position)" to clarify that residue positions refer to the provided FASTA sequence

*Fixed*

- RNA alignment in ``conkit-validate``: BLOSUM62 substitution matrix (protein-only alphabet, no U) replaced with simple match/mismatch scoring for non-protein moltype, fixing ``ValueError`` for any RNA sequence containing U
- Renumbering fallthrough: ``write_renumbered_version_of_chain_in_struct`` now raises ``ValueError`` when no chain aligns, rather than silently returning an empty result
- Scaler files resaved with plain ``pickle`` to avoid NumPy 2.5 deprecation warnings from joblib's ``NumpyArrayWrapper`` when loading

**[0.14.0]**

*Added*

- RNA model validation in ``conkit-validate``: ``--moltype RNA`` selects RNA-specific Random Forest classifiers (note: referred to as SVM in code for historical reasons; a Random Forest was found to outperform an SVM for RNA); supports AlphaFold 3 C1′ distograms as prediction input
- DNATCO integration for RNA secondary-structure annotation: ``--dnatco_exe`` provides the DNATCO node executable, enabling RNA-specific classifier features
- ``include_hetatms`` flag in PDB/mmCIF reader to include modified nucleotides (e.g. PSU, M2G) and other HETATM records in distogram construction
- RNA-aware residue handling: modified nucleotide codes of three or more characters are now accepted by the ``Contact`` residue setter
- CMO (map_align) and SVM/RF scoring can now be run independently in ``conkit-validate`` via separate command line flags
- False-positive indicators in ``conkit-validate``: contact counting per residue, pLDDT reading from predicted structures, and GESAMT Q-score calculation — these help identify regions flagged as errors due to missing experimental density or domain orientation differences rather than genuine modelling errors
- ``--min_err_size`` and ``--svm_threshold`` flags in ``conkit-validate`` to control the minimum run length and score threshold for called errors
- Generalised contact-definition flags in ``conkit-validate`` and the underlying ``set_contact_definition`` utility
- New command line tool ``conkit-summarize``: takes multiple structures or contact maps for the same molecule and computes an averaged contact probability map; intended for identifying contacts conserved across a class of conformations (e.g. amyloid polymorphs)
- New command line tool ``conkit-trim``: trims a structure to retain only residues that participate in contacts scoring above a given threshold on a provided contact probability map

*Changed*

- Packaging modernised: project metadata, entry points, and dependencies consolidated into ``pyproject.toml``; ``setup.py`` reduced to Cython extensions only; ``setup.cfg`` removed
- Supported Python versions updated to 3.9–3.12
- Replaced deprecated ``Bio.pairwise2`` with ``Bio.Align.PairwiseAligner`` throughout
- wRMSD calculation from structures now applies the same distance cutoff that arises naturally when computing wRMSDs from distograms, making the two comparable and improving protein classifier performance when a structure is provided instead of a distogram
- Annotation bars in ``ModelValidationFigure`` now stack dynamically: rows absent because a feature is disabled no longer leave visual gaps; the y-axis lower limit adjusts to the bars actually drawn
- DNATCO subprocess output now written to an isolated temporary directory, avoiding working-directory side effects

*Fixed*

- Removed deprecated numpy type aliases ``np.int`` and ``np.bool`` (replaced with ``int``/``bool``)
- Fixed ``scikit-learn`` version incompatibility in bundled classifiers and scalers

**[0.13.3]**

*Added*

-  ``conkit.core.ContactMap.match_naive`` method for contact map match when no sequence alignment is required
- Examples on how to use ``conkit-validate`` in documentation at conkit.org
- Examples on how file conversions for ditances in documentation at conkit.org
- Examples on how to plot residue ditances in documentation at conkit.org

*Changed*

- Update ``requirements.txt`` to include versions of biopython and sklearn compatible with CCP4 8.0
- Update requirements list in documentation at conkit.org

*Fixed*

- Resolve contact map match when one of the input maps is empty

**[0.13.2]**

*Fixed*

- Further fixes to pip install package

**[0.13.1]**

*Fixed*

- Minor fix for pip install and cython extension

**[0.13]**

*Added*

- Added support for distance prediction files
- Added new visualisation plots for distance files
- Added new command line tool for model validation conkit-validate

*Changed*

- Remove support for Python3.6
- Add support for Python3.9

**[0.12]**

*Fixed*

- Resolve plotting of small contact maps

*Changed*

- Remove support for Python2.7
- Remove support for Python3.5
- Add support for Python3.8

**[0.11.3]**

*Fixed*

- Test cases ensure file removal regardless of failure

*Changed*

- Code formatting to adapt [Black](https://black.readthedocs.io/en/stable/) formatting

*Added*

- [``map_align``](https://github.com/sokrypton/map_align) contact file parser
- AppVeyor and TravisCI runs against Python3.7

**[0.11.2]**

*Fixed*

- Bug fix to avoid rare ``ZeroDivision``

**[0.11.1]**

*Changed*

- ``conkit/core/ext/c_sequencefile.pyx`` removed ``print`` statement

**[0.11]**

*Added*

- ``conkit.io`` routines now accept keyword arguments
- SAINT2 and ROSETTA distance restraints can now be written, ``format`` keywords are ``saint2`` and ``rosetta``
- ``StructureSelector`` added to score protein structures by contact satisfaction

**[0.10.2]**

- ``MANIFEST.ini`` file required by PyPi

**[0.10.1]**

- Critical bug fix in installation procedure and Cython-code compilation

**[0.10]**

*Added*

- Support for Python 3.7
- ``Cython`` added as dependency and ``SciPy`` removed
- ``conkit.misc.deprecate`` decorator for easier tagging
- ``ContactMap.match`` provides keyword to ``add_false_negatives`` found in the reference but not in contact map
- ``ContactMap.remove_false_negatives`` allows convenient removal of false negatives
- ``ContactMap.recall`` to calculate the recall of a contact map
- ``SequenceFile.summary`` for quick alignment summaries
- ``A2mParser`` to read HH-suite A2M alignment files
- Automatic ``sphinx-apidoc`` generation for up-to-date index
- ``ClustalParser`` to read CLUSTAL formatted files

*Changed*

- ``SequenceFile.calculate_freq`` backend changed from ``numpy`` to ``Cython`` for faster computation
- ``SequenceFile.calculate_weights`` backend changed from ``numpy`` to ``Cython`` for faster computation
- ``SequenceFile.filter`` backend changed from ``numpy`` to ``Cython`` for faster computation
- ``SequenceFile.filter_gapped`` backend changed from ``numpy`` to ``Cython`` for faster computation
- ``SequenceFile.calculate_weights`` renamed to ``SequenceFile.get_weights``
- ``SequenceFile.compute_freq`` renamed to ``SequenceFile.get_frequency``
- ``ContactMap.singletons`` backend changed from ``numpy`` to ``Cython`` for faster computation
- ``Bandwidth`` backend changed from ``numpy`` to ``Cython`` for faster computation
- ``ContactMap.short_range_contacts`` renamed to ``ContactMap.short_range``
- ``ContactMap.medium_range_contacts`` renamed to ``ContactMap.medium_range``
- ``ContactMap.long_range_contacts`` renamed to ``ContactMap.long_range``
- ``ContactMap.calculate_scalar_score`` renamed to ``ContactMap.set_scalar_score``
- ``ContactMap.calculate_contact_density`` renamed to ``ContactMap.get_contact_density``
- ``ContactMap.calculate_jaccard_index`` renamed to ``ContactMap.get_jaccard_index``
- ``ContactMatchState`` provides options for true positive, true negative, false positive and false negative, which can be added to contacts in the map at will
- ``Contact.is_match`` and ``Contact.define_match`` renamed to attribute ``Contact.true_positive``
- ``Contact.is_mismatch`` and ``Contact.define_mismatch`` renamed to attribute ``Contact.false_positive``
- ``Contact.is_unknown`` and ``Contact.define_unknown`` renamed to attribute ``Contact.status_unknown``
- ``Entity``, ``Gap`` and ``Residue`` classes made public

*Fixed*

- Bug fix in ``SequenceFile.filter`` to remove ``Sequence`` entries reliably
- Bug fix in ``ContactMapMatrixFigure`` when ``gap`` variable was less than 1

*Removed*

- Python 3.4 support

**[0.9]**

*Added*

- ``conkit.plot`` subpackage refactored to allow ``matplotlib`` access of ``Figure`` instances. This provides
  functionality similar to ``seaborn``, so ``matplotlib.Axes`` can be provided into which a plot is drawn.
- ``ContactMap.as_list`` function to represent the contact map as a 2D-list of residue indexes
- ``conkit.misc.normalize`` function to apply Feature scaling normalization
- ``CONTRIB.rst`` file to list all contributors
- ``SequenceFile.diversity`` property defined by :math:`\sqrt{N}/L`
- ``ContactMap.reindex`` to reindex a contact map given a new starting index
- ``ContactMap.singletons`` returns a copy of the contact map with singleton contacts, i.e. ones without neighbors
- ``Sequence.seq_encoded`` to allow turning a sequence into an encoded list
- ``Sequence.encoded_matrix`` to give the entire alignment as encoded matrix
- ``SequenceFile.filter_gapped`` to filter sequences with a certain threshold of gaps
- ``SequenceFile.to_string`` and ``ContactMap.to_string`` methods
- ``ContactMapMatrixFigure`` added to illustrate prediction signal of entire ``ContactMap``
- Added support for ``nebcon`` contact prediction format

*Changed*

- Changed API interface for ``conkit.plot`` in accordance to necessary changes for above
- ``ContactMapFigure`` now accepts ``lim`` parameters for axes limits
- ``ContactMapFigure`` and ``ContacctMapChordFigure`` improved to better space marker size
- Typos corrected in documentation 
- ``THREE_TO_ONE`` and ``ONE_TO_THREE`` dictionaries modified to ``Enum`` objects
- ``SequeneFile.neff`` renamed to ``SequenceFile.meff``
- ``ContactMapChordFigure.get_radius_around_circle`` moved to ``conkit.plot.tools.radius_around_circle``
- ``AmiseBW.curvature`` renamed to ``AmiseBW.gauss_curvature``

*Fixed*

- ``A3mParser`` keyword argument mismatch sorted

**[0.8.4]**

*Added*

- ``Entity.top`` property to always return the first child in the list
- ``ContactMap.find`` function accepts ``strict`` keyword argument to find contact pairs with both residues in ``register``
- ``PdbParser`` takes a distance cutoff of ``0`` to include all Cb-Cb contacts in the protein structure
- ``ContactMatchState`` enumerated type for definitions of state constants for contact
- ``SequenceAlignmentState`` enumerated type for definitions of state constants for each sequence file 
- ``NcontParser`` added to extract contact pairs identified by NCONT (CCP4 Software Suite) 

*Changed*

- Optimized some functions and comparisons according to the recommended Python optimization instructions 
- ``ContactMap.match`` does __not__ modifiy ``other`` by default anymore. Specify ``match_other=True`` as kwarg!
- ``ContactMap.calculate_kernel_density`` renamed to ``ContactMap.calculate_contact_density`` 
- ``ContactDensityFigure`` draws domain boundary lines instead of symbols

**[0.8.3]**

*Added*

- ``requirements.txt`` file re-added for easier dependency installation
- ``LinearBW`` calculator added for linear bandwidth calculation in analysis
- ``seq_ascii`` property to ``Sequence`` for encoded sequence
- ``ascii_matrix`` property to ``SequenceFile`` for encoded alignment 
- ``SequenceFile`` and ``ContactFile`` classes have new ``empty`` properties
- ``flib`` format for ``ContactFile`` classes to allow easier conversions for the Flib-Coevo fragment picking library

*Changed*

- Distance definitions accept floating point values
- ``_BandwidthCalc`` class renamed to ``BandwidthBase``
- Abstractified ``BandwidthBase``, and ``Parser`` with all subparser classes 
- Refactored ``conkit/io/__init__.py`` to avoid duplication of code

*Fixed*

- ``PconsParser`` class accepts negative ``raw_score`` values
- ``SequenceFile.neff`` returns ``float`` instead of ``int``
- ``CCMpredParser.read()`` returns empty ``ContactFile`` when matrix file empty

**[0.8.2]**

Added*

- Test function skipping added for ``SequenceFile.filter()`` when SciPy not installed

*Changed*

- Renamed conkit/io/tests files for filenames to agree with modules in conkit/io
- Performance of ``write()`` in parsers improved by construction of string and single call to ``write()`` of filehandle

*Fixed*

- Critical bug fix for automated opening of filehandle in Python2.7 

**[0.8.1]**

*Changed*

- Revoked catching of ``SystemExit(0)`` exception in scripts when invoked with ``--help`` flag 

*Fixed*

- Bug fix relating to Python3 automatic opening of file handles - Thanks to Miguel Correa for reporting this bug

**[0.8]**

*Added*

- Logging message coloring according to message level
- ``filter()`` function added for redundancy/distant homolog removal from ``SequenceFile``
- License text added to each module
- ``io`` sub-package caches modules and imports upon request

*Changed*

- Default value in ``calculate_meff()`` and ``calculate_weights()`` changed from 0.7 to 0.8 [more commonly used in literature]
- ``core`` classes extracted to individual module files

*Fixed*

- Bug fix with PyPi installation where ``requirements.txt`` not found; fix includes removal of ``requirements.txt`` and addition of ``install_requires`` to ``setup.py`` instead. - Thanks to Miguel Correa for reporting this bug
