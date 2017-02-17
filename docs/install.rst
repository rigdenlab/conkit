.. _installation:

Installation
============

This page contains information on how to install ConKit on your local machine or server. A number of different options exist; however, note that the installation via the :ref:`Python Package Index` remains to be the recommended and simplest.

Python Package Index
~~~~~~~~~~~~~~~~~~~~
The easiest way to install ConKit is via `easy_install` or `pip`. To do so, simply run the following command on your command line and you are ready to go.

To install using `easy_install`:

.. code-block:: bash

   $> easy_install conkit

To install using `pip`:

.. code-block:: bash

   $> pip install conkit

.. note::

   The executable scripts are automatically installed in your default ``bin`` directory. There is nothing more that you need to do.

Source Code
~~~~~~~~~~~

If you would like to install ConKit from source, download the `latest version <https://github.com/rigdenlab/conkit/releases>`_ from the GitHub repository. Then, use the following commands to install ConKit.

.. code-block:: bash

   $> git clone https://github.com/rigdenlab/conkit
   $> cd conkit

Once downloaded, you might want to check that all functions are available on your system. Run the following command:

.. code-block:: bash

   $> python setup.py test

If this has completed successfully, you are good to go and you can now install ConKit.

.. code-block:: bash

   $> sudo python setup.py build install

ConKit is now successfully installed in your system's default Python.

.. note::

   Similarly to the :ref:`Python Package Index` install, the executable scripts are automatically installed.

External software
~~~~~~~~~~~~~~~~~

.. note::

   If you install ConKit via the :ref:`Python Package Index`, the dependencies are automatically installed for you!

Required dependencies
+++++++++++++++++++++
Python 2.7, 3.4 or 3.5
  `Download Python <https://www.python.org/downloads/>`_
NumPy 1.8.2 (or later)
  `Download NumPy <http://www.scipy.org/scipylib/download.html>`_
BioPython 1.64 (or later)
  `Download BioPython <http://biopython.org/wiki/Documentation>`_
Matplotlib 1.3.1 (or later)
  `Download matplotlib <http://matplotlib.org/users/installing.html>`_
setuptools
  `Documentation <https://setuptools.readthedocs.io/en/latest/>`_    

Optional dependencies
+++++++++++++++++++++
SciPy 0.16 (or later)
  `Download SciPy <http://www.scipy.org/scipylib/download.html>`_
HHblits
   `Download HHblits <https://github.com/soedinglab/hh-suite>`_
HHblits Database
   `Download HHblits Database <http://wwwuser.gwdg.de/~compbiol/data/hhsuite/databases/hhsuite_dbs>`_
CCMpred
   `Download CCMpred <https://github.com/soedinglab/CCMpred>`_
