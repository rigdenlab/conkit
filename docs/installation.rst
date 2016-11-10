.. _installation:

Installation
============

PyPi
----
The easiest way to install ConKit is via ``easy_install`` or ``pip``. To do so, simply run the following command on your command line and your are ready to go.

.. code-block:: bash

   $> pip install conkit


Source Code
-----------

If you want to install ConKit from source, download the `latest version`_ from the GitHub repository. Then, use the following commands to install ConKit.

Download the latest version directly from GitHub or via your command line using ``git clone``

.. code-block:: bash

   $> git clone https://github.com/fsimkovic/conkit conkit
   $> cd conkit

One downloaded, you might want to check that all functions are available on your system. Run the following command:

.. code-block:: bash

   $> python setup.py test

If this has completed successfully, you are good to go and you can now install ConKit.

.. code-block:: bash

   $> sudo python setup.py build install

ConKit is now successfully installed in your Python's site-packages directory.

.. note::

   The ConKit project has recently been added to the CCP4 and CCP-EM devtools distributions. This means that they will most likely be shipped with the full releases in the near future.


.. _latest version: https://github.com/fsimkovic/conkit/releases
