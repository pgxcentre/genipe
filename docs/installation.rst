Installation
=============


Quick navigation
-----------------

1. :ref:`install-requirements`
2. :ref:`install-virt`
3. :ref:`install-test`
4. :ref:`install-update`


.. _install-requirements:

Requirements
-------------

The tools require a standard Python 3 installation with the following modules:

* ``numpy`` version 1.8.2 or latest
* ``jinja2`` version 2.7.3 or latest
* ``pandas`` version 0.15.2 or latest
* ``setuptools`` version 12.0.5 or latest

The following modules are optional (required for statistical analysis and data
management):

* ``matplotlib`` version 1.4.2
* ``scipy`` version 0.15.1 or latest (required by ``statsmodels``)
* ``statsmodels`` version 0.6.1 or latest
* ``lifelines`` version 0.7.0 or latest
* ``biopython`` version 1.65 or latest
* ``pyfaidx`` version 0.3.7 or latest

.. note::

   Only Python versions 3.3 and higher are supported for now.


.. _install-virt:

Virtual environment
--------------------

Using pyvenv
^^^^^^^^^^^^^

The following commands should successfully create a virtual environment and
activate it, as long as Python 3 was previously installed on the machine.

.. code-block:: bash

   # Creating the environment
   pyvenv $HOME/gwip_pyvenv

   # Activating the environment
   source $HOME/gwip_pyvenv/bin/activate

Then, install the :py:mod:`gwip` package (and all its dependencies) using the
following command:

.. code-block:: bash

   pip install gwip

If required, the optional dependencies can be installed using the following
command:

.. code-block:: bash

   pip install statsmodels
   pip install lifelines
   pip install biopython
   pip install pyfaidx
   pip install matplotlib


Using Miniconda
^^^^^^^^^^^^^^^^

The following commands should successfully download and install Miniconda,
create and activate a new python virtual environment. This installation method
doesn't required a previous Python 3 installation.

.. code-block:: bash

   # Installing Miniconda
   wget http://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
   bash miniconda.sh -b -p $HOME/miniconda

   # Creating the environment
   $HOME/miniconda/bin/conda create -q -n gwip_pyvenv python=3.4

   # Activating the environment
   source $HOME/miniconda/bin/activate gwip_pyvenv

Then, install the :py:mod:`gwip` package (and all its dependencies) using the
following command:

.. code-block:: bash

   conda install gwip -c http://statgen.org/wp-content/uploads/Softwares/gwip

If required, the optional dependencies can be installed using the following
command:

.. code-block:: bash

   conda install -y scipy
   conda install -y statsmodels
   conda install -y biopython
   conda install -y matplotlib
   pip install --no-deps pyfaidx
   pip install --no-deps lifelines


.. _gwip-pyvenv-activation:

Virtual environment activation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Before using the :py:mod:`gwip` module for any analysis, the Python virtual
environment needs to be activated. Depending of your installation type (*i.e*
*pyvenv* or *miniconda*), the activation process will differ.


Pyvenv
"""""""

If the module was installed into a *pyvenv* environment, perform the following
command to activate it.

.. code-block:: bash

   source $HOME/gwip_pyvenv/bin/activate


Miniconda
""""""""""

If the module was installed into a *Miniconda* environment, perform the
following command to activate it.

.. code-block:: bash

   source $HOME/miniconda/bin/activate gwip_pyvenv


.. _install-test:

Testing the installation
-------------------------

To test the installation, make sure that the virtual environment is activated.
Then, launch python and use the following commands:

.. code-block:: python

   >>> import gwip
   >>> gwip.test()
   ......................ss.ss.......................ss...ss...s.s.........
   ----------------------------------------------------------------------
   Ran 72 tests in 107.268s
   
   OK (skipped=10)


.. _install-update:

Updating the package
---------------------

If there is a new :py:mod:`gwip` release, perform one of the following command
(depending of the installation method). Don't forget to first activate the
python virtual environment.


Pyvenv
^^^^^^^

.. code-block:: bash

   pip install -U gwip


Miniconda
^^^^^^^^^^

.. code-block:: bash

   conda update gwip -c http://statgen.org/wp-content/uploads/Softwares/gwip

