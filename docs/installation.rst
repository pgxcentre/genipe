Installation
=============

.. contents:: Quick navigation


.. _install-requirements:

Requirements
-------------

The tools require a standard Python 3.4 (or latest) installation with the
following modules:

* ``numpy`` version 1.8.2 or latest
* ``jinja2`` version 2.7.3 or latest
* ``pandas`` version 0.17.0 or latest
* ``setuptools`` version 12.0.5 or latest

The following modules are optional (required for statistical analysis and data
management):

* ``matplotlib`` version 1.4.2
* ``scipy`` version 0.15.1 or latest (required by ``statsmodels``)
* ``statsmodels`` version 0.6.1 or latest
* ``lifelines`` version 0.7.0 or latest
* ``biopython`` version 1.65 or latest
* ``pyfaidx`` version 0.3.7 or latest
* ``drmaa`` version 0.7.6 or latest

.. note::

   Only Python version 3.4 and higher are supported.


.. _install-virt:

Installing in a virtual environment
------------------------------------


.. _install-pyvenv:

Using python's pyvenv
^^^^^^^^^^^^^^^^^^^^^^

The following commands should successfully create a virtual environment and
activate it, as long as Python 3.4 (or latest) was previously installed on the
machine.

.. code-block:: bash

   # Creating the environment
   pyvenv $HOME/genipe_pyvenv

   # Activating the environment
   source $HOME/genipe_pyvenv/bin/activate

Then, install the :py:mod:`genipe` package (and all its dependencies) using the
following command:

.. code-block:: bash

   pip install genipe

If required, the optional dependencies can be installed using the following
command:

.. code-block:: bash

   pip install scipy
   pip install statsmodels
   pip install lifelines
   pip install biopython
   pip install pyfaidx
   pip install matplotlib
   pip install drmaa


.. _install-miniconda:

Using Miniconda
^^^^^^^^^^^^^^^^

The following commands should successfully download and install Miniconda,
create and activate a new python virtual environment. This installation method
doesn't require a previous Python 3 installation.

.. code-block:: bash

   # Installing Miniconda
   wget http://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
   bash miniconda.sh -b -p $HOME/miniconda

   # Creating the environment
   $HOME/miniconda/bin/conda create -q -n genipe_pyvenv python=3

   # Activating the environment
   source $HOME/miniconda/bin/activate genipe_pyvenv

Then, install the :py:mod:`genipe` package (and all its dependencies) using the
following command:

.. code-block:: bash

   conda install genipe -c http://statgen.org/wp-content/uploads/Softwares/genipe

.. note::

   It is possible to add the channel to conda's configuration (so that you
   won't need to use the ``-c`` option for installing or updating). To do so,
   perform the following command:

   .. code-block:: bash

      conda config --add channels http://statgen.org/wp-content/uploads/Softwares/genipe

   Once this command is executed, you can always ommit
   ``-c http://statgen.org/...`` in the ``conda`` commands (for installing or
   updating).

If required, the optional dependencies can be installed using the following
command:

.. code-block:: bash

   conda install -y scipy
   conda install -y statsmodels
   conda install -y biopython
   conda install -y matplotlib
   conda install -y drmaa
   pip install --no-deps pyfaidx
   pip install --no-deps lifelines


.. _genipe-pyvenv-activation:

Virtual environment activation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Before using the :py:mod:`genipe` module for any analysis, the Python virtual
environment needs to be activated. Depending of your installation type (*i.e*
*pyvenv* or *miniconda*), the activation process will differ.


Pyvenv
"""""""

If the module was installed into a *pyvenv* environment, perform the following
command to activate it.

.. code-block:: bash

   source $HOME/genipe_pyvenv/bin/activate


Miniconda
""""""""""

If the module was installed into a *Miniconda* environment, perform the
following command to activate it.

.. code-block:: bash

   source $HOME/miniconda/bin/activate genipe_pyvenv


.. _install-test:

Testing the installation
-------------------------

The :py:mod:`genipe` module has been tested with the most recent versions of
the requirements. To test the installation, make sure that the virtual
environment is activated. Then, launch Python and use the following python
commands:

.. code-block:: python

   >>> import genipe
   >>> genipe.test()


.. _install-update:

Updating the package
---------------------

If there is a new :py:mod:`genipe` release, perform one of the following
commands (depending of the installation method). Don't forget to first activate
the python virtual environment.


Using python's pyvenv
^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   pip install -U genipe


Using Miniconda
^^^^^^^^^^^^^^^^

.. code-block:: bash

   conda update genipe -c http://statgen.org/wp-content/uploads/Softwares/genipe

.. note::

   If you have configured ``conda`` to use the :py:mod:`genipe` channel (see
   the note above), the following command can be executed to update the
   package:

   .. code-block:: bash

      conda update genipe

