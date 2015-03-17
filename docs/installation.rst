Installation
============

The tools require a standard Python 3 installation with the following modules:

* ``numpy`` version 1.8.2 or latest
* ``jinja2`` version 2.7.3 or latest
* ``pandas`` version 0.15.2 or latest
* ``setuptools`` version 12.0.5 or latest
* ``pyfaidx`` version 0.3.7 or latest

The following modules are optional (required for statistical analysis):

* ``lifelines`` version 0.7.0 or latest
* ``statsmodels`` version 0.6.1 or latest

.. note::

   Only Python versions 3.3 and higher are supported for now.


Download
--------

Since the module is still in development, it is recommended to simply clone the
directory directly from `github <https://github.com/pgxcentre/gwip>`_. Once
official releases are available, they will be located
`here <https://github.com/pgxcentre/gwip/releases>`_.

.. code-block:: console

   $ git clone https://github.com/pgxcentre/gwip.git $HOME/gwip

For simplicity of use, we recommend to install the :py:mod:`gwip` module using
a Python virtual environment. There are two possible ways to create such an
environment: using Python's ``pyvenv`` command, or using
`Miniconda <http://conda.pydata.org/miniconda.html>`_. The latter is easier to
use since it doesn't require any compilation.


Using pyvenv
------------

The following commands should successfully create a virtual environment and
activate it, as long as Python 3 was previously installed on the machine.

.. code-block:: console

   $ pyvenv $HOME/gwip_pyvenv
   $ source $HOME/gwip_pyvenv/bin/activate


Then, install the package (and all its dependencies) using the following
command:

.. code-block:: console

   $ pip install $HOME/gwip


Using Miniconda
---------------

The following commands should successfully download and install Miniconda,
create and activate a new python virtual environment. This installation method
doesn't required a previous Python 3 installation.

.. code-block:: console

   $ wget http://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
   $ bash miniconda.sh -b -p $HOME/miniconda
   $ $HOME/miniconda/bin/conda create -q -n gwip_pyvenv python=3.4
   $ source $HOME/miniconda/bin/activate gwip_pyvenv

It is recommended to install the dependencies using ``conda`` instead of
``pip`` when available (before installing :py:mod:`gwip`).

.. code-block:: console

   $ conda install numpy
   $ conda install jinja2
   $ conda install pandas
   $ conda install matplotlib
   $ conda install setuptools
   $ conda install statsmodels

Finally, using ``pip``, install the missing dependencies and :py:mod:`gwip`:

.. code-block:: console

   $ pip install --no-deps pyfaidx
   $ pip install --no-deps lifelines
   $ pip install --no-deps $HOME/gwip


Testing the installation
------------------------

To test the installation, make sure that the virtual environment is activated.

.. code-block:: console
   
   $ python --version
   Python 3.4.3

Then, launch python and use the following commands:

.. code-block:: python

    >>> import gwip
    >>> gwip.test()
    .....................ss.ss........................ss...s.s.........
    ----------------------------------------------------------------------
    Ran 67 tests in 76.401s

    OK (skipped=8)

