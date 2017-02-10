
.. contents:: Quick navigation
   :depth: 3


Imputed Stats - Options and Parameters
=======================================


Available statistical models
-----------------------------

.. table::

    +--------------+----------------------------------------------------------+
    | Name         | Description                                              |
    +==============+==========================================================+
    | ``cox``      | Cox's proportional hazard model (survival regression).   |
    +--------------+----------------------------------------------------------+
    | ``linear``   | Linear regression (ordinary least squares).              |
    +--------------+----------------------------------------------------------+
    | ``logistic`` | Logistic regression (GLM with binomial distribution).    |
    +--------------+----------------------------------------------------------+
    | ``mixedlm``  | Linear mixed effect model (random intercept).            |
    +--------------+----------------------------------------------------------+
    | ``skat``     | SKAT analysis.                                           |
    +--------------+----------------------------------------------------------+


Common options
---------------


General options
^^^^^^^^^^^^^^^^

.. table::

    +-----------------------+-------------------------------------------------+
    | Option                | Description                                     |
    +=======================+=================================================+
    | ``-h``, ``--help``    | Show this help message and exit.                |
    +-----------------------+-------------------------------------------------+
    | ``-v``, ``--version`` | Show program's version number and exit.         |
    +-----------------------+-------------------------------------------------+
    | ``--debug``           | Set the logging level to debug.                 |
    +-----------------------+-------------------------------------------------+


Input files
^^^^^^^^^^^^

.. table::

    +--------------------------+----------------------------------------------+
    | Option                   | Description                                  |
    +==========================+==============================================+
    | ``--impute2 FILE``       | The output from IMPUTE2.                     |
    +--------------------------+----------------------------------------------+
    | ``--sample FILE``        | The sample file (the order should be the same|
    |                          | as in the IMPUTE2 files).                    |
    +--------------------------+----------------------------------------------+
    | ``--pheno FILE``         | The file containing phenotypes and co        |
    |                          | variables.                                   |
    +--------------------------+----------------------------------------------+
    | ``--extract-sites FILE`` | A list of sites to extract for analysis      |
    |                          | (optional).                                  |
    +--------------------------+----------------------------------------------+


Output options
^^^^^^^^^^^^^^^

.. table::

    +----------------+--------------------------------------------------------+
    | Option         | Description                                            |
    +================+========================================================+
    | ``--out FILE`` | The prefix for the output files. [``imputed_stats``]   |
    +----------------+--------------------------------------------------------+


General options
^^^^^^^^^^^^^^^^

.. table::

    +--------------------------+----------------------------------------------+
    | Option                   | Description                                  |
    +==========================+==============================================+
    | ``--nb-process INT``     | The number of process to use. [``1``]        |
    +--------------------------+----------------------------------------------+
    | ``--nb-lines INT``       | The number of line to read at a time.        |
    |                          | [``1000``]                                   |
    +--------------------------+----------------------------------------------+
    | ``--chrx``               | The analysis is performed for the non        |
    |                          | pseudo-autosomal region of the chromosome X  |
    |                          | (male dosage will be divided by 2 to get     |
    |                          | values [0, 0.5] instead of [0, 1]) (males are|
    |                          | coded as 1 and option '``--gender-column``'  |
    |                          | should be used).                             |
    +--------------------------+----------------------------------------------+
    | ``--gender-column NAME`` | The name of the gender column (use to exclude|
    |                          | samples with unknown gender (*i.e.* not 1,   |
    |                          | male, or 2, female). If gender not available,|
    |                          | use 'None'. [``Gender``]                     |
    +--------------------------+----------------------------------------------+


Dosage options
^^^^^^^^^^^^^^^

.. table::

    +------------------+------------------------------------------------------+
    | Option           | Description                                          |
    +==================+======================================================+
    | ``--scale INT``  | Scale dosage so that values are in [0, n] (possible  |
    |                  | values are 1 (no scaling) or 2). [``2``]             |
    +------------------+------------------------------------------------------+
    | ``--prob FLOAT`` | The minimal probability for which a genotype should  |
    |                  | be considered. [``>=0.9``]                           |
    +------------------+------------------------------------------------------+
    | ``--maf FLOAT``  | Minor allele frequency threshold for which marker    |
    |                  | will be skipped. [``<0.01``]                         |
    +------------------+------------------------------------------------------+


Phenotype options
^^^^^^^^^^^^^^^^^^

.. table::

    +--------------------------+----------------------------------------------+
    | Option                   | Description                                  |
    +==========================+==============================================+
    | ``--covar NAME``         | The co variable names (in the phenotype      |
    |                          | file), separated by coma.                    |
    +--------------------------+----------------------------------------------+
    | ``--categorical NAME``   | The name of the variables that are           |
    |                          | categorical (note that the gender is always  |
    |                          | categorical). The variables are separated by |
    |                          | coma.                                        |
    +--------------------------+----------------------------------------------+
    | ``--missing-value NAME`` | The missing value in the phenotype file.     |
    +--------------------------+----------------------------------------------+
    | ``--sample-column NAME`` | The name of the sample ID column (in the     |
    |                          | phenotype file). [``sample_id``]             |
    +--------------------------+----------------------------------------------+
    | ``--interaction NAME``   | Add an interaction between the genotype and  |
    |                          | this variable.                               |
    +--------------------------+----------------------------------------------+


Cox's proportional hazard model options
----------------------------------------

.. table::

    +--------------------------+----------------------------------------------+
    | Option                   | Description                                  |
    +==========================+==============================================+
    | ``--time-to-event NAME`` | The time to event variable (in the pheno     |
    |                          | file).                                       |
    +--------------------------+----------------------------------------------+
    | ``--event NAME``         | The event variable (1 if observed, 0 if not  |
    |                          | observed).                                   |
    +--------------------------+----------------------------------------------+


Linear regression options
--------------------------

.. table::

    +-----------------------+-------------------------------------------------+
    | Option                | Description                                     |
    +=======================+=================================================+
    | ``--pheno-name NAME`` | The phenotype.                                  |
    +-----------------------+-------------------------------------------------+


Logistic regression options
----------------------------

.. table::

    +-----------------------+-------------------------------------------------+
    | Option                | Description                                     |
    +=======================+=================================================+
    | ``--pheno-name NAME`` | The phenotype.                                  |
    +-----------------------+-------------------------------------------------+


Linear mixed effects options
-----------------------------

.. table::

    +-------------------------+-----------------------------------------------+
    | Option                  | Description                                   |
    +=========================+===============================================+
    | ``--pheno-name NAME``   | The phenotype.                                |
    +-------------------------+-----------------------------------------------+
    | ``--use-ml``            | Fit the standard likelihood using maximum     |
    |                         | likelihood (ML) estimation instead of REML    |
    |                         | (default is REML).                            |
    +-------------------------+-----------------------------------------------+
    | ``--p-threshold FLOAT`` | The p-value threshold for which the real      |
    |                         | MixedLM analysis will be performed.           |
    |                         | [``<0.0001``]                                 |
    +-------------------------+-----------------------------------------------+


SKAT options
-------------

.. table::

    +------------------------------------------+------------------------------+
    | Option                                   | Description                  |
    +==========================================+==============================+
    | ``--snp-sets FILE``                      | A file indicating a snp_set  |
    |                                          | and an optional weight for   |
    |                                          | every variant.               |
    +------------------------------------------+------------------------------+
    | ``--outcome-type {continuous,discrete}`` | The variable type for the    |
    |                                          | outcome. This will be passed |
    |                                          | to SKAT.                     |
    +------------------------------------------+------------------------------+
    | ``--skat-o``                             | By default, the regular SKAT |
    |                                          | is used. Setting this flag   |
    |                                          | will use the SKAT-O algorithm|
    |                                          | instead.                     |
    +------------------------------------------+------------------------------+
    | ``--pheno-name NAME``                    | The phenotype.               |
    +------------------------------------------+------------------------------+

