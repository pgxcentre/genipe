
.. contents:: Quick navigation
   :depth: 2


Impute2 Extractor - ``impute2-extractor``
==========================================


Extract imputed markers located in a specific genomic region.


General options
----------------

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
------------

.. table::

    +--------------------+----------------------------------------------------+
    | Option             | Description                                        |
    +====================+====================================================+
    | ``--impute2 FILE`` | The output from IMPUTE2.                           |
    +--------------------+----------------------------------------------------+


Indexation options
-------------------

.. table::

    +-------------+-----------------------------------------------------------+
    | Option      | Description                                               |
    +=============+===========================================================+
    | ``--index`` | Only perform the indexation.                              |
    +-------------+-----------------------------------------------------------+


Output options
---------------

.. table::

    +----------------------------------+--------------------------------------+
    | Option                           | Description                          |
    +==================================+======================================+
    | ``--out PREFIX``                 | The prefix of the output files.      |
    |                                  | [``impute2_extractor``]              |
    +----------------------------------+--------------------------------------+
    | ``--format FORMAT [FORMAT ...]`` | The output format. Can specify either|
    |                                  | 'impute2' for probabilities (same as |
    |                                  | impute2 format, i.e. 3 values per    |
    |                                  | sample), 'dosage' for dosage values  |
    |                                  | (one value between 0 and 2 by        |
    |                                  | sample), 'calls' for hard calls, or  |
    |                                  | 'bed' for Plink binary format (with  |
    |                                  | hard calls). [``impute2``]           |
    +----------------------------------+--------------------------------------+
    | ``--long``                       | Write the output file in the long    |
    |                                  | format (one line per sample per      |
    |                                  | marker). This option is only         |
    |                                  | compatible with the 'calls' and      |
    |                                  | 'dosage' format (option '-- format').|
    +----------------------------------+--------------------------------------+
    | ``--prob FLOAT``                 | The probability threshold used when  |
    |                                  | creating a file in the dosage or call|
    |                                  | format. [``0.9``]                    |
    +----------------------------------+--------------------------------------+


Extraction options
-------------------

.. table::

    +-----------------------------+-------------------------------------------+
    | Option                      | Description                               |
    +=============================+===========================================+
    | ``--extract FILE``          | File containing marker names to extract.  |
    +-----------------------------+-------------------------------------------+
    | ``--genomic CHR:START-END`` | The range to extract (*e.g.* 22 1000000   |
    |                             | 1500000). Can be use in combination with  |
    |                             | '``--rate``', '``--maf``' and             |
    |                             | '``--info``'.                             |
    +-----------------------------+-------------------------------------------+
    | ``--maf FLOAT``             | Extract markers with a minor allele       |
    |                             | frequency equal or higher than the        |
    |                             | specified threshold. Can be use in        |
    |                             | combination with '--rate', '--info' and   |
    |                             | '--genomic'.                              |
    +-----------------------------+-------------------------------------------+
    | ``--rate FLOAT``            | Extract markers with a completion rate    |
    |                             | equal or higher to the specified          |
    |                             | threshold. Can be use in combination with |
    |                             | '``--maf``', '``--info``' and             |
    |                             | '``--genomic``'.                          |
    +-----------------------------+-------------------------------------------+
    | ``--info FLOAT``            | Extract markers with an information equal |
    |                             | or higher to the specified threshold. Can |
    |                             | be use in combination with '``--maf``',   |
    |                             | '``--rate``' and '``--genomic``'.         |
    +-----------------------------+-------------------------------------------+

