
.. contents:: Quick navigation
   :depth: 2


Main Pipeline - Options and Parameters
=======================================


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
    | ``--thread THREAD``   | Number of threads [``1``].                      |
    +-----------------------+-------------------------------------------------+


Input options
--------------

.. table::

    +----------------------+--------------------------------------------------+
    | Option               | Description                                      |
    +======================+==================================================+
    | ``--bfile PREFIX``   | The prefix of the binary pedfiles (input data).  |
    +----------------------+--------------------------------------------------+
    | ``--reference FILE`` | The human reference to perform an initial strand |
    |                      | check (useful for genotyped markers not in the   |
    |                      | IMPUTE2 reference files) (optional).             |
    +----------------------+--------------------------------------------------+


Output options
---------------

.. table::

    +-------------------------------+-----------------------------------------+
    | Option                        | Description                             |
    +===============================+=========================================+
    | ``--chrom CHROM [CHROM ...]`` | The chromosomes to process.             |
    +-------------------------------+-----------------------------------------+
    | ``--output-dir DIR``          | The name of the output directory.       |
    |                               | [``genipe``]                            |
    +-------------------------------+-----------------------------------------+
    | ``--bgzip``                   | Use bgzip to compress the impute2 files.|
    +-------------------------------+-----------------------------------------+


HPC options
------------

.. table::

    +-------------------------+-----------------------------------------------+
    | Option                  | Description                                   |
    +=========================+===============================================+
    | ``--use-drmaa``         | Launch tasks using DRMAA.                     |
    +-------------------------+-----------------------------------------------+
    | ``--drmaa-config FILE`` | The configuration file for tasks (use this    |
    |                         | option when launching tasks using DRMAA). This|
    |                         | file should describe the walltime and the     |
    |                         | number of nodes/processors to use for each    |
    |                         | task.                                         |
    +-------------------------+-----------------------------------------------+
    | ``--preamble FILE``     | This option should be used when using DRMAA on|
    |                         | a HPC to load required module and set         |
    |                         | environment variables. The content of the file|
    |                         | will be added between the 'shebang' line and  |
    |                         | the tool command.                             |
    +-------------------------+-----------------------------------------------+


SHAPEIT options
----------------

.. table::

    +-----------------------------+-------------------------------------------+
    | Option                      | Description                               |
    +=============================+===========================================+
    | ``--shapeit-bin BINARY``    | The SHAPEIT binary if it's not in the     |
    |                             | path.                                     |
    +-----------------------------+-------------------------------------------+
    | ``--shapeit-thread INT``    | The number of thread for phasing. [``1``] |
    +-----------------------------+-------------------------------------------+
    | ``--shapeit-extra OPTIONS`` | SHAPEIT extra parameters. Put extra       |
    |                             | parameters between single or normal quotes|
    |                             | (*e.g.* ``--shapeit-extra '-- states 100  |
    |                             | --window 2'``).                           |
    +-----------------------------+-------------------------------------------+


Plink options
--------------

.. table::

    +------------------------+------------------------------------------------+
    | Option                 | Description                                    |
    +========================+================================================+
    | ``--plink-bin BINARY`` | The Plink binary if it's not in the path.      |
    +------------------------+------------------------------------------------+


IMPUTE2 autosomal reference
----------------------------

.. table::

    +--------------------------------+----------------------------------------+
    | Option                         | Description                            |
    +================================+========================================+
    | ``--hap-template TEMPLATE``    | The template for IMPUTE2's haplotype   |
    |                                | files (replace the chromosome number by|
    |                                | '``{chrom}``', *e.g.*                  |
    |                                | '``1000GP_{chrom}.hap.gz``').          |
    +--------------------------------+----------------------------------------+
    | ``--legend-template TEMPLATE`` | The template for IMPUTE2's legend files|
    |                                | (replace the chromosome number by      |
    |                                | '``{chrom}``', *e.g.*                  |
    |                                | '``1000GP_{chrom}.legend.gz``').       |
    +--------------------------------+----------------------------------------+
    | ``--map-template TEMPLATE``    | The template for IMPUTE2's map files   |
    |                                | (replace the chromosome number by      |
    |                                | '``{chrom}``', *e.g.*                  |
    |                                | '``genetic_map_chr{chrom}.txt``').     |
    +--------------------------------+----------------------------------------+
    | ``--sample-file FILE``         | The name of IMPUTE2's sample file.     |
    +--------------------------------+----------------------------------------+


IMPUTE2 chromosome X reference
-------------------------------

.. table::

    +--------------------------+----------------------------------------------+
    | Option                   | Description                                  |
    +==========================+==============================================+
    | ``--hap-nonPAR FILE``    | The IMPUTE2's haplotype file for the         |
    |                          | non-pseudoautosomal region of chromosome 23. |
    +--------------------------+----------------------------------------------+
    | ``--hap-PAR1 FILE``      | The IMPUTE2's haplotype file for the first   |
    |                          | pseudoautosomal region of chromosome 23.     |
    +--------------------------+----------------------------------------------+
    | ``--hap-PAR2 FILE``      | The IMPUTE2's haplotype file for the second  |
    |                          | pseudoautosomal region of chromosome 23.     |
    +--------------------------+----------------------------------------------+
    | ``--legend-nonPAR FILE`` | The IMPUTE2's legend file for the            |
    |                          | non-pseudoautosomal region of chromosome 23. |
    +--------------------------+----------------------------------------------+
    | ``--legend-PAR1 FILE``   | The IMPUTE2's legend file for the first      |
    |                          | pseudoautosomal region of chromosome 23.     |
    +--------------------------+----------------------------------------------+
    | ``--legend-PAR2 FILE``   | The IMPUTE2's legend file for the second     |
    |                          | pseudoautosomal region of chromosome 23.     |
    +--------------------------+----------------------------------------------+
    | ``--map-nonPAR FILE``    | The IMPUTE2's map file for the               |
    |                          | non-pseudoautosomal region of chromosome 23. |
    +--------------------------+----------------------------------------------+
    | ``--map-PAR1 FILE``      | The IMPUTE2's map file for the first         |
    |                          | pseudoautosomal region of chromosome 23.     |
    +--------------------------+----------------------------------------------+
    | ``--map-PAR2 FILE``      | The IMPUTE2's map file for the second        |
    |                          | pseudoautosomal region of chromosome 23.     |
    +--------------------------+----------------------------------------------+


IMPUTE2 options
----------------

.. table::

    +---------------------------------------+---------------------------------+
    | Option                                | Description                     |
    +=======================================+=================================+
    | ``--impute2-bin BINARY``              | The IMPUTE2 binary if it's not  |
    |                                       | in the path.                    |
    +---------------------------------------+---------------------------------+
    | ``--segment-length BP``               | The length of a single segment  |
    |                                       | for imputation. [``5e+06``]     |
    +---------------------------------------+---------------------------------+
    | ``--filtering-rules RULE [RULE ...]`` | IMPUTE2 filtering rules         |
    |                                       | (optional).                     |
    +---------------------------------------+---------------------------------+
    | ``--impute2-extra OPTIONS``           | IMPUTE2 extra parameters. Put   |
    |                                       | the extra parameters between    |
    |                                       | single or normal quotes (*e.g.* |
    |                                       | ``--impute2-extra '-buffer 250  |
    |                                       | -Ne 20000'``).                  |
    +---------------------------------------+---------------------------------+


IMPUTE2 merger options
-----------------------

.. table::

    +-------------------------+-----------------------------------------------+
    | Option                  | Description                                   |
    +=========================+===============================================+
    | ``--probability FLOAT`` | The probability threshold for no calls.       |
    |                         | [``<0.9``]                                    |
    +-------------------------+-----------------------------------------------+
    | ``--completion FLOAT``  | The completion rate threshold for site        |
    |                         | exclusion. [``<0.98``]                        |
    +-------------------------+-----------------------------------------------+
    | ``--info FLOAT``        | The measure of the observed statistical       |
    |                         | information associated with the allele        |
    |                         | frequency estimate threshold for site         |
    |                         | exclusion. [``<0.00``]                        |
    +-------------------------+-----------------------------------------------+


Automatic report options
-------------------------

.. table::

    +------------------------------------+------------------------------------+
    | Option                             | Description                        |
    +====================================+====================================+
    | ``--report-number NB``             | The report number.                 |
    |                                    | [``genipe automatic report``]      |
    +------------------------------------+------------------------------------+
    | ``--report-title TITLE``           | The report title. [``genipe:       |
    |                                    | Automatic genome-wide imputation``]|
    +------------------------------------+------------------------------------+
    | ``--report-author AUTHOR``         | The report author. [``Automatically|
    |                                    | generated by genipe``]             |
    +------------------------------------+------------------------------------+
    | ``--report-background BACKGROUND`` | The report background section (can |
    |                                    | either be a string or a file       |
    |                                    | containing the background.         |
    |                                    | [``General background``]           |
    +------------------------------------+------------------------------------+

