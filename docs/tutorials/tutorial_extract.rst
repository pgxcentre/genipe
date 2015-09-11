Site Extraction Tutorial
=========================


Quick navigation
-----------------

1. :ref:`extract-tut-input-files`
2. :ref:`extract-tut-execute`
3. :ref:`extract-tut-output-files`
4. :ref:`extract-tut-usage`


Site extraction
----------------

Genome-wide imputation dataset might be huge. Often, it is required to extract
a subset of imputed sites (*e.g.* specific markers, genomic location, or
markers with a specific minor allele frequency, information value or completion
rate). Also, different format might be required, depending of the underlying
analysis (*e.g.* hard calls or dosage values). We provide an easy tool to
perform site extraction of multiple *impute2* files using either marker
identification number, or genomic location and/or minor allele frequency and/or
call rate and/or information value.

We suppose that you have followed the main :ref:`genipe-tut-page`. The
following command will create the working directory for this tutorial.

.. code-block:: bash

   mkdir -p $HOME/genipe_tutorial/extraction


.. _extract-tut-input-files:

Input files
^^^^^^^^^^^^

After running the :py:mod:`genipe` pipeline, all the required files for the
extraction tools are automatically created in the ``final_impute2`` directories
(see the :ref:`genipe-tut-output-files-final_impute2` section in the main
:ref:`genipe-tut-page`).

The files that are required in these directories depends of what kind of
extraction is required (by name, or by genomic location and/or by minor allele
frequency and/or by calling rate and/or by information value).

Once the required *impute2* files are provided to the tool, the other required
files will be automatically fetched (if required).


.. _extract-tut-execute:

Executing the extraction
^^^^^^^^^^^^^^^^^^^^^^^^^

The first time the tool is used on a set of *impute2* files, indexation will
automatically occur (to speed of the analysis for future extraction). There are
two ways to extract markers: using their identification number (``--extract``),
or using their properties (``--genomic``, ``--maf``, ``--rate`` and/or
``--info``).

.. note::

   It is possible to extract from multiple *impute2* files at the same time (by
   specifying multiple input files).


Extraction by ID
"""""""""""""""""

To extract markers using their identification number, you need a file
containing the list of marker to extract (one marker per line).

.. code-block:: bash

   cd $HOME/genipe_tutorial/extraction

   echo "rs76139713:51137523:C:T" > marker_list.txt
   echo "rs372879164:17037188:A:G" >> marker_list.txt

This ``marker_list.txt`` file will contain the following:

.. code-block:: text

   rs76139713:51137523:C:T
   rs372879164:17037188:A:G

Then, the following command (using the ``--extract`` option) will extract those
two markers from the *impute2* file.

.. code-block:: bash

   impute2-extractor \
       --impute2 ../genipe/chr22/final_impute2/chr22.imputed.impute2.gz \
       --extract marker_list.txt

.. note::

   To gather a list of marker identification numbers, refer to the file
   ``chr22.imputed.map``, which contains the list of all sites in the *impute2*
   file.


Extraction by characteristics
""""""""""""""""""""""""""""""

There are four ways to extract markers according to their characteristics. The
first way is to specify the genomic location of the markers to extract (*i.e.*
the ``--genomic`` option). The second way is to specify a minor allele
frequency threshold (*i.e.* the ``--maf`` option). The third way is to specify
a call rate threshold (*i.e.* the ``--rate`` option). The fourth and final way
is to specify an information value threshold (*i.e.* the ``--info`` option).
Those four ways can be used at the same time (*e.g.* to get markers in a
specific genomic range and a specific call rate).

For example, to extract markers with a MAF :math:`\geq` 0.05 located in the
*CYP2D6* gene, perform the following command:

.. code-block:: bash

   cd $HOME/genipe_tutorial/extraction

   impute2-extractor \
       --impute2 ../genipe/chr22/final_impute2/chr22.imputed.impute2.gz \
       --genomic chr22:42522501-42526883 \
       --maf 0.05 \
       --out cyp2d6_common

To gather all markers with a MAF :math:`\geq` 0.05 and a call rate :math:`\geq`
0.99, perform the following command:

.. code-block:: bash

   impute2-extractor \
       --impute2 ../genipe/chr22/final_impute2/chr22.imputed.impute2.gz \
       --maf 0.05 \
       --rate 0.99 \
       --out common_complete


.. _extract-tut-output-files:

Output files
^^^^^^^^^^^^^

The output files will depend on the output format selected (the ``--format``
option). You can specify either ``impute2``, ``dosage`` and/or ``calls``, for
the *impute2* format (*i.e.* three probabilities per sample), the *dosage*
format (*i.e.* one value between 0 and 2 per sample), and hard calls (*i.e.*
genotypes).

``.impute2`` file
""""""""""""""""""

This file is generated when the ``impute2`` format is used. It has the same
format as the original *impute2* file.

The general structure of the file contains the following columns (which are
space delimited): the chromosome, the name of the marker, its position and its
two alleles. The subsequent columns correspond to the probabilities of each
genotype (hence, there are three columns per sample). The first value
correspond to the probability of being homozygous of the first allele. The
second value correspond to the probability of being heterozygous. Finally, the
third value correspond to the probability of being homozygous of the second
allele.

The following example shows two lines of the ``.impute2`` file.

.. code-block:: text

   22 rs7289830 16058758 C A 0 0 1 0 0 1 0 1 0 ...
   22 rs6423472 16087621 A G 0 1 0 1 0 0 0 1 0 ...


``.dosage`` file
"""""""""""""""""

This file contains the dosage computed from the *impute2* probabilities. The
general structure of the file contains the following columns (which are
tabulation separated): the chromosome, the position on the chromosome, its
name, its minor and major allele and the dosage value. The dosage values vary
between 0 and 2 (inclusively), where values close to 0 represent a higher
chance of been homozygous of the major allele, values close to 1 represent a
higher chance of been heterozygous, and values close to 2 represent a higher
chance of been homozygous of the minor allele.

The following example shows two lines of the ``.dosage`` file.

.. code-block:: text

   22	16058758	rs7289830	C	A	0.0	0.0	1.0	...
   22	16087621	rs6423472	A	G	1.0	2.0	1.0	...

.. note::

   Dosage values computed from probabilities that are below the quality
   threshold (specified by the ``--prob`` option) will have a missing value of
   ``nan``.


``.calls`` file
""""""""""""""""

This file contains the hard calls computed from the *impute2* probabilities. It
has the same format as a transposed pedfile (from *Plink*). The general
structure of the file contains the following columns (which are tabulation
separated): the chromosome, the marker name, the genetic position, the genomic
location, and the hard calls.

The following example shows two lines of the ``.calls`` file.

.. code-block:: text

   22	rs7289830	0	16058758	A A	A A	C A	...
   22	rs6423472	0	16087621	A G	A A	A G	...

.. note::

   Hard calls computed from probabilities that are below the quality threshold
   (specified by the ``--prob`` option) will have a missing value of ``0 0``.


.. _extract-tut-usage:

Usage
^^^^^^

The following command will display the documentation for the extraction
analysis in the console:

.. code-block:: console

   $ impute2-extractor --help
   usage: impute2-extractor [-h] [-v] [--debug] --impute2 FILE [FILE ...]
                            [--out PREFIX] [--format FORMAT [FORMAT ...]]
                            [--prob FLOAT] [--extract FILE]
                            [--genomic CHR:START-END] [--maf FLOAT]
                            [--rate FLOAT] [--info FLOAT]

   Extract imputed markers located in a specific genomic region. This script is
   part of the 'genipe' package, version 1.2.2).

   optional arguments:
     -h, --help            show this help message and exit
     -v, --version         show program's version number and exit
     --debug               set the logging level to debug

   Input Files:
     --impute2 FILE [FILE ...]
                           The output from IMPUTE2.

   Output Options:
     --out PREFIX          The prefix of the output files. [impute2_extractor]
     --format FORMAT [FORMAT ...]
                           The output format. Can specify either 'impute2' for
                           probabilities (same as impute2 format, i.e. 3 values
                           per sample), 'dosage' for dosage values (one value
                           between 0 and 2 by sample), or 'calls' for hard calls.
                           ['impute2']
     --prob FLOAT          The probability threshold used when creating a file in
                           the dosage or call format. [0.9]

   Extraction Options:
     --extract FILE        File containing marker names to extract.
     --genomic CHR:START-END
                           The range to extract (e.g. 22 1000000 1500000). Can be
                           use in combination with '--rate', '--maf' and '--
                           info'.
     --maf FLOAT           Extract markers with a minor allele frequency equal or
                           higher than the specified threshold. Can be use in
                           combination with '--rate', '--info' and '--genomic'.
     --rate FLOAT          Extract markers with a completion rate equal or higher
                           to the specified threshold. Can be use in combination
                           with '--maf', '--info' and '--genomic'.
     --info FLOAT          Extract markers with an information equal or higher to
                           the specified threshold. Can be use in combination
                           with '--maf', '--rate' and '--genomic'.

.. note::

   When using the ``--index`` option, only the indexation (of files without an
   index) will be performed.

