.. genipe documentation master file, created by
   sphinx-quickstart on Mon Mar 16 12:58:45 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Genome-Wide Imputation Pipeline
================================


Introduction
-------------

The :py:mod:`genipe` (GENome-wide Imputation PipelinE) module provides an easy
an efficient way of performing genome-wide imputation analysis using the three
commonly used softwares `PLINK <http://pngu.mgh.harvard.edu/~purcell/plink/>`_,
`SHAPEIT <https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html>`_ and
`IMPUTE2 <https://mathgen.stats.ox.ac.uk/impute/impute_v2.html>`_.

A quality metrics report is automatically generated at the end of the
imputation process to easily assess the quality of the analysis. The report is
compiled into a PDF. For information on how to compile the report, refer to the
:ref:`genipe-tut-compile-report` section in the main :ref:`genipe-tut-page`.

Finally, it also provides a useful standalone tool to perform statistical
analysis on imputed (dosage) data (such as linear, logistic or survival
regressions, or `SKAT <http://www.hsph.harvard.edu/skat/>`_ analysis of rare
variants).

.. toctree::
   :maxdepth: 2

   installation
   tutorials
   input_files
   output_files
   implementation


.. _genipe-usage:

Usage
------

.. code-block:: console

   $ genipe-launcher --help
   usage: genipe-launcher [-h] [-v] [--debug] [--thread THREAD] --bfile PREFIX
                          [--reference FILE] [--output-dir DIR] [--bgzip]
                          [--use-drmaa] [--drmaa-config FILE] [--preamble FILE]
                          [--shapeit-bin BINARY] [--shapeit-thread INT]
                          [--plink-bin BINARY] [--impute2-bin BINARY]
                          [--segment-length BP] --hap-template TEMPLATE
                          --legend-template TEMPLATE --map-template TEMPLATE
                          --sample-file FILE [--filtering-rules RULE [RULE ...]]
                          [--probability FLOAT] [--completion FLOAT]
                          [--report-number NB] [--report-title TITLE]
                          [--report-author AUTHOR]

   Execute the genome-wide imputation pipeline. This script is part of the
   'genipe' package, version 1.1.0.

   optional arguments:
     -h, --help            show this help message and exit
     -v, --version         show program's version number and exit
     --debug               set the logging level to debug
     --thread THREAD       number of threads [1]

   Input Options:
     --bfile PREFIX        The prefix of the binary pedfiles (input data).
     --reference FILE      The human reference to perform an initial strand check
                           (useful for genotyped markers not in the IMPUTE2
                           reference files) (optional).

   Output Options:
     --output-dir DIR      The name of the output directory. [genipe]
     --bgzip               Use bgzip to compress the impute2 files.

   HPC Options:
     --use-drmaa           Launch tasks using DRMAA.
     --drmaa-config FILE   The configuration file for tasks (use this option when
                           launching tasks using DRMAA). This file should
                           describe the walltime and the number of
                           nodes/processors to use for each task.
     --preamble FILE       This option should be used when using DRMAA on a HPC
                           to load required module and set environment variables.
                           The content of the file will be added between the
                           'shebang' line and the tool command.

   SHAPEIT Options:
     --shapeit-bin BINARY  The SHAPEIT binary if it's not in the path.
     --shapeit-thread INT  The number of thread for phasing. [1]

   Plink Options:
     --plink-bin BINARY    The Plink binary if it's not in the path.

   IMPUTE2 Options:
     --impute2-bin BINARY  The IMPUTE2 binary if it's not in the path.
     --segment-length BP   The length of a single segment for imputation. [5e+06]
     --hap-template TEMPLATE
                           The template for IMPUTE2's haplotype files (replace
                           the chromosome number by '{chrom}', e.g.
                           '1000GP_Phase3_chr{chrom}.hap.gz').
     --legend-template TEMPLATE
                           The template for IMPUTE2's legend files (replace the
                           chromosome number by '{chrom}', e.g.
                           '1000GP_Phase3_chr{chrom}.legend.gz').
     --map-template TEMPLATE
                           The template for IMPUTE2's map files (replace the
                           chromosome number by '{chrom}', e.g.
                           'genetic_map_chr{chrom}_combined_b37.txt').
     --sample-file FILE    The name of IMPUTE2's sample file.
     --filtering-rules RULE [RULE ...]
                           IMPUTE2 filtering rules (optional).

   IMPUTE2 Merger Options:
     --probability FLOAT   The probability threshold for no calls. [<0.9]
     --completion FLOAT    The completion rate threshold for site exclusion.
                           [<0.98]

   Automatic Report Options:
     --report-number NB    The report number. [genipe automatic report]
     --report-title TITLE  The report title. [genipe: Automatic genome-wide
                           imputation]
     --report-author AUTHOR
                           The report author. [Automatically generated by genipe]


About
------

This project was initiated at the
`Beaulieu-Saucier Pharmacogenomics Centre <http://www.pharmacogenomics.ca>`_ of
the `Montreal Heart Institute <https://www.icm-mhi.org>`_. The aim was to speed
up (and automatize) the imputation process for the whole genome.

.. image:: _static/logo_pgx.png
   :width: 350px
   :align: center

.. image:: _static/logo_statgen.png
   :width: 200px
   :align: center

.. image:: _static/logo_icm.jpg
   :width: 350px
   :align: center
