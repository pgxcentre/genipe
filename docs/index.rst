.. genipe documentation master file, created by
   sphinx-quickstart on Mon Mar 16 12:58:45 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Genome-Wide Imputation Pipeline
================================


Introduction
-------------

The :py:mod:`genipe` (GENome-wide Imputation PipelinE) module provides an easy
and efficient way of performing genome-wide imputation analysis using the three
commonly used tools `PLINK <http://pngu.mgh.harvard.edu/~purcell/plink/>`_,
`SHAPEIT <https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html>`_ and
`IMPUTE2 <https://mathgen.stats.ox.ac.uk/impute/impute_v2.html>`_.

A database keeps track of all executed steps and enables post-failure relaunch
of the pipeline where it last stopped, saving processing time and resources.

A report is automatically generated at the end of the imputation process to
easily assess important quality metrics. This report can be rendered to PDF
using LaTeX. For information on rendering, refer to :ref:`this section
<genipe-tut-compile-report>` of the main :ref:`pipeline tutorial
<genipe-tut-page>`.

Finally, it also provides a useful standalone tool to perform statistical
analysis on imputed (dosage) data (such as linear, logistic, repeated
measurements, survival analysis, or `SKAT
<http://www.hsph.harvard.edu/skat/>`_). For more information about execution
time, see the :ref:`stats-exec-time` page.

.. toctree::
   :maxdepth: 2

   installation
   tutorials
   execution_time
   parameters
   input_files
   output_files
   implementation


.. _genipe-usage:

Usage
------

.. code-block:: console

   $ genipe-launcher --help
   usage: genipe-launcher [-h] [-v] [--debug] [--thread THREAD] --bfile PREFIX
                          [--reference FILE] [--chrom CHROM [CHROM ...]]
                          [--output-dir DIR] [--bgzip] [--use-drmaa]
                          [--drmaa-config FILE] [--preamble FILE]
                          [--shapeit-bin BINARY] [--shapeit-thread INT]
                          [--shapeit-extra OPTIONS] [--plink-bin BINARY]
                          [--hap-template TEMPLATE] [--legend-template TEMPLATE]
                          [--map-template TEMPLATE] --sample-file FILE
                          [--hap-nonPAR FILE] [--hap-PAR1 FILE] [--hap-PAR2 FILE]
                          [--legend-nonPAR FILE] [--legend-PAR1 FILE]
                          [--legend-PAR2 FILE] [--map-nonPAR FILE]
                          [--map-PAR1 FILE] [--map-PAR2 FILE]
                          [--impute2-bin BINARY] [--segment-length BP]
                          [--filtering-rules RULE [RULE ...]]
                          [--impute2-extra OPTIONS] [--probability FLOAT]
                          [--completion FLOAT] [--info FLOAT]
                          [--report-number NB] [--report-title TITLE]
                          [--report-author AUTHOR]
                          [--report-background BACKGROUND]

   Execute the genome-wide imputation pipeline. This script is part of the
   'genipe' package, version 1.4.2.

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
     --chrom CHROM [CHROM ...]
                           The chromosomes to process. It is possible to write
                           'autosomes' to process all the autosomes (from
                           chromosome 1 to 22, inclusively).
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
     --shapeit-extra OPTIONS
                           SHAPEIT extra parameters. Put extra parameters between
                           single or normal quotes (e.g. --shapeit-extra '--
                           states 100 --window 2').

   Plink Options:
     --plink-bin BINARY    The Plink binary if it's not in the path.

   IMPUTE2 Autosomal Reference:
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

   IMPUTE2 Chromosome X Reference:
     --hap-nonPAR FILE     The IMPUTE2's haplotype file for the non-
                           pseudoautosomal region of chromosome 23.
     --hap-PAR1 FILE       The IMPUTE2's haplotype file for the first
                           pseudoautosomal region of chromosome 23.
     --hap-PAR2 FILE       The IMPUTE2's haplotype file for the second
                           pseudoautosomal region of chromosome 23.
     --legend-nonPAR FILE  The IMPUTE2's legend file for the non-pseudoautosomal
                           region of chromosome 23.
     --legend-PAR1 FILE    The IMPUTE2's legend file for the first
                           pseudoautosomal region of chromosome 23.
     --legend-PAR2 FILE    The IMPUTE2's legend file for the second
                           pseudoautosomal region of chromosome 23.
     --map-nonPAR FILE     The IMPUTE2's map file for the non-pseudoautosomal
                           region of chromosome 23.
     --map-PAR1 FILE       The IMPUTE2's map file for the first pseudoautosomal
                           region of chromosome 23.
     --map-PAR2 FILE       The IMPUTE2's map file for the second pseudoautosomal
                           region of chromosome 23.

   IMPUTE2 Options:
     --impute2-bin BINARY  The IMPUTE2 binary if it's not in the path.
     --segment-length BP   The length of a single segment for imputation. [5e+06]
     --filtering-rules RULE [RULE ...]
                           IMPUTE2 filtering rules (optional).
     --impute2-extra OPTIONS
                           IMPUTE2 extra parameters. Put the extra parameters
                           between single or normal quotes (e.g. --impute2-extra
                           '-buffer 250 -Ne 20000').

   IMPUTE2 Merger Options:
     --probability FLOAT   The probability threshold for no calls. [<0.9]
     --completion FLOAT    The completion rate threshold for site exclusion.
                           [<0.98]
     --info FLOAT          The measure of the observed statistical information
                           associated with the allele frequency estimate
                           threshold for site exclusion. [<0.00]

   Automatic Report Options:
     --report-number NB    The report number. [genipe automatic report]
     --report-title TITLE  The report title. [genipe: Automatic genome-wide
                           imputation]
     --report-author AUTHOR
                           The report author. [Automatically generated by genipe]
     --report-background BACKGROUND
                           The report background section (can either be a string
                           or a file containing the background. [General
                           background]


Citing genipe
--------------

If you use :py:mod:`genipe` in any published work, please cite the paper
describing the tool.

   Lemieux Perreault LP, Legault MA, Asselin G, Dubé MP:
   **genipe: an automated genome-wide imputation pipeline with automatic
   reporting and statistical tools.**
   *Bioinformatics* 2016, **32** (23): 3661-3663
   [DOI:`10.1093/bioinformatics/btw487 <http://dx.doi.org/10.1093/bioinformatics/btw487>`_].


About
------

This project was conducted at the
`Beaulieu-Saucier Pharmacogenomics Centre <http://www.pharmacogenomics.ca>`_ of
the `Montreal Heart Institute <https://www.icm-mhi.org>`_. The aim was to speed
up (and automatize) the imputation process of the whole genome and facilitate
downstream data analysis.

.. image:: _static/logo_pgx.png
   :width: 350px
   :align: center

.. image:: _static/logo_statgen.png
   :width: 200px
   :align: center

.. image:: _static/logo_icm.jpg
   :width: 350px
   :align: center
