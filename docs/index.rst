.. gwip documentation master file, created by
   sphinx-quickstart on Mon Mar 16 12:58:45 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Genome-Wide Imputation Pipeline
================================

Introduction
-------------

The :py:mod:`gwip` (Genome-Wide Imputation Pipeline) module provides an easy an
efficient way of performing genome-wide imputation analysis using the three
commonly used softwares `PLINK <http://pngu.mgh.harvard.edu/~purcell/plink/>`_,
`SHAPEIT <https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html>`_ and
`IMPUTE2 <https://mathgen.stats.ox.ac.uk/impute/impute_v2.html>`_. It also
provides a useful standalone tool to perform statistical analysis on imputed
(dosage) data (such as linear, logistic or survival regressions, or
`SKAT <http://www.hsph.harvard.edu/skat/>`_ analysis of rare variants).

.. toctree::
   :maxdepth: 2

   installation
   tutorial
   input_files
   output_files

