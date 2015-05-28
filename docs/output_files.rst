Output Files
=============

The :py:mod:`genipe` pipeline generates output files for each of the autosomal
chromosome (from 1 to 22). To get a description of those file, refer to the
:ref:`genipe-tut-output-files` section of the tutorial.

In summary, here is the structure of the output files. Again, refer to the
:ref:`genipe-tut-output-files` section of the tutorial for more information.

.. code-block:: text

   genipe/
   │
   ├── chr1/
   │   ├── chr1.1_5000000.impute2
   │   ├── chr1.1_5000000.impute2_info
   │   ├── chr1.1_5000000.impute2_info_by_sample
   │   ├── chr1.1_5000000.impute2_summary
   │   ├── chr1.1_5000000.impute2_warnings
   │   ├── ...
   │   ├── chr1.final.bed
   │   ├── chr1.final.bim
   │   ├── chr1.final.fam
   │   ├── chr1.final.log
   │   ├── chr1.final.phased.haps
   │   ├── chr1.final.phased.ind.me
   │   ├── chr1.final.phased.ind.mm
   │   ├── chr1.final.phased.log
   │   ├── chr1.final.phased.sample
   │   ├── chr1.final.phased.snp.me
   │   ├── chr1.final.phased.snp.mm
   │   ├── ...
   │   │
   │   └── final_impute2/
   │       ├── chr1.imputed.alleles
   │       ├── chr1.imputed.completion_rates
   │       ├── chr1.imputed.good_sites
   │       ├── chr1.imputed.impute2.gz
   │       ├── chr1.imputed.imputed_sites
   │       ├── chr1.imputed.log
   │       ├── chr1.imputed.maf
   │       ├── chr1.imputed.map
   │       └── chr1.imputed.sample
   │
   ├── .../
   │
   ├── chromosome_lengths.txt
   ├── frequency_pie.pdf
   ├── genipe.log
   ├── markers_to_exclude.txt
   ├── markers_to_flip.txt
   │
   ├── missing
   │   ├── missing.imiss
   │   ├── missing.lmiss
   │   └── missing.log
   │
   ├── report
   │   ├── frequency_pie.pdf
   │   ├── Makefile
   │   ├── references.bib
   │   ├── references.bst
   │   └── report.tex
   │
   └── tasks.db

