Input Files
============

Multiple input files are required for :py:mod:`genipe` including the initial
genotypes and the reference files (for imputation). For a thorough description
of these files, see the :ref:`genipe-tut-input-files` section of the tutorial.

Here is a typical project directory structure. Again, refer to the
:ref:`genipe-tut-input-files` section of the tutorial for more information.

.. code-block:: text

   $HOME/genipe_tutorial/
   │
   ├── 1000GP_Phase3/
   │   ├── 1000GP_Phase3_chr1.hap.gz
   │   ├── 1000GP_Phase3_chr2.hap.gz
   │   ├── ...
   │   ├── 1000GP_Phase3_chr1.legend.gz
   │   ├── 1000GP_Phase3_chr2.legend.gz
   │   ├── ...
   │   ├── 1000GP_Phase3.sample
   │   ├── genetic_map_chr1_combined_b37.txt
   │   ├── genetic_map_chr2_combined_b37.txt
   │   └── ...
   │
   ├── bin/
   │   ├── impute2
   │   ├── plink
   │   └── shapeit
   │
   ├── data/
   │   ├── hapmap_CEU_r23a_hg19.bed
   │   ├── hapmap_CEU_r23a_hg19.bim
   │   └── hapmap_CEU_r23a_hg19.fam
   │
   ├── genipe_config.ini  # OPTIONAL (--use-drmaa, --drmaa-config)
   │
   ├── hg19/
   │   ├── hg19.fasta
   │   └── hg19.fasta.fai
   │
   └── preamble.txt     # OPTIONAL (--use-drmaa, --preamble)

