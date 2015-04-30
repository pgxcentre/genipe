Input Files
===========

Multiple input files are required for :py:mod:`gwip` in order to perform the
genome-wide imputation. Those includes the original genotypes of the study
cohort, and the reference files (for imputation). For a description of the
required files, see the :ref:`gwip-tut-input-files` section of the tutorial.

Here is a typical project directory structure. Again, refer to the
:ref:`gwip-tut-input-files` section of the tutorial for more information.

.. code-block:: text

   $HOME/gwip_tutorial/
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
   ├── gwip_config.ini  # OPTIONAL (--use-drmaa, --drmaa-config)
   │
   ├── hg19/
   │   ├── hg19.fasta
   │   └── hg19.fasta.fai
   │
   └── preamble.txt     # OPTIONAL (--use-drmaa, --preamble)

