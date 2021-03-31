Table of Contents
_________________

1. Common information
.. 1. Why?
.. 2. Software
2. TACO MODE
.. 1. ALL SAMPLES TOGETHER (TACO_MODE=1)
.. 2. SAMPLES BY GROUPS (TACO_MODE=2)
3. Step execution
4. Tracks creating





1 Common information
====================

1.1 Why?
~~~~~~~~

  TACO produces meta-assembly from GTF files obtained as input. The
  resulted Meta-assembly (GTF and bed files) includes all possible
  isophorms with abundance and FPKM value for each isophorm
  separatelly. For further analysis the GTF as well as the bigBed files
  will be stored in the Job_Summary directory.


1.2 Software
~~~~~~~~~~~~

  To convert each SAMPLE_ID BAM file to GTF the pipeline will be utilize
  *StringTie* assembler. If for some reasons you would like to use
  something different (ex. Cufflinks), please pay attention to
  *run_assembler.qsub* file, it contains all information what you need.

  *TACO* is used for gathering individual transcriptomes. Please look
   inside the *run_taco.qsub* script to change the default TACO's
   parameters if it is required.


2 TACO MODE
===========

  Entry point for this step is *Run_Jobs.sh* script, it exports all
  parameters from Pipeline_Setup.conf. The main parameter is *TACO_MODE*
  which allows to change assembling behavior of whole step. There are
  three options:

  - TACO_MODE=0 # skip this step without any analysis
  - TACO_MODE=1 # combine all samples together and consider them as one
    group.
  - TACO_MODE=2 # consider different groups of samples separately as
    they are described in Sample_Label.txt in column (Group) and create
    separate tracks for them.

  Output for option ALL (TACO_MODE=1) and BY_GROUPS (TACO_MODE=2) a
  little different, because in first case TACO_MODE option can produce
  long names for output files and I use different ways to fix this
  problem.


2.1 ALL SAMPLES TOGETHER (TACO_MODE=1)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  When we need to combine samples in one group we use names which look
  like as follows - G186_M1M2M3 where G186 project name and M1M2M3
  samples inside the group. But TACO_MODE=1 generates really long names,
  because we combine all samples together. To prevent such behavior
  instead of long name I use 7 first characters of md5 *HASH* of this
  long name. In addition I create description file which shows who
  created this file, when and real long name.

  Using HASH the pipeline try to detect *BIGBED* on the server
  (${VM_DIR_UCSC}/INDEXED_PROJECTS/project/...). If the *BIGBED* exist
  then the pipeline copies all required files from the server (HASH.gtf,
  HASH.bb, HASH_description.txt) to the *Job_Summary/HASH* directory. If
  not the pipeline recalculates all required files and put them in the
  *Job_Summary/HASH* directory.


2.2 SAMPLES BY GROUPS (TACO_MODE=2)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  In this case the pipeline does absolutely the same job as in previous
  paragraph, but instead of HASH it uses original combined names for
  output (ex. G186_M1M2M3M4), because usually one group contains no more
  the 4-5 samples. The description file will not be created for the same
  reasons.


3 Step execution
================

  The step consist of two jobs which are dependent to each other. First
  of all the pipeline calculates individual GTF files
  (run_assembler.qsub). Until all individual GTF files are counted we
  cannot run the TACO meta-assembler. To postpone the TACO job until all
  jobs are completed the following option of the SGE will be used (qsub
  -hold_jid).


4 Tracks creating
=================

  Tracks for this step will be created before the cluster tasks are
  completed. The final path for the tracks will be
  ${VM_DIR_UCSC}/PERSONAL/${BU_USER}/${DATASET_LABEL}/TACO_Track_Lines.
