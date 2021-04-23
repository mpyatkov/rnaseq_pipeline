Table of Contents
_________________

1. Why?
2. Software
3. Variables in config file
4. Step description
5. Tracks creating
6. How to compare two arbitrary groups using taco_refcomp?





1 Why?
======

  TACO produces meta-assembly from the GTF files obtained as input. The
  resulted Meta-assembly (GTF and bed files) includes all possible
  isophorms with abundance and FPKM value for each isophorm
  separatelly. For further analysis the GTF as well as the bigBed files
  will be stored in the Job_Summary directory. In addition the pipeline
  compares TACO meta-assembly and reference GTF (which can be set by
  TACO_REFCOMP_REFERENCE variable in the Pipeline_Setup.conf) using
  *taco_refcomp* program.


2 Software
==========

  The pipeline will be utilize *StringTie* assembler to convert each
  SAMPLE_ID BAM file to the GTF file. If, for some reasons, you would
  like to use something different (ex. Cufflinks), please pay attention
  to the *run_assembler.qsub* script, it contains all information what
  you need.

  *TACO* is used for gathering the individual transcriptomes. Please
   look inside the *run_taco.qsub* script to change the default TACO's
   parameters if it is required.

  *TACO_REFCOMP* is used to compare reference and test GTF files.


3 Variables in config file
==========================

  Entry point for this step is *Run_Jobs.sh* script, it exports all
  parameters from Pipeline_Setup.conf. TACO step consist of 3 substeps
  where the outputs of one are inputs to other. The main parameter is
  *TACO_ENABLE* which allows users to activate/deactivate TACO step. The
  other options allow to make small adjestments:
  + *TACO_STRINGTIE_REFERENCE* - path to the reference GTF file for
     StringTie program, which extract all isoforms from the individual
     BAM and output them as the GTF file.
  + *TACO_REFCOMP_REFERENCE* - path to the reference GTF file which will
     be used as reference for the *taco_refcomp* software.
  + *TACO_KEEP_GTF* - 0/1 option which allows users to specify if they
     would like to keep GTF file after taco processing
  + *TACO_REFCOMP_KEEP_GTF* - 0/1 option which allows users to specify
     if they would like to keep GTF file after taco_refcomp processing.

  The last two options set by default to 0 (off), because the output GTF
  files usually HUGE in size.


4 Step description
==================

  The step executes 2 times and produces 2 Job_Summary directories
  (Job_Summary_RefOn, Job_Summary_RefOff). The only difference between
  two executions that in case of RefOn the step will use
  *TACO_STRINGTIE_REFERENCE* for StringTie, in the case of RefOff,
  StringTie will not use any reference GTF and creates de-novo
  assembly. The main reason for this is that StringTie without reference
  (RefOff) can find potentially new isoforms. Thus the following
  description will be the same for both modes (RefOn, RefOff), but keep
  in mind that we run it twice.

  The step consist of two 2 jobs which are dependent to each
  other. First of all the pipeline calculates individual GTF files
  (run_assembler.qsub). Until all individual GTF files are counted we
  cannot run the TACO meta-assembler.

  On the second step TACO meta-assemble combine all individual GTF files
  by groups which specified in the Sample_Labels.txt. In addition to
  exsisting groups the pipeline creates group which contains all samples
  together, the name for such group is *ALL*. The resulted number of
  groups will be N + 1, where N is the groups specified in the
  Samples_Labels.txt + 1 group with all samples together. The output for
  the TACO meta-assembly step is GTF and bigBed files.

  At the end of the step *taco_refcomp* produces pairwise comparisons
  between TACO_REFCOMP_REFERENCE and meta-assemblies for each group
  which were obtained previosly.

  TACO and job uses the "qsub -hold_jid" command to put own execution
  off until the previous jobs are completed.

  Output for group ALL and groups specified in Samples_Labels.txt is a
  little different, because in first case TACO can produce long names
  for output files and I use different ways to fix this problem. When we
  need to combine samples in one group we use names which look like as
  follows - G186_M1M2M3 where G186 project name and M1M2M3 samples
  inside the group. But for group ALL scripts generates really long
  names, because we combine all samples together. In addition, I create
  a file with description that includes information about who created
  the file, when, and what samples were used to create it.

  For the groups which specified in the Sample_Labels.txt the pipeline
  does absolutely the same job as in previous paragraph, but instead of
  HASH it uses original combined names for output (ex. G186_M1M2M3M4),
  because usually one group contains no more the 4-5 samples.

  The TACO step does not copy anything from the server. It copies bigBed
  file and description file to server if these files do not exist.


5 Tracks creating
=================

  In addition TACO step creates tracks
  (${DATASET_LABEL}-TACO_Tracks.txt) for each bigBed file and copy it to
  Job_Summary as well as
  ${VM_DIR_UCSC}/PERSONAL/${BU_USER}/${DATASET_LABEL}/UCSC_Track_Lines/.


6 How to compare two arbitrary groups using taco_refcomp?
=========================================================

  First of all, user should get the GTF files for groups. For that, user
  should specify option TACO_KEEP_GTF to 1. After that, user need to
  recalculate TACO to obtain GTF files in the same directory where
  bigBed files are located. On the next step it is required to create
  cluster job using the following template:

  ,----
  | #!/bin/bash -l
  | 
  | #$ -P PROJECT # need to specify correct project
  | #$ -o log.txt
  | #$ -N taco_refcomp
  | #$ -pe omp 8
  | #$ -l mem_per_core=2G
  | #$ -l h_rt=01:00:00
  | 
  | module load miniconda
  | conda activate --stack /projectnb/wax-es/routines/condaenv/isoforms
  | 
  | set -o errexit
  | set -o pipefail
  | set -o nounset
  | 
  | REFERENCE_GTF=/path/to/gtf/which/will/be/reference
  | TESTED_GTF=/path/to/gtf/which/will/be/test
  | TACO_REFCOMP_OUTPUT=taco_refcomp
  | taco_refcomp -p ${NSLOTS} -o taco_refcomp -r ${REFERENCE_GTF}  -t ${TESTED_GTF}
  `----

  and run it using qsub. After ~10-15 minutes taco_refcomp script will
  create directory containing GTF and TSV files with necessary
  information.

  The other source of information is an official site:
  <https://tacorna.github.io/> which contains additional options for
  TACO and TACO_REFCOMP programs.
