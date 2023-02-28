Table of Contents
_________________

1. Quick start
.. 1. Install the most recent version of the pipeline
.. 2. Generating config files
.. 3. Index files
.. 4. Before starting the pipeline
.. 5. Starting the pipeline
.. 6. Checking the pipeline status





1 Quick start
=============

1.1 Install the most recent version of the pipeline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  This and all other steps assume that users are using the *scc1* or
  *scc4* clusters on the SCC.

  Run the following command inside the directory that you are going to
  use for the RNAseq analysis:

  ,----
  | git clone https://github.com/mpyatkov/rnaseq_pipeline.git PROJECT_NAME
  `----

  This command creates the PROJECT_NAME directory. PROJECT_NAME can be
  anything, but providing a sensible name will help you know which
  analyses are included in your output folders.

  Inside the *PROJECT_NAME* directory users will find several files and
  a directory named *Scripts*

  - make_cleanup.sh - this script removes large files that are no longer
    needed (BAM files)
  - make_archive.sh - this script compresses the directory with pipeline
    analysis results for future use; it also asks users about the
    directory cleaning
  - README.txt - this file, with pipeline instructions
  - VERSIONS_# - file which contains the pipeline change log
  - Scripts - directory with all scripts required for analysis


1.2 Generating config files
~~~~~~~~~~~~~~~~~~~~~~~~~~~

  Before starting the pipeline, the user needs to configure it.

  *WARNING: DO NOT USE QUOTES WHEN MAKING ENTRIES IN THE CONFIG FILES*
  *WARNING: DO NOT USE SPACES IN VARIABLES IN CONFIG FILES, USE "_"
  *INSTEAD*

  Go to inside the PROJECT_NAME/Scripts/00_Setup_Pipeline/
  directory. There are two files inside this directory
  (01_Pipeline_Setup.py, 02_Run_Pipeline.sh)

  - 01_Pipeline_Setup.py - this script facilitates work with the
    configuration files. The pipeline usually uses it inside the bash
    scripts.
  - 02_Run_Pipeline.sh - the script gives the user several options for
    starting the pipeline.

  First of all, the user needs to run the following command to obtain a
  set of configuration files specific to the user and the user's
  folders:
  ,----
  | ./01_Pipeline_Setup.py # or equivalent ./01_Pipeline_Setup.py --generate
  `----

  The above command produces 4 config files:

  1) *Pipeline_Setup.conf* - contains all information about specific
      configuration details of the pipeline; it includes optimal values
      for many variables, which do not usually need to be changed for
      routine runs of the pipeline. However, the user should check the
      [USER] section, because in 99% of cases the config generator sets
      optimal values for all variables.

  2) *Sample_Labels.txt* - This is the most important setup file in the
     pipeline. This file contains information about sample specific
     options. The main columns are:
     - *Group* - each sample must be associated with only one group
     - *Condition_Name* - this serves as a longer, more descriptive name
       *for the group Group*. All samples in a group should have the
       *same *Condition_Name*
     - *Sample_ID* - unique ID for each sample. In the pipeline, sets of
        samples are combined into projects. Special files, called
        "index" files, are used to organize information about samples
        and projects. Within these files, the samples are uniquely
        associated with a pair of FASTQ files for left and right
        sequence reads. Before running the pipeline, the user needs to
        know which index files have already been created and what set of
        samples they refer to. For external public datasets, the user
        needs to provide their own index files. All information about
        how to create index files is presented in the section below.
     - *Description* - a longer, more complete text description of the
        sample
     - *Color* - the color to be used in the UCSC browser

     IMPORTANT NOTE: None of the values in *Sample_Labels.txt* may begin
     with a number, except for *Color*

  3) *Comparisons.txt* - contains information about which groups are to
      be compared for identification of DEGs (differentially expressed
      genes). The main rule here is:
      *Condition_2(Treatment)/Condition_1(Control)*. Information from
      this step is required to generate directories for Differential
      Expression analysis.

  4) *venn_comparisons.txt* (Optional) - Allows the user to specify
     which groups of samples will be analyzed to generate different Venn
     comparisons. The format of this file pretty straightforward,
     example:

  ,----
  | venn_comparisons
  | 1;2
  | 1;2;3
  | 1;3
  | 2;3;4
  | 1;2;3;4
  `----

  Where "1;2" means compare the DE output and build Venn diagrams for
  comparisons 1 and 2 (numbers correspond to Comparison_Number column
  from the *Comparison.txt* configuration file), "1;2;3" - compare the
  DE output for comparisons 1,2,3 etc. The maximum number of DE
  comparisons that can be compared with each other is four (as in
  1;2;3;4). There is no limit to the number of lines entered in
  *venn_comparisons.txt*

  To start the pipeline, in most cases, users only need to fill in
  *Sample_Labels.txt* and *Comparisons.txt*. Complete these
  configuration files and proceed to the next step.


1.3 Index files
~~~~~~~~~~~~~~~

  To analyze public RNA-SEQ DATASETS, sometimes you will need to create
  a special file that will associate unique SAMPLE_ID with a pair of
  FASTQ files. Such unique SAMPLE_IDs allow us to organize and reuse
  already calculated results as follows: after the first analysis, the
  FASTQ files will be mapped to a genome, resulting in a BAM file. The
  resulting BAM FILE will be compressed to A CRAM file and stored on the
  WaxmanLab server. The next time the pipeline is run to re-analyze
  these same samples, it will skip the step of mapping reads to the
  reference genome and reuse the already computed BAM files by unpacking
  them from the corresponding CRAM files.

  Index files contain information about each SAMPLE_ID and FASTQ file,
  and are usually csv files in a special file structure.  By default,
  the pipeline already uses at least one index file, specified in the
  FASTQ_DEFAULT_INDEX variable inside the Pipeline_Setup.conf.  Use who
  needs to use their own index file should put a copy of that index file
  inside the Scripts/00_Setup_Pipeline/ directory, with name suffixed as
  follows "index.csv" (ex. G186_index.csv). When the pipeline detects an
  index file in the 00_Setup_Pipeline/ directory, it will use that index
  file, overriding the FASTQ_DEFAULT_INDEX. Thus, the pipeline will
  first try to find info about each SAMPLE_ID inside the user defined
  index file, before it goes to FASTQ_DEFAULT_INDEX.

  The format for index files is as follows:

  1. First line contains three columns separated by commas:
     - keyword PRJ separates projects from each other
     - PROJECT_NAME
     - PROJECT_PATH - maximal full path to the directory with fastq
       files
  2. Other lines also contain three columns to represent the SAMPLE_ID,
     RELATIVE_PATH_TO_READS1 and RELATIVE_PATH_TO_READS2
     - SAMPLE_ID - unique identifier which contains info about project
       and short ID after the "_" character. The pipeline generates an
       error if the identifier has more than one "_" character. Example
       of a good name for sample id is 'G186_M1' or something similar.
     - PATH to the file with first reads. Path should be relative to the
       PROJECT_PATH. If all users FASTQ files are located in one
       directory, the paths to the reads will be just the file names. If
       FASTQ files are already separated by pairs and stored in
       different directories, then the user will need to add info about
       relative location to the reads file name (see examples below).

  To separate the projects each project should start with a PRJ keyword:

  Index files can contain as many projects as you wish, the main rule
  here is that all SAMPLE_ID must be unique.

  ,----
  | PRJ, PROJECT_NAME1, PROJECT_PATH1
  | SAMPLE_ID1, FASTQ_WITH_READS1, FASTQ_WITH_READS2
  | SAMPLE_ID2, FASTQ_WITH_READS1, FASTQ_WITH_READS2
  | SAMPLE_ID3, FASTQ_WITH_READS1, FASTQ_WITH_READS2
  | ...
  | 
  | PRJ, PROJECT_NAME2, PROJECT_PATH2
  | SAMPLE_ID4, FASTQ_WITH_READS1, FASTQ_WITH_READS2
  | SAMPLE_ID5, FASTQ_WITH_READS1, FASTQ_WITH_READS2
  | SAMPLE_ID6, FASTQ_WITH_READS1, FASTQ_WITH_READS2
  `----

  In the example below, each pair of FASTQ files is stored in its own
  directory, so we need to add the directory to identify the relative
  path when specifying read names:

  ,----
  | PRJ, G186, /net/waxman-server/mnt/data/volume2/Waxman_Illumina_HiSeq_Raw_Data_02/G186/usftp21.novogene.com/raw_data
  | G186_M1, /G186_M1/G186_M1_CKDL200158629-1a-D701-AK1780_HC2WKBBXX_L4_1.fq.gz, /G186_M1/G186_M1_CKDL200158629-1a-D701-AK1780_HC2WKBBXX_L4_2.fq.gz
  | G186_M2, /G186_M2/G186_M2_CKDL200158629-1a-D701-AK1543_HC2WKBBXX_L4_1.fq.gz, /G186_M2/G186_M2_CKDL200158629-1a-D701-AK1543_HC2WKBBXX_L4_2.fq.gz
  | ...
  | 
  `----

  For SINGLE-END data, just use the same fastq file for READ1 and
  READ2. The pipeline will automatically change the necessary
  parameters for downstream programs.

  ,----
  | PRJ, GNNN, /path/to/single-end/fastq/files/directory/
  | GNNN_M1, GNNN_M1_single_end.fq.gz, GNNN_M1_single_end.fq.gz
  | GNNN_M2, GNNN_M2_single_end.fq.gz, GNNN_M2_single_end.fq.gz
  | ...
  | 
  `----


1.4 Before starting the pipeline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  When running the pipeline on *scc4*, the session will be interrupted
  after 15 minutes of inactivity, which will cause a problem when the
  pipeline needs several hours to finish running. To prevent unexpected
  interruption, the user should use the *screen* or *tmux* - terminal
  multiplexers which will keep session alive. The following set of
  commands is enough to insure a successful pipeline run:

  1. tmux (run tmux session and connect to it)
  2. tmux attach (connect to tmux session if not connected)
  3. inside the tmux session press "Ctrl-B D" to detach from session


1.5 Starting the pipeline
~~~~~~~~~~~~~~~~~~~~~~~~~

  There are several options for starting the pipeline, to check all of
  them, run the following command:

  ,----
  | ./02_Run_Pipeline.sh
  `----

  The output will contain the following lines:

  ,----
  | <option> = START - light version of FULL, without recalculation of output for samples that already exist
  | <option> = FULL - full pipeline run with recalculation of all steps
  | <option> = DEONLY recalculate DE and summary directories (09abcd,12,13,14)
  | <option> = TRACKS creates UCSC tracks (steps 01,02,04 and 11 without full recalculation if that possible)
  | <option> = start_step (example. 05) start to run pipeline from specific step
  `----

  In the most cases START option is the optimal way how to start the
  pipeline, because it allows you to reuse already precalculated BAM
  files obtained from previous pipeline users.  But when the user is
  working with new Fastq files, and needs information about quality
  checks (such as FASTQ, Picard RnaSeqMetrics, etc.), then the FULL
  option would be a Better choice. The FULL option informs the pipeline
  that all intermediate files (BAM, counts, etc) should be recalculated
  from the FASTQ source files, including all quality control checks.

  To start the pipeline users should run the following command:

  ,----
  | tmux 
  | ./02_Run_Pipeline.sh START
  | Press "Ctrl-B D" to detach from current tmux session
  `----


1.6 Checking the pipeline status
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  Next time when users need to check the pipeline status they need to
  attach to the last tmux session using the next command:

  ,----
  | tmux attach
  `----

  In addition, users can check the cluster queue to get information
  about current execution step of the pipeline:

  ,----
  | qstat -u BU_USER_NAME
  `----

  Finally, after the pipeline run is completed, the user can enter the
  Terminal View mode, with the command (right) Control + b, followed by
  the use of Page Up and page Down, to view all of the pipeline notes
  output to the screen.

  To exit the Command View mode, enter Escape.
