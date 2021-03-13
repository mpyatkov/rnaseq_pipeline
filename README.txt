Table of Contents
_________________

1. Quick start
.. 1. Install the last version of pipeline
.. 2. Generate config files
.. 3. Index files
.. 4. Before start the pipeline
.. 5. Start the pipeline





1 Quick start
=============

1.1 Install the last version of pipeline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  This and all other steps assume that users are using *scc1* or *scc4*
  clusters

  Run the following command in the directory that users are going to use
  for RNAseq analysis:

  ,----
  | git clone https://github.com/mpyatkov/rnaseq_pipeline.git PROJECT_NAME
  `----

  PROJECT_NAME can be any, but it would be better if the name was
  associated with the project


1.2 Generate config files
~~~~~~~~~~~~~~~~~~~~~~~~~

  inside the *PROJECT_NAME* directory users will find several files and
  a directory *Scripts*

  - make_cleanup.sh - a script that allows users to free a directory
    from large files that are no longer needed
  - make_archive.sh - a script that allows users to compress directory
    with results for future analysis
  - README.txt - this file
  - VERSIONS_# - file which contains the pipeline change log
  - Scripts - directory with all scripts required for analysis

  To start the pipeline first of all users need to configure it.

  Go to the Scripts/00_Setup_Pipeline Inside this directory users will
  find 2 files (01_Pipeline_Setup.py, 02_Run_Pipeline.sh)

  - 01_Pipeline_Setup.py - a script that allows to work with
    configuration files
  - 02_Run_Pipeline.sh - a script for starting the pipeline with
    different options

  First of all, users need to run the following command to get the
  configuration files:

  ,----
  | ./01_Pipeline_Setup.py
  `----

  The command produces 4 config files:

  1) *Pipeline_Setup.conf* - contains all information about locations of
      samples, project, and other specific details. To run the pipeline
      it will be enough to check out only [USER] section, because in 99%
      config generator already set correct values for all variables.

  2) *Sample_Labels.txt* - The most important file in the pipeline. This
      file contains information about how to configure samples. The main
      columns are:
     - *Group* - each sample should be associated with only one group
     - *Condition_Name* - should be considered as a long name for the
        *Group*. Samples from one group should have same
        *Condition_Name*
     - *Sample_ID* - ID for sample. The new version of the pipeline uses
        index files. Each index file contains a set of unique sample_ids
        that are associated with the corresponding fastq files. Before
        start the pipeline users should know which index files are
        created and what set of samples their contain. There are
        situtations when users need to create their own subset of
        samples, for that just read section about index files below.
     - *Description* - the full description of the sample
     - *Color* - the color to be used in the UCSC browser

  3) *Comparisons.txt* - contains information about which groups should
      be compared. The main rule here:
      *Condition_2(Treatment)/Condition_1(Control)*. Information from
      this step is required to generate directories for Differential
      Expression analysis.

  4) *venn_comparisons.txt* (Optional) - Allows to specify which groups
      of samples will be subject to different venn comparisons. By
      default this file contains only header, but the format of this
      file pretty straightforward, example:

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
  from the Comparison.txt configuration file), "1;2;3" - compare the DE
  output for comparisons 1,2,3 etc. The maximum number of DE comparisons
  that can be compared with each other is 4.

  *WARNING: DO NOT USE QUOTES IN CONFIG FILES* WARNING: DO NOT USE
  *SPACES IN VARIABLES IN CONFIG FILES, USE "_" INSTEAD*

  To start the pipeline, in most cases, users only need to fill in
  *Sample_Labels.txt* and *Comparisons.txt*. Complete these
  configuration files and proceed to the next step.


1.3 Index files
~~~~~~~~~~~~~~~

  To analyze the public data, sometimes you will need to create a
  special file that will associate unique SAMPLE_ID with a pair of FASTQ
  files. Such unique SAMPLE_IDs allow us to organize and reuse already
  calculated results as follows: after the first analysis, the FASTQ
  files will be mapped to a genome, resulting in a BAM file. The
  resulting BAM will be compressed to CRAM file and stored on the
  server. Next time the pipeline can skip the step of mapping reads to
  the reference genome and reuse the already computed BAM file by
  unpacking it from the CRAM file.

  The files contained information about SAMPLE_ID and FASTQ files aka
  'index files' are usual csv files which have special structure.  By
  default, the pipeline already uses at least one index file specified
  in the FASTQ_DEFAULT_INDEX variable inside the Pipeline_Setup.conf.
  To use their own index the users should put index file inside the
  Scripts/00_Setup_Pipeline/ directory with name suffixed as follows
  "index.csv" (ex. G186_index.csv), in this case the pipeline
  automatically detect index file and override FASTQ_DEFAULT_INDEX. The
  overriding means that the pipeline will try to find info about
  SAMPLE_ID inside the user defined index first.

  The format for index files is pretty simple:

  1. First line contains three columns separated by commas:
     - keyword PRJ separate projects to each other
     - PROJECT_NAME
     - PROJECT_PATH - maximal full path to the directory with fastq
       files
  2. Other lines contains also three columns to represent the SAMPLE_ID,
     RELATIVE_PATH_TO_READS1 and RELATIVE_PATH_TO_READS2
     - SAMPLE_ID - unique identifier which contains info about project
       and short ID after the "_" character. The pipeline raises error
       if identifier has more than one "_" character. Example of good
       name for sample id is 'G186_M1' or something similar.
     - PATH to the file with first reads. Path should be relative to the
       PROJECT_PATH. If all users FASTQ files located in one directory
       the paths to the reads will be just file names. If FASTQ files
       already separated by pairs and stored in different directories,
       users need to add info about relative location to the reads file
       name (see examples below).
  To separate the projects each project should start with a PRJ keyword:

  Index files can keep as many projects as you wish, the main rule here
  all SAMPLE_ID should be unique.

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

  In this example each pair of FASTQ files is stored in its own
  directory, so we add the each directory to relative path when
  specifying read names:

  ,----
  | PRJ, G186, /net/waxman-server/mnt/data/volume2/Waxman_Illumina_HiSeq_Raw_Data_02/G186/usftp21.novogene.com/raw_data
  | G186_M1, /G186_M1/G186_M1_CKDL200158629-1a-D701-AK1780_HC2WKBBXX_L4_1.fq.gz, /G186_M1/G186_M1_CKDL200158629-1a-D701-AK1780_HC2WKBBXX_L4_2.fq.gz
  | G186_M2, /G186_M2/G186_M2_CKDL200158629-1a-D701-AK1543_HC2WKBBXX_L4_1.fq.gz, /G186_M2/G186_M2_CKDL200158629-1a-D701-AK1543_HC2WKBBXX_L4_2.fq.gz
  | ...
  | 
  `----


1.4 Before start the pipeline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  In the case of *scc4*, the session will be interrupted after 15
  minutes of inactivity. But for the pipeline it is required to have
  active session for several hours. To prevent unexpected interruption
  user should use the *screen* or *tmux* - terminal multiplexers which
  will keep session alive. The following set of commands is enough for
  the successful pipeline run.

  1. tmux (run tmux session and connect to it)
  2. tmux attach (connect to tmux session if not connected)
  3. inside the tmux session press "Ctrl-B D" to detach from session


1.5 Start the pipeline
~~~~~~~~~~~~~~~~~~~~~~

  There are several options for starting the pipeline, to check all of
  them, run the following command:

  ,----
  | ./02_Run_Pipeline.sh
  `----

  The output will contain the following lines:

  ,----
  | <option> = START - light version of FULL without recalculation already existed samples
  | <option> = FULL - full pipeline run with recalculation of all steps
  | <option> = DEONLY recalculate DE and summary directories (09abcd,12,13,14)
  | <option> = TRACKS creates UCSC tracks (01,02,04 and 11 without full recalculation if that possible)
  | <option> = start_step (example. 05) start from specific step
  `----

  In the most cases START option is optimal way how to start the
  pipeline, because it allows to reuse already precalculated BAM files
  obtained from previous pipeline users.  But if users are going to work
  with fairly new data and need information about quality checks such as
  (FASTQ, Picard RnaSeqMetrics, etc.), then the FULL option would be a
  possible choice. The FULL option informs the pipeline that all
  intermediate files (BAM, counts, etc) should be recalculated from
  FASTQ sources with all quality control checks.

  To start the pipeline users should run the following command:

  ,----
  | tmux 
  | ./02_Run_Pipeline.sh START
  | Press "Ctrl-B D" to detach from current tmux session
  `----

  Next time when users need to check the pipeline status they need to
  attach to the last tmux session using the next command:

  ,----
  | tmux attach
  `----

  In addition users can check the cluster queue to get information about
  current execution step of the pipeline:

  ,----
  | qstat -u BU_USER_NAME
  `----
