Table of Contents
_________________

1. Quick start
.. 1. Install the last version of pipeline
.. 2. Generate config files
.. 3. Before start the pipeline
.. 4. Start the pipeline




1 Quick start
=============

1.1 Install the last version of pipeline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  This and all other steps assume that you are using *scc4* cluster.

  Run the following command in the directory that you are going to use
  for RNAseq analysis:

  ,----
  | git clone https://github.com/mpyatkov/rnaseq_pipeline.git PROJECT_NAME
  `----

  PROJECT_NAME can be any, but it would be better if the name was
  associated with the project


1.2 Generate config files
~~~~~~~~~~~~~~~~~~~~~~~~~

  inside the *PROJECT_NAME* directory you will find several files and a
  directory *Scripts*

  - make_cleanup.sh - a script that allows users to free a directory
    from large files that are no longer needed
  - make_archive.sh - a script that allows users to compress directory
    with results for future analysis
  - README.txt - this file
  - VERSIONS_# - file which contains the pipeline change log
  - Scripts - directory with all scripts required for analysis

  To start the pipeline first of all you need to configure it.

  Go to the Scripts/00_Setup_Pipeline Inside this directory you will
  find 2 files (01_Pipeline_Setup.py, 02_Run_Pipeline.sh)

  - 01_Pipeline_Setup.py - a script that allows to work with
    configuration files
  - 02_Run_Pipeline.sh - a script for starting the pipeline with
    different options

  First of all, you need to run the following command to get the
  configuration files:

  ,----
  | ./01_Pipeline_Setup.py
  `----

  The command produces 4 config files:

  1) *Pipeline_Setup.conf* - contains all information about locations of
      samples, project, and other specific details. To run the pipeline
      it will be enough to check out only [USER] section in 99% config
      generator already fill set correct values for all variables.

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
        start the pipeline you should know what index files is created
        and what set of samples their contain.
     - *Description* - the full description of the sample
     - *Color* - the color to be used in the UCSC browser

  3) *Comparisons.txt* - contains information about which groups should
      be compared. The main rule here:
      *Condition_2(Treatment)/Condition_1(Control)*. Information from
      this step is required to generate directories for Differential
      Expression analysis.

  4) *venn_comparisons.txt* (Optional) - Allows to specify which groups
      of samples will be subject to different venn comparisons (numbers
      correspond to Comparison_Number column from the Comparison.txt):

  *WARNING: DO NOT USE QUOTES IN CONFIG FILES* WARNING: DO NOT USE
  *SPACES IN VARIABLES IN CONFIG FILES, USE "_" INSTEAD*

  To start the pipeline in most case is required to fill only
  *Sample_Labels.txt* and *Comparisons.txt* Fill out this configuration
  files and go to the next step.


1.3 Before start the pipeline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  In the case of *scc4*, the session will be interrupted after 15
  minutes of inactivity. But for the pipeline it is required to have
  active session for several hours. To prevent unexpected interruption
  user should use the *screen* or *tmux* - terminal multiplexers. The
  following set of commands is enough for successful run the pipeline.

  1. tmux (run tmux session and connect to it)
  2. tmux attach (connect to tmux session if not connected)
  3. inside the tmux session press "Ctrl-B D" to detach from session


1.4 Start the pipeline
~~~~~~~~~~~~~~~~~~~~~~

  To start the pipeline you should run the following command.

  ,----
  | ./02_Run_Pipeline.sh START
  `----

  After that user can detach current tmux session using "Ctrl+B D" and
  check the current status using *qstat* command

  ,----
  | qstat -u BU_USER_NAME
  `----
