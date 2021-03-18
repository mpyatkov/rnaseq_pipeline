#!/usr/bin/env python3

# Input:
# 3 config files listed below
# Output:
# option1: export - export env. variables from PIPELINE_CONFIG
# option2: samples - return bash array of triples (DIR,ID,DESCR)
#          from SAMPLES_CONFIG
# option3: generate Diff.Ex directories from SAMPLES and COMPARISONS

from collections import namedtuple
import os
import configparser
from distutils.dir_util import copy_tree
import fileinput
import argparse
import glob
import shutil
import sys
import string

# hide tracebacks
sys.tracebacklimit = 0

PIPELINE_CONFIG = "Pipeline_Setup.conf"
SAMPLES_CONFIG = "Sample_Labels.txt"
COMPARISON_CONFIG = "Comparisons.txt"
VENN_CONFIG = "venn_comparisons.txt"

DEFAULT_VENN_CONFIG = """
venn_comparisons
# Example1 of format 1;2
# Example2 of format 1;2;4;5
""".strip()

DEFAULT_COMPARISON_CONFIG = """
Comparison_Number; Condition_1; Condition_2
1;                 A;           B
2;                 A;           C
""".strip()

DEFAULT_SAMPLES_CONFIG = """
Group;  Condition_Name; Sample_ID;  Description;                 Color
A;      Male;           G186_M1;    Male_Liver_4wk_Control;      255,0,0
A;      Male;           G186_M2;    Male_Liver_4wk_Control;      255,0,0
A;      Male;           G186_M3;    Male_Liver_6wk_Control;      255,0,0
A;      Male;           G186_M4;    Male_Liver_6wk_Control;      255,0,0
B;      Male_STAT5ca;   G186_M5;    Male_Liver_4wk_STAT5ca_AAV;  0,255,0
B;      Male_STAT5ca;   G186_M6;    Male_Liver_4wk_STAT5ca_AAV;  0,255,0
B;      Male_STAT5ca;   G186_M7;    Male_Liver_6wk_STAT5ca_AAV;  0,255,0
B;      Male_STAT5ca;   G186_M8;    Male_Liver_6wk_STAT5ca_AAV;  0,255,0
B;      Male_STAT5ca;   G186_M9;    Male_Liver_6wk_STAT5ca_AAV;  0,255,0
C;      Male_TCPO_2wk;  G186_M10;   Male_Liver_2wk_TCPO;         0,0,255
C;      Male_TCPO_2wk;  G186_M11;   Male_Liver_2wk_TCPO;         0,0,255
C;      Male_TCPO_2wk;  G186_M12;   Male_Liver_2wk_TCPO;         0,0,255
C;      Male_TCPO_2wk;  G186_M13;   Male_Liver_2wk_TCPO;         0,0,255
""".strip()

DEFAULT_PIPELINE_CONFIG = """
[USER]
# Setup BU login user name
BU_USER=CHANGE_USER_NAME

# Setup project name. Choices: (wax-dk,waxmanlab,wax-es)
PROJECT=wax-es

# Setup dataset label
DATASET_LABEL=DEFAULT_LABEL

# The location for temporary files (BAM, counts, fastq) In 99% of
# cases you do not need to change it, the default is already the
# optimal value
DATASET_DIR=PATH_TO_PROJECT

[STEPS]
# Bioanalyzer length (from Bioanalyzer tracings)
BIOANALYZER_LEN=330
# Usually the adaptor length is 60bp (60bp on each end: 120bp total)
ADAPTOR_LEN=60
# Input the adaptor length (bp) for one end
READ_LEN=150

# Depending on the library type, we can take into account information
# about specific strands. At the moment, the laboratory uses the dUTP
# preparation method, which corresponds to "fr-firststrand". Previous
# researchs could use other types of libraries, thus if you are unsure
# about library type just select STRANDEDNESS = 3. It will perform
# step 01 _..., which uses RSEQC (infer_experiment.py) and HISAT to
# detect a read strandedness from the small subset of reads for each
# sample.
# possible values:
# 0 - "unstranded"
# 1 - "firststrand"
# 2 - "secondstrand"
# 3 - "auto" (default)
STRANDEDNESS=3

# HTSeq only feature which allows to count different intersections for
# reads it will be used only for step 09a
# https://htseq.readthedocs.io/en/master/count.html
# 0="union"
# 1="intersection-strict"
# 2="intersection-nonempty"
MODE=2

# if you need tracks for samples and BigWig files use BIGWIG_ENABLE
# option
# 1 - enable bigwig files and tracks
# 0 - disable bigwig files and tracks
BIGWIG_ENABLE=1

# This configuration file contains all specific settings for GTF files
# used in the pipeline, like: name, enable/disable double counting,
# etc. For all options please visit GTF_FILES_DIR directory (see
# below) and check 9abcd.csv default config. The GTF_FILES_CONFIG
# parameter below combines GTF_FILES_DIR and DEFAULT_GTF_CONFIG
# Options: 9d.csv, 9abcd.csv
DEFAULT_GTF_CONFIG=9d.csv

# Default aligner, options:
# 0 - TopHat
# 1 - STAR
DEFAULT_ALIGNER=1

# DEFAULT_FC and DEFAULT_FDR options below impact only on
# PCA/correllation plots (Step 13)
DEFAULT_FC=2
DEFAULT_FDR=0.05

[SYSTEM]
# The location of the GTF files (you will need to have access to
# wax-es to run pipeline)
GTF_FILES_DIR=/projectnb/wax-es/routines/GTF_Files_default

# GTF configuration file which contains all information and options
# for GTF files used in this pipeline. If you would like to use custom
# configuration file, just copy the default file in the current
# directory and assign the path to the custom file. Example:
# GTF_FILES_CONFIG=custom.csv
GTF_FILES_CONFIG=${SYSTEM:GTF_FILES_DIR}/${STEPS:DEFAULT_GTF_CONFIG}

# The location of the Bowtie2 indexes, at present time, this only
# includes indexes for mouse mm9 assembly
BOWTIE2INDEX_DIR=/projectnb/wax-es/routines/BowtieIndex

# Directory with full genomes, at present time, this only
# includes indexes for mouse mm9 assembly
FASTA_DIR=/projectnb2/wax-es/routines/FASTA

# The location of the conda packages. This directory also contains
# scripts and configs and qsub file for setting up environments
CONDA_DIR=/projectnb/wax-es/routines/condaenv

# The location of the HISAT indexes, used to detect STRANDEDNESS in
# STEP 02, at present time, this only includes indexes for mouse mm9
# assembly
HISAT2INDEX_DIR=/projectnb/wax-es/routines/hisat2index

# The location of the STAR indexes, at present time, this only
# includes indexes for mouse mm9 assembly (ExonCollapsed_76k GTF file)
STARINDEX_DIR=/projectnb/wax-es/routines/starindex_EC76K

# The location of the global FASTQ file index for all laboratory NGS
# experiments. You can extract information about sample using the
# following command: ./01_Pipeline_Setup.py --get_sample_info
# SAMPLE_ID Each SAMPLE_ID should be unique in this file.  User should
# check index.csv to make sure that the required files are
# listed. User may add their own index file as specified in
# README.txt, to override use of FASTQ_DEFAULT_INDEX
FASTQ_DEFAULT_INDEX=/projectnb/wax-es/routines/index.csv

# special directory which will contain BIGWIG files. 
# $VM_DIR_UCSC/INDEXED_PROJECTS subdirectory for cram/bigwig/combined_bigwig
# $VM_DIR_UCSC/PERSONAL/USERNAME/DATASET_LABEL subdirectory for track
# lines
VM_DIR_UCSC=/net/waxman-server/mnt/data/waxmanlabvm_home/TRACKS/

# default time limit
TIME_LIMIT=36:00:00
""".strip()


class SampleConfig:
    separator = ";"
   
    def __init__(self, configpath):
        generate_example(configpath, DEFAULT_SAMPLES_CONFIG)
        self.samples = self.__read_config(configpath, self.separator)
       
    def __read_config(self, configpath, separator):

        with open(configpath, "r") as config:
            header_and_data = map(lambda x: x.strip(), config.readlines())
            header_and_data = list(filter(lambda x: len(x.strip()) != 0 and not x.strip().startswith("#"), header_and_data))
            # print(header_and_data)
            header = header_and_data[0].replace(separator, ",")
            Sample = namedtuple('Sample', header)
            result = []
            for sample in header_and_data[1:]:
                sample = self.check_number_of_columns(sample.split(separator))
                sample = list(map(lambda s: self.check(s.strip()), sample))
                result.append(Sample(*sample))
            return result


    # return value or raise error
    def check(self,w):
        if w[0] in string.digits and ',' not in w:
            raise ValueError(f"ERROR: all values in Sample_Labels.txt must start with a character, not a number: {w}")

        pattern=string.ascii_letters+string.digits+"_,"
        for c in w:
            if c not in pattern:
                raise ValueError(f"ERROR: Sample_Labels.txt value contains unexpected character '{c}' in word '{w}'. Only '_' can be use for separation")
        
        return w

    def check_number_of_columns(self,l):
        if len(list(l)) != 5:
            raise ValueError(f"ERROR: incorrect number of columns in line: {';'.join(l)}")
        return l
        
    def samplesByGroup(self,group):
        result = []
        for sample in self.samples:
            if sample.Group == group:
                result.append(sample)
        return result

    def groups(self):
        return sorted(set(map(lambda s: s.Group, self.samples)))

    # def groupsToBash(self):
    #     return " ".join(self.groups())
    
    def samplesToBash(self, group=None, color=False):
        # Sample_ID, Description
        result = []
        if group:
            SAMPLES = self.samplesByGroup(group)
        else:
            SAMPLES=self.samples
            
        for sample in SAMPLES:
            # result.append(sample.Sample_DIR)
            result.append(sample.Sample_ID)
            if group:
                result.append(sample.Condition_Name)
            else:
                result.append(sample.Description)
            if color:
                result.append(sample.Color)
        return " ".join(result)

    def samplesWithColorToBash(self):
        #only Sample_ID, Description, Color
        result = []
        for sample in self.samples:
            # result.append(sample.Sample_DIR)
            result.append(sample.Sample_ID)
            result.append(sample.Description)
            result.append(sample.Color)
        return " ".join(result)

    def __str__(self):
        return "\n".join(map(lambda x: str(x), self.samples))


class ComparisonsConfig:
    separator = ";"

    def __init__(self, configpath):
        generate_example(configpath, DEFAULT_COMPARISON_CONFIG)
        self.comparisons = self.__read_config(configpath, self.separator)

    def __read_config(self, configpath, separator):
        with open(configpath, "r") as config:
            header_and_data = map(lambda x: x.strip(), config.readlines())

            # ignore empty lines
            header_and_data = list(filter(lambda s: len(s.strip()) > 0 and not s.strip().startswith("#"), header_and_data))

            header = header_and_data[0].replace(separator, ",")
            Comparison = namedtuple('Comparison', header)
            result = []
            for comparison in header_and_data[1:]:
                comparison = comparison.split(separator)

                # TODO: too hart 
                # if len(comparison) != 3:
                #     continue

                # strip trailing spaces
                comparison = list(map(lambda s: s.strip(), comparison))
                result.append(Comparison(*comparison))

            return result

    def __str__(self):
        return "\n".join(map(lambda x: str(x), self.comparisons))

    def groups(self):
        groups = set()
        for comparison in self.comparisons:
            groups = groups.union(set(comparison.Condition_1.split(",")))
            groups = groups.union(set(comparison.Condition_2.split(",")))
        return groups


class VennConfig:
    # TODO: check if comparison not presented in comparisons.txt 
    # TODO: check if format is not correct
    separator = ";"

    def __init__(self, configpath):
        generate_example(configpath, DEFAULT_VENN_CONFIG)
        self.venn = self.__read_config(configpath, self.separator)

    def __read_config(self, configpath, separator):
        with open(configpath, "r") as config:
            header_and_data = map(lambda x: x.strip(), config.readlines())
            # ignore empty lines
            header_and_data = list(filter(lambda s: len(s.strip()) > 0 and not s.strip().startswith("#"), header_and_data))
            header = header_and_data[0].replace(separator, ",")
            VennLine = namedtuple('VennLine', header)
            result = []
            for vennline in header_and_data[1:]:

                vennline = vennline.split(separator)

                # strip trailing spaces
                vennline = list(map(lambda s: s.strip(), vennline))
                result.append(VennLine(vennline))

            return result

    def __str__(self):
        return "\n".join(map(lambda x: str(x), self.venn))

    def get_number_of_venns(self):
        return len(self.venn)

    # return bash array by ix
    def get_venn_by_index(self, ix):
        return " ".join(self.venn[ix][0])

    
class GTFconfig():
    separator = ";"

    def __init__(self, config_path, current_path):
        # if configpath is not created use current_path+configpath
        path = config_path
        if os.path.exists(config_path):
            path = config_path
        elif os.path.basename(config_path) == config_path and os.path.exists(current_path+"/"+config_path):
            path = current_path+"/"+config_path
        else:
            print(f"Can not find GTF config file: {config_path}")
            exit(1)

        self.gtf_files = self.__read_config(path, self.separator)

    def __read_config(self, configpath, separator):
        with open(configpath, "r") as config:
            header_and_data = map(lambda x: x.strip(), config.readlines())

            # ignore empty lines
            header_and_data = list(filter(lambda s: len(s.strip()) > 0 and not s.strip().startswith("#"), header_and_data))
            
            header = header_and_data[0].replace(separator,",")
            GTFLine = namedtuple('GTFline', header)
            result = []
            for gtf_line in header_and_data[1:]:
                gtf_line = gtf_line.split(separator)

                # strip trailing spaces
                gtf_line = list(map(lambda s: s.strip(), gtf_line))
                result.append(GTFLine(*gtf_line))

            return result

    def __str__(self):
        return "\n".join(map(lambda x: str(x), self.gtf_files))

    def export_by_name_and_counter(self, gtfname, counter):

        exports = []
        # TODO: should be only one line ???
        # What if somebody would like map both (single and multimapping) we can use different output directories for such mapping

        # get gtfline
        try:
            gtfline = list(filter(lambda x: x.ANNOTATION_FILE == gtfname and x.COUNTER == counter, self.gtf_files))[0]
        except IndexError:
            print(f"Such ANNOTATION_FILE = {gtfname} or COUNTER = {counter} are not exist in GTF config file")
            exit(1)

        # TODO: Add path to dir with GTF files
        # create list of exports
        exports = "\n".join([f"export {e}={getattr(gtfline, e)}" for e in gtfline._fields])
        return exports

    def gtfNameCounterToBash(self):
        result = []
        for gtf in self.gtf_files:
            result.append(gtf.ANNOTATION_FILE)
            result.append(gtf.COUNTER)
        return " ".join(result)

    def gtf_by_de_index(self, deix):
        filtered_gtf = filter(lambda g: g.DE_INDEX == deix, self.gtf_files)
        unique_gtf = set(map(lambda x: (x.ANNOTATION_FILE, x.COUNTER), filtered_gtf))
        flatten_set=list(sum(unique_gtf, ()))
        return " ".join(flatten_set)

    def counter_by_de_index(self, deix):
        filtered_gtf = filter(lambda g: g.DE_INDEX == deix, self.gtf_files)
        unique_counter = set(map(lambda x: x.COUNTER, filtered_gtf))
        list_of_counters = list(unique_counter)
        # TODO: if list contains more than 1 value --> error
        return list_of_counters[0]


class EnvInterpolation(configparser.ExtendedInterpolation):
    """Interpolation which expands environment variables in values."""

    def before_get(self, parser, section, option, value, defaults):
        value = super().before_get(parser, section, option, value, defaults)
        return os.path.expandvars(value)


class SystemConfig:

    def __init__(self, config_path):
        generate_example(config_path, DEFAULT_PIPELINE_CONFIG)
        self.current_path = os.path.dirname(config_path)

        # read config with case sensitive options
        config = configparser.RawConfigParser(interpolation=EnvInterpolation())
        config.optionxform = lambda option: option
        config.read(config_path)
        self.config = config

    def setup_env_variables(self):
        exports = []

        for section in self.config.sections():
            keys = list(self.config[section].keys())
            for key in keys:
                exports.append((key, self.config[section][key]))

        exports.append(('SETUP_PIPELINE_DIR', self.current_path))
        exports = [f"export {e[0]}={e[1]}" for e in exports]
        return "\n".join(exports)


class DiffExpression():
    """
    Generate folders for step 9a,9b,9c according to configuration and
    TEMPLATES directories
    """
    def __init__(self, templates, comp_config, sample_config, sys_config, gtf_config):
        # check template path
        self.templates = templates
        self.comp_config = comp_config
        self.sample_config = sample_config
        self.sys_config = sys_config
        self.gtf_config = gtf_config

    def create_dir(self, dirpath):
        try:
            os.mkdir(dirpath)
        except OSError:
            print(f"Creation of the directory {dirpath} failed ")
            exit(1)

    def remove_dir(self, dirpath):
        try:
            shutil.rmtree(dirpath)
        except OSError:
            print(f"Removing of the directory {dirpath} failed. Try to remove it manually.")
            exit(1)

    def copy_dir(self, src, dest):
        try:
            copy_tree(src, dest)
        except OSError:
            print(f"Copying of the directory {src} failed ")
            exit(1)

    def write_condition(self, path, sample_config, compar_config, condition_num):
        """
        Write files Condition_{1,2}
        return Condition_Name
        """

        cond_value = getattr(compar_config, condition_num)
        cond_samples = sample_config.samplesByGroup(cond_value)
        condition_name = cond_samples[0].Condition_Name
        
        with open(path+"/"+condition_num+".txt", "w") as Condition:
            header = "Sample_ID\tDescription\n"
            Condition.write(header)
            samples = [f"{s.Sample_ID}\t{s.Condition_Name}" for s in cond_samples]
            samples = "\n".join(samples)+"\n"
            Condition.write(samples)
        return condition_name

    def clean(self, de_path="./"):
        de_dirs = glob.glob(de_path+"/09*")
        for d in de_dirs:
            print(f"Removing {d}")
            self.remove_dir(d)

    def replace_in_file(self, fname, list_of_replacements):
        with fileinput.input(fname, inplace=True) as runjobs:
            for line in runjobs:
                found = False
                for rp in list_of_replacements:
                    if rp[0] in line:
                        found = True
                        newline = line.replace("TEMPLATE", rp[1])
                        print(newline, end='')
                if not found:
                    print(line,  end='')

    # checking that the groups in the comparisons are a subset of sample_labels.txt
    def check_consistency(self):
        return self.comp_config.groups().issubset(self.sample_config.groups())
                    
    def generate(self, output_location="./"):
        # output_location start from SETUP_PIPELINE_DIR
        # set(map(lambda x: x.DE_INDEX,self.gtf_config.gtf_files()))
        for CMP in self.comp_config.comparisons:
            for TMPL in self.templates:
                for DEIX, GTFID, COUNTER in set(map(lambda x: (x.DE_INDEX, x.OUTPUT_DIR.split("_")[0], x.COUNTER), self.gtf_config.gtf_files)):
                    DIR_NAME = os.path.basename(TMPL).replace('TEMPLATE_', '')
                    DIR_NAME = DIR_NAME.replace('NN', CMP.Comparison_Number)
                    DIR_NAME = DIR_NAME.replace('II', DEIX)
                    DIR_NAME = DIR_NAME.replace('GTFID', GTFID)

                    # SCRIPT_PATH = os.getenv('SETUP_PIPELINE_DIR')

                    SCRIPT_PATH = os.getcwd()

                    # TODO: change it to OUTPUT path (SCRIPT_PATH+"/"+)

                    DIR_PATH = output_location + "/" + DIR_NAME

                    print(f"Creating {DIR_PATH}")
                    # copy template dirs from TMPL to DIR_PATH
                    self.copy_dir(TMPL, DIR_PATH)

                    # go to inside diff.ex task directory
                    os.chdir(DIR_PATH)

                    # write Condition_{1,2} files, return condition_name
                    cond_name_1 = self.write_condition(DIR_PATH, sample_config, CMP, "Condition_1")
                    cond_name_2 = self.write_condition(DIR_PATH, sample_config, CMP, "Condition_2")

                    # inplace changing TEMPLATE files

                    replacement_list = [
                        ("CONDITION_1_NAME=TEMPLATE", cond_name_1),
                        ("CONDITION_2_NAME=TEMPLATE", cond_name_2),
                        ("DE_INDEX=TEMPLATE", DEIX),
                        ("COMPAR_NUM=TEMPLATE", str(CMP.Comparison_Number)),
                        ("COUNT_PROGRAM=TEMPLATE", COUNTER)
                    ]

                    self.replace_in_file("Run_Jobs.sh", replacement_list)
                    self.replace_in_file("Summarize_Jobs.sh", replacement_list)

                    # with fileinput.input("setup_DiffExp.sh", inplace=True) as runjobs:

                    #     for line in runjobs:
                    #         if cond_1_tmpl in line:
                    #             newline = line.replace("TEMPLATE", f'"{cond_name_1}"')
                    #             print(newline, end='')
                    #         elif cond_2_tmpl in line:
                    #             newline = line.replace("TEMPLATE", f'"{cond_name_2}"')
                    #             print(newline, end='')
                    #         elif compar_num in line:
                    #             newline = line.replace("TEMPLATE", str(CMP.Comparison_Number))
                    #             print(newline, end='')
                    #         else:
                    #             print(line,  end='')

                    # go back to script directory
                    os.chdir(SCRIPT_PATH)


def generate_example(config_path, default_config):
    if not os.path.exists(config_path):
        print(f"Generating example config file: {os.path.basename(config_path)}")

        import getpass
        import os.path as op
        
        # change user in Pipeline_Setup.conf
        default_config=default_config.replace("CHANGE_USER_NAME", getpass.getuser())

        # change DATASET_DIR in Pipeline_Setup.conf
        PROJPATH_BASE=op.abspath(op.join(__file__, op.pardir, op.pardir, op.pardir))
        # add SAMPLES subdir for BAM files
        PROJPATH=PROJPATH_BASE+"/SAMPLES"
        default_config = default_config.replace("PATH_TO_PROJECT", PROJPATH)

        # change project name
        PROJECT_LABEL=os.path.basename(PROJPATH_BASE)
        default_config = default_config.replace("DEFAULT_LABEL", PROJECT_LABEL)
        
        with open(config_path, "w") as conf:
            conf.write(default_config)
            

# function which parse FASTQ index files
# INPUT: index_file - file with FASTQ indexes
# INPUT: sample_id (ex. G186_M1)
# OUTPUT: None (if found nothing) otherwise (PROJECT_NAME, READ1,
# READ2), can produce exception
# FASTQ index simple text format allows to store information about
# sample:
# PRJ,PROJECT_NAME1,/PROJECT_PATH1
# SAMPLE_ID1, FASTQ_WITH_READS1, FASTQ_WITH_READS2
# SAMPLE_ID2, FASTQ_WITH_READS1, FASTQ_WITH_READS2
# PRJ,PROJECT_NAME2,/PROJECT_PATH2
# SAMPLE_ID3, FASTQ_WITH_READS1, FASTQ_WITH_READS2
# SAMPLE_ID4, FASTQ_WITH_READS1, FASTQ_WITH_READS2
# ...
# each SAMPLE_ID must have unique name

def find_in_index(index_file, sample_id_param):

    sample_dict = {}
    with open(index_file, "r") as f:
        project_name = ""
        project_path = ""
        for line in f:
            # TODO: check if line is appropriate for format
            if len(line.split(",")) != 3:
                continue
            
            # if header line
            if line.startswith("PRJ"):
                header = line.strip().split(",")
                header = list(filter(lambda x: x != '', header))
                project_name, project_path = header[1], header[2]
                if not project_path.endswith('/'):
                    project_path = project_path+"/"
                continue
            
            # if sample line
            sample = list(filter(lambda y: y != '', map(lambda x: x.strip(), line.strip().split(','))))
            sample_id, read1, read2 = sample[0],sample[1],sample[2]

            # set some boundaries for the sample_id
            # only one underscore for sample id (ex. G186_M1)
            if len(sample_id.split("_")) != 2:
                raise ValueError(f"ERROR: {sample_id} sample name should contain only one underscore character, like G186_M1 (PROJECT_SHORTNAME)")
            
            # only unique sample_id in dict
            if sample_id in sample_dict:
                raise ValueError(f"ERROR: {sample_id} already in the index file ({index_file}). All SAMPLE_IDs should be unique")
            sample_dict[sample_id] = [project_name, project_path+read1, project_path+read2]

    if sample_id_param not in sample_dict:
        return None
    else:
        # return values, but strip before
        def only_printable(w):
            return "".join(list(filter(lambda x: x in string.ascii_letters + string.digits+"-_/.()", w)))
        return list(map(lambda x: only_printable(x), sample_dict[sample_id_param]))

    
# This function take local and global FASTQ index files and try to
# find required SAMPLE_ID inside. If local index file is presented
# then function tries to find SAMPLE_ID in this file first. If
# SAMPLE_ID was not found in local index then try to find it in the
# global index. If sample was not found in both indexes raise
# ValueError
# INPUT: sample_id
# INPUT: local index file path
# INPUT: global index file path
# OUTPUT: (PROJECT_NAME, /PATH/TO/READ1.f*q.gz, /PATH/TO/READ2.f*q.gz) or exception

def get_sample_info(SAMPLE_ID, USER_DEFINED_INDEX, DEFAULT_INDEX):

    # if USER_DEFINED_INDEX exist try to find SAMPLE_ID otherwise continue
    if USER_DEFINED_INDEX is not None and os.path.exists(USER_DEFINED_INDEX):
        result = find_in_index(USER_DEFINED_INDEX, SAMPLE_ID)
        # return result only if sample in local index
        if result is not None:
            sys.stderr.write(f"INFO: {SAMPLE_ID} found in {USER_DEFINED_INDEX}\n")
            return result
        
    # if DEFAULT_INDEX exist try to find SAMPLE_ID otherwise generate error 
    if os.path.exists(DEFAULT_INDEX):
        result = find_in_index(DEFAULT_INDEX, SAMPLE_ID)
        if result is not None:
            sys.stderr.write(f"INFO: {SAMPLE_ID} found in {DEFAULT_INDEX}\n")
            return result
        else:
            raise ValueError(f"ERROR: cannot find {SAMPLE_ID} sample in {DEFAULT_INDEX}")
    else:
        raise ValueError(f"ERROR: {DEFAULT_INDEX} file does not exist")
    
            
if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-e", "--export",
                        help="export all env.variables",
                        action="store_true")

    parser.add_argument("-g", "--generate",
                        help="generate Diff.Ex directories from templates",
                        action="store_true")

    parser.add_argument("-s", "--samples",
                        help="return bash array of triples (DIR,ID,DESCR)",
                        action="store_true")

    parser.add_argument("--samples_by_group",
                        nargs='+',
                        metavar=('group', 'color_enable'),
                        help="returns bash array of samples by group name (ex. A, B, ..), if specify second parameter like 'color_enable' the function returns 4-element bash array per sample with color in the last column - (DIR,ID,DESCR,COLOR)")
    
    parser.add_argument("-c", "--samples_with_color",
                        help="return bash array of (DIR,ID,DESCR,COLOR)",
                        action="store_true")

    parser.add_argument("--groups",
                        help="return bash array of all groups from Sample_Labels.txt",
                        action="store_true")

    parser.add_argument("-f", "--gtf_annotation_and_counter",
                        help="return bash array of pairs (GTF_ANNOTATION_FILE, COUNTER)",
                        action="store_true")

    parser.add_argument("-n", "--export_gtf_by_name_and_counter",
                        nargs=2,
                        metavar=('gtf_fname', 'counter'),
                        help="export GTF variables by gtf filename and counter")

    parser.add_argument("-i", "--gtf_by_DE_INDEX",
                        nargs=1,
                        metavar=('DE_INDEX'),
                        help="return bash array of gtf files which groupped by DE_INDEX (9a,9b,...)")

    parser.add_argument("-t", "--counter_by_DE_INDEX",
                        nargs=1,
                        metavar=('DE_INDEX'),
                        help="return counter by DE_INDEX (9a,9b,...)")

    parser.add_argument("--venn_number", 
                        help=f"return the number of Venn comparisons ({VENN_CONFIG})", 
                        action="store_true")
    
    parser.add_argument("--venn_comparisons_by_ix",
                        nargs=1,
                        metavar=('VENN_INDEX'),
                        help=f"return bash array of comparison numbers by venn index ({VENN_CONFIG})")

    parser.add_argument("--get_sample_info",
                        nargs=1,
                        metavar=('SAMPLE_ID'),
                        help="return (PROJECT, READ1, READ2) or SAMPLE_NOT_FOUND by SAMPLE_ID from index files")
    
    # gtf_by_annotation_and_counter
    args = parser.parse_args()

    # directory with configs and python script
    current_path = os.path.dirname(os.path.abspath(__file__))

    system_config = SystemConfig(current_path+"/"+PIPELINE_CONFIG)
    sample_config = SampleConfig(current_path + "/" + SAMPLES_CONFIG)
    comparison_config = ComparisonsConfig(current_path + "/" + COMPARISON_CONFIG)
    venn_config = VennConfig(current_path + "/" + VENN_CONFIG)
    # print(venn_config.get_venn_by_index(0))
    # print(venn_config.get_number_of_venns())

    # TODO: if config is not exist make error exit
    DEFAULT_GTF_CONFIG = system_config.config["SYSTEM"]["GTF_FILES_CONFIG"]
    # DEFAULT_GTF_CONFIG = current_path + "/GTFconfig.csv"
    
    gtf_config = GTFconfig(DEFAULT_GTF_CONFIG, current_path)

    if args.export:
        print(system_config.setup_env_variables())
        exit(0)
    elif args.generate:
        
        # by default all templates located in level up dir
        TEMPLATE_PATHS = glob.glob("../TEMPLATE*")

        diffex = DiffExpression(TEMPLATE_PATHS,
                                comparison_config,
                                sample_config,
                                system_config,
                                gtf_config)
        
        DE_DIR_PATH = os.path.dirname(TEMPLATE_PATHS[0])

        # checking groups consistency, exit if not
        if not diffex.check_consistency():
            print(f"Groups in config files ({DEFAULT_SAMPLES_CONFIG} and {DEFAULT_COMPARISON_CONFIG}) are not consistent")
            exit(1)
        
        # Remove DE directories before generation
        diffex.clean(DE_DIR_PATH)

        # Output dir by default the same as TEMPLATES_PATH
        diffex.generate(DE_DIR_PATH)
        print("Diff.ex directories are generated")
        exit(0)
        
    elif args.samples:
        print(sample_config.samplesToBash())
        exit(0)

    elif args.samples_by_group:
        if len(args.samples_by_group) == 1:
            group = args.samples_by_group[0]
            color = False
        else:
            group = args.samples_by_group[0]
            color = True
        print(sample_config.samplesToBash(group, color))
        exit(0)

    elif args.samples_with_color:
        print(sample_config.samplesWithColorToBash())
        exit(0)

    elif args.groups:
        print(" ".join(sample_config.groups()))
        exit(0)

    elif args.gtf_annotation_and_counter:
        print(gtf_config.gtfNameCounterToBash())
        exit(0)
        
    elif args.export_gtf_by_name_and_counter:
        a = args.export_gtf_by_name_and_counter
        print(gtf_config.export_by_name_and_counter(a[0], a[1]))
        exit(0)

    elif args.gtf_by_DE_INDEX:
        print(gtf_config.gtf_by_de_index(args.gtf_by_DE_INDEX[0]))

    elif args.counter_by_DE_INDEX:
        print(gtf_config.counter_by_de_index(args.counter_by_DE_INDEX[0]))

    elif args.venn_number:
        print(venn_config.get_number_of_venns())

    elif args.venn_comparisons_by_ix:
        a = args.venn_comparisons_by_ix
        print(venn_config.get_venn_by_index(int(a[0])))

    elif args.get_sample_info:
        sample_id = args.get_sample_info[0]

        # check how many indexes
        from glob import glob
        any_index = glob(current_path+"/*index.csv")
        user_defined_index = None
        if len(any_index) == 1:
            user_defined_index = any_index[0]
        elif len(any_index) > 1:
            import os.path as p
            fnames = list(map(lambda x: p.basename(x)), any_index)
            raise ValueError(f"ERROR: too many index files in the setup directory: {fnames}, only 0 or 1 should be provided")
        else:
            user_defined_index = None
            
        default_index = system_config.config["SYSTEM"]["FASTQ_DEFAULT_INDEX"]
        sample_info = get_sample_info(sample_id, user_defined_index, default_index)

        print(" ".join(sample_info))
        # try:
        #     res = " ".join(get_sample_info(sample_id, user_defined_index, default_index))
        # except Exception as e:
        #     print(e)
        #     raise e
        
        # print(res)
        exit(0)
    else:
        # parser.print_help()
        print("All config files are created!")


