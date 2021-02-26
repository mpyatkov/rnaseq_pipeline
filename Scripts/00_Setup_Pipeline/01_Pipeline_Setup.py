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

PIPELINE_CONFIG = "Pipeline_Setup.conf"
SAMPLES_CONFIG = "Sample_Labels.txt"
COMPARISON_CONFIG = "Comparisons.txt"
VENN_CONFIG = "venn_comparisons.txt"

DEFAULT_VENN_CONFIG = """
venn_comparisons
""".strip()

DEFAULT_COMPARISON_CONFIG = """
Comparison_Number; Condition_1; Condition_2
1;                 A;           B
2;                 A;           C
""".strip()

DEFAULT_SAMPLES_CONFIG = """
Group;  Condition_Name;      Sample_DIR;  Sample_ID;  Description;                 Color
A;      E0771_Untreated_T72; G180_M1;     G180_M1;    E0771_Untreated_72h;         255,0,0
A;      E0771_Untreated_T72; G180_M2;     G180_M2;    E0771_Untreated_72h;         255,0,0
A;      E0771_Untreated_T72; G180_M3;     G180_M3;    E0771_Untreated_72h;         255,0,0
B;      E0771_4HC_T72;       G180_M4;     G180_M4;    E0771_5uM_4HC_72h;           0,255,0
B;      E0771_4HC_T72;       G180_M5;     G180_M5;    E0771_5uM_4HC_72h;           0,255,0
B;      E0771_4HC_T72;       G180_M6;     G180_M6;    E0771_5uM_4HC_72h;           0,255,0
C;      E0771_4HC_Ab_T72;    G180_M7;     G180_M7;    E0771_5uM_4HC_AntiIFNAR_72h; 0,0,255
C;      E0771_4HC_Ab_T72;    G180_M8;     G180_M8;    E0771_5uM_4HC_AntiIFNAR_72h; 0,0,255
C;      E0771_4HC_Ab_T72;    G180_M9;     G180_M9;    E0771_5uM_4HC_AntiIFNAR_72h; 0,0,255
""".strip()

DEFAULT_PIPELINE_CONFIG = """
[USER]
# Setup dataset directory
DATASET_DIR=PATH_TO_PROJECT

# Setup dataset label
DATASET_LABEL=DEFAULT_LABEL

# Setup BU login user name
BU_USER=CHANGE_USER_NAME

# Setup project name. Choices: (wax-dk,waxmanlab,wax-es)
PROJECT=wax-es

[STEPS]
# Bioanalyzer length (from Bioanalyzer tracings)
BIOANALYZER_LEN=330
# Usually the adaptor length is 60bp (60bp on each end: 120bp total)
ADAPTOR_LEN=60
# Input the adaptor length (bp) for one end
READ_LEN=150

# Depending on the library type, we can take into account information about
# specific strands. At the moment, the laboratory uses the dUTP preparation
# method, which corresponds to "fr-firststrand". Previous researchs could use
# other types of libraries, thus if you are unsure about library type just select
# STRANDEDNESS = 3. It will perform step 01 _..., which uses RSEQC
# (infer_experiment.py) and HISAT to detect a read strandedness from the small
# subset of reads for each sample.
# possible values:
# 0 - "unstranded"
# 1 - "firststrand"
# 2 - "secondstrand"
# 3 - "auto" (default)
STRANDEDNESS=3

# htseq only feature which allows to count different intersections for reads
# https://htseq.readthedocs.io/en/master/count.html
# 0="union"
# 1="intersection-strict"
# 2="intersection-nonempty"
MODE=2

# if you need tracks for samples and BigWig files use BIGWIG_ENABLE=1 in the
# option below
BIGWIG_ENABLE=1

# This configuration file contains all specific settings for GTF files
# used in the pipeline, like: name, enable/disable double counting,
# etc. For all options please visit GTF_FILES_DIR directory and check
# 9abcd.csv default config. The GTF_FILES_CONFIG parameter below
# combines GTF_FILES_DIR and DEFAULT_GTF_CONFIG
# Options: 9d.csv, 9abcd.csv
DEFAULT_GTF_CONFIG=9d.csv

# Default aligner, options:
# 0 - TopHat
# 1 - STAR
DEFAULT_ALIGNER=1

# DEFAULT_FC and DEFAULT_FDR options below impact only on PCA/correllation
# plots (Step 13)
DEFAULT_FC=2
DEFAULT_FDR=0.05

[SYSTEM]

# The location of the GTF files (by default you should have access to wax-es)
GTF_FILES_DIR=/projectnb/wax-es/routines/GTF_Files_default

# GTF configuration file which contains all information and options for GTF
# files used in this pipeline. If you would like to use custom configuration
# file just copy the default file in the current directory and assign the path
# to this file. Example: GTF_FILES_CONFIG=custom.csv
GTF_FILES_CONFIG=${SYSTEM:GTF_FILES_DIR}/${STEPS:DEFAULT_GTF_CONFIG}

# The location of the Bowtie2 indexes, at the moment it is only indexes for mouse mm9
# assembly
BOWTIE2INDEX_DIR=/projectnb/wax-es/routines/BowtieIndex

# Directory with full genomes
FASTA_DIR=/projectnb2/wax-es/routines/FASTA

# The location of the conda packages. This directory also contains scripts and
# configs and qsub file for setting up environments
CONDA_DIR=/projectnb/wax-es/routines/condaenv

# The location of the HISAT indexes, at the moment it is only indexes for mouse mm9
# assembly
HISAT2INDEX_DIR=/projectnb/wax-es/routines/hisat2index

# The location of the STAR indexes, at the moment it is only indexes for mouse mm9
# assembly (ExonCollapsed_76k GTF file)
STARINDEX_DIR=/projectnb/wax-es/routines/starindex_EC76K

# special directory which will contain FASTQC reports
VM_DIR_FASTQC=/net/waxman-server/mnt/data/waxmanlabvm_home/${USER:BU_USER}/FASTQC/${USER:DATASET_LABEL}

# special directory which will contain BIGWIG files
# subdirectory $VM_DIR_UCSC/COMMON for cram/bw/combined bw
# subdirectory $VM_DIR_UCSC/PERSONAL/USERNAME/DATASET_LABEL for track lines
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
            header_and_data = list(filter(lambda x: len(x) != 0, header_and_data))
            # print(header_and_data)
            header = header_and_data[0].replace(separator, ",")
            Sample = namedtuple('Sample', header)
            result = []
            for sample in header_and_data[1:]:
                sample = sample.split(separator)
                sample = list(map(lambda s: s.strip(), sample))
                result.append(Sample(*sample))
            return result

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
    
    def samplesToBash(self, group=None):
        # only Sample_DIR, Sample_ID, Description
        result = []
        if group:
            SAMPLES = self.samplesByGroup(group)
        else:
            SAMPLES=self.samples
            
        for sample in SAMPLES:
            result.append(sample.Sample_DIR)
            result.append(sample.Sample_ID)
            if group:
                result.append(sample.Condition_Name)
            else:
                result.append(sample.Description)
        return " ".join(result)

    def samplesWithColorToBash(self):
        #only Sample_DIR, Sample_ID, Description, Color
        result = []
        for sample in self.samples:
            result.append(sample.Sample_DIR)
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
            header_and_data = list(filter(lambda s: len(s) > 0, header_and_data))

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
            header_and_data = list(filter(lambda s: len(s) > 0, header_and_data))
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
            header_and_data = list(filter(lambda s: len(s) > 0, header_and_data))
            
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
            header = "Sample_DIR\tSample_ID\tDescription\n"
            Condition.write(header)
            samples = [f"{s.Sample_DIR}\t{s.Sample_ID}\t{s.Condition_Name}" for s in cond_samples]
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
        PROJPATH=op.abspath(op.join(__file__, op.pardir, op.pardir, op.pardir))
        default_config = default_config.replace("PATH_TO_PROJECT", PROJPATH)

        # change project name
        PROJECT_LABEL=os.path.basename(PROJPATH)
        default_config = default_config.replace("DEFAULT_LABEL", PROJECT_LABEL)
        
        with open(config_path, "w") as conf:
            conf.write(default_config)


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
                        nargs=1,
                        metavar=('group'),
                        help="return bash array of samples by group name (ex. A, B, ..)")
    
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
        print(sample_config.samplesToBash(args.samples_by_group[0]))
        # sample_config.samplesByGroup(args.samples_by_group[0])
        # print(" ".join(sample_config.samplesByGroup(args.samples_by_group[0])))
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
    else:
        # TODO: print only information that configuration files were generated
        # parser.print_help()
        print("All config files are created!")


