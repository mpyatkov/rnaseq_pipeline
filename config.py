from collections import namedtuple
import os
import sys
import configparser
from distutils.dir_util import copy_tree
from functools import reduce
import operator

class SampleConfig:
    separator = ";"
    def __init__(self, configpath):
        self.samples = self.__read_config(configpath, self.separator)
        
    def __read_config(self, configpath, separator):

        with open(configpath, "r") as config:
            header_and_data = list(map(lambda x: x.strip(), config.readlines()))
            header = header_and_data[0].replace(separator,",")
            Sample = namedtuple('Sample', header)
            result = []
            for sample in header_and_data[1:]:
                sample = sample.split(separator)
                result.append(Sample(*sample))
            return result

    def samplesByGroup(self,group):
        result = []
        for sample in self.samples:
            if sample.Group == group:
                result.append(sample)
        return result
                

    def groups(self):
        return set(map(lambda s: s.Group, self.samples))

    def samplesToBash(self):
        #only Sample_DIR, Sample_ID, Description
        result = [() for s in self.samples]
        result = []
        for sample in self.samples:
            result.append(sample.Sample_DIR)
            result.append(sample.Sample_ID)
            result.append(sample.Description)
        return " ".join(result)
        
    def __str__(self):
        return "\n".join(map(lambda x: str(x), self.samples))

class ComparisonsConfig:
    separator = ";"
    def __init__(self, configpath):
        self.comparisons = self.__read_config(configpath, self.separator)

    def __read_config(self, configpath, separator):
        with open(configpath, "r") as config:
            header_and_data = list(map(lambda x: x.strip(), config.readlines()))
            header = header_and_data[0].replace(separator,",")
            Comparison = namedtuple('Comparison', header)
            result = []
            for comparison in header_and_data[1:]:
                comparison = comparison.split(separator)
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
    
class SystemConfig:
    # read config file
    # export variables
    # copy directories

    def __init__(self, filename = "/Pipeline_Setup.conf"):
        self.current_path = os.path.dirname(os.path.abspath(__file__))
        os.environ['SCRIPT_DIR'] = self.current_path

        # read config with case sensitive option
        config = configparser.RawConfigParser()
        config.optionxform = lambda option: option
        config.read(os.getenv('SCRIPT_DIR')+"/"+filename)
        self.config = config

    def setup_env_variables(self):
        # for section in self.config.sections():
        #     x = list(self.config[section].keys())
        #     for y in x:
        #         print(y,self.config[section][y])
        # return f"export SCRIPT_DIR={self.current_path}"

        return f"export SCRIPT_DIR={self.current_path}\n"+f"export SCRIPT_DIR1={self.current_path}/asdf"

class DiffExpression():
    """
    Generate folders for step 9a,9b,9c according to configuration and
    TEMPLATES directories
    """
    def __init__(self, templates, comp_config, sample_config, sys_config):
        # check template path
        self.templates = templates
        self.comp_config = comp_config
        self.sample_config = sample_config
        self.sys_config = sys_config

    def create_dir(self, dirpath):
        try:
            os.mkdir(path)
        except OSError:
            print (f"Creation of the directory {path} failed ")
            exit(1)

    def remove_dir(self, dirpath):
        try:
            os.rmdir(path)
        except OSError:
            print (f"Deletion of the directory {path} failed ")
            exit(1)

    def copy_dir(self, src, dest):
        try:
            destination = copy_tree(src, dest)  
        except OSError:
            print (f"Copying of the directory {src} failed ")
            exit(1)

    def write_condition(self, path, sample_config, compar_config, condition_num):
        condition = [sample_config.samplesByGroup(s) for s in getattr(compar_config, condition_num).split(",")]
        condition = reduce(operator.iconcat, condition, [])
        condition.sort(key=lambda s: s.Sample_DIR)
        with open(path+"/"+condition_num+".txt", "w") as Condition:
            header = "Sample_DIR\tSample_ID\tDescription\n"
            print(header)
            Condition.write(header)
            samples = [f"{s.Sample_DIR}\t{s.Sample_ID}\t{s.Description}" for s in condition]
            samples = "\n".join(samples)
            Condition.write(samples)

    def generate(self):
        
        # get first template name
        # TODO: move dict to comparison class 
        counter_dict = {"HTSEQ":1, "FEATURECOUNT":1, "LNCRNA":1}
        
        for cm in self.comp_config.comparisons:
            # print(cm)
            # print(self.templates)
            TEMPLATE = list(filter(lambda t: cm.Counter in t, self.templates))[0]
            
            DIR_NAME = os.path.basename(TEMPLATE).replace('TEMPLATE_','')
            DIR_NAME = DIR_NAME.replace('NN', str(counter_dict.get(cm.Counter)))
            
            counter_dict[cm.Counter]+=1
            
            SCRIPT_PATH = os.getenv('SCRIPT_DIR')
            DIR_PATH = SCRIPT_PATH+"/"+DIR_NAME # TODO: change it to OUTPUT path
            print(SCRIPT_PATH)
            print(DIR_PATH)

            # create 
            self.copy_dir(TEMPLATE, DIR_PATH)
            
            # go to inside diff.ex task directory
            os.chdir(DIR_PATH)

            self.write_condition(DIR_PATH, sample_config, cm, "Condition_1")
            self.write_condition(DIR_PATH, sample_config, cm, "Condition_2")

            # go back to main directory
            os.chdir(SCRIPT_PATH)

                        
            # create CONDITION_{1,2}
            # change Run_Jobs....

            
            
    def regenerate():
        # remove all folders
        # generate
        pass

if __name__ == "__main__":
    system_config = SystemConfig()
    sample_config = SampleConfig(os.getenv('SCRIPT_DIR')+"/NSample_Labels.txt")
    comparison_config = ComparisonsConfig(os.getenv('SCRIPT_DIR')+"/group_comparison.txt")
    diffex = DiffExpression(["../TEMPLATE_09a_DiffExp_NN_HTSEQ",
                             "../TEMPLATE_09b_DiffExp_NN_FEATURECOUNTS",
                             "../TEMPLATE_09c_DiffExp_NN_LNCRNA"],
                            comparison_config,
                            sample_config,
                            system_config)

    diffex.generate()

    


    # if sys.argv[1] == "setup":
    #     a = SystemConfig()
    #     print(a.setup_env_variables())
    # else:
    #     print(sample_config.samplesToBash())



    
    # current_path = os.path.dirname(os.path.abspath(__file__))
    # sample_config = SampleConfig(current_path+"/NSample_Labels.txt")
    # print(f"export SCRIPT_DIR={current_path}")
    
    # print(sample_config)
    # print(sample_config.samplesByGroup('D'))
    # grs = sample_config.groups()


    # comparison_config = ComparisonsConfig("./group_comparison.txt")
    # print(comparison_config)
    # grc = comparison_config.groups()
    
    # print(grc == grs)

    # print(sample_config.samples[0][1])
