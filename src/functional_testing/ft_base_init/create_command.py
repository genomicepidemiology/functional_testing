
import os
from resfinder import run_resfinder

class SoftwareCommand(list):
    def __init__(self, DirController):
        self.DirController = DirController
           
    def add_list(self, new_list):
        """
        :param list: List. Contains either one or two elements
        """
        self.extend(new_list)
    
    def base_command(self, filepath):
        """
        :param path_resfinder: path to resfinder module
        :param filepath: List. Contains either one or two elements
        :species: Default is None. If not None than the species name will be added. Otherwise default of program will be chosen
        """

        filename0 = filepath[0]
        while True:
            filename0 = os.path.dirname(filename0)
            if os.path.isdir(os.path.join(filename0, "src")):
                path_resfinder = os.path.join(filename0, "src", "resfinder", "run_resfinder.py")
                break
            elif filename0 == "/":
                raise FileNotFoundError("No module found. You must have have a directory named src in one of the sub directories of your input files.")
            else:
                pass
        
        self.append("python")
        
        self.append(run_resfinder.__file__)
        
    def add_species(self, species):
        self.append("-s")
        self.append(species)
        
    def add_file_input(self, filepath):
        match  = [end for end in ["fastq", "fq"] if end in os.path.basename(filepath[0])]

        if len(match) > 0:
            self.append("-ifq")
            for i in filepath:
                self.append(i)
            
        else:
            filepath = filepath[0]
            self.append("-ifa")
            self.append(filepath)

    
    def add_mut_setting(self, setting):
        if "-" not in setting:
            setting = "-" + setting
        if not setting in ["-c", "-acq", "-u"]: 
            print("setting for mutation has not been provided. Mutation parameter is automatically set to -c.")
            setting = "-c"
        else:
            pass

            
        if setting == "-u":
            self.append(setting)
            self.append('-acq')
        else:
            self.append(setting)
        
    def add_output(self):
        temp_dir = self.DirController.create_test_temp()
        self.append("-o")
        self.append(temp_dir)
        
    def add_arg(self, arg, value):
        """
        :param arg: Argument for module
        :param value: Value for argument
        """
        self.append("-" + arg)
        if value != "":
            self.append(value)
        
    def add_json_output(self, species, basename = None, output_dir = None):
        self.species_dir = self.DirController.create_org_dir(species, basename)
        self.append("-j")
        if output_dir is None:
            matching_elements = [elem2 for elem1 in ["-c", "-acq", "-u"]  for elem2 in self if elem1 == elem2]
            if matching_elements:
                setting = "".join(matching_elements)
                output_json_file = os.path.join(self.species_dir, f"{basename}_{setting}.json")
            else:
                output_json_file = self.species_dir +".json"
            self.append(output_json_file)
        else:
            output_json_file = os.path.join(output_dir, basename + ".json")
            self.append(output_json_file)
        return output_json_file
    



