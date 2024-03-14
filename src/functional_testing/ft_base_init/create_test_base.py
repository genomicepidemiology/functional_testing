import subprocess
from functional_testing.ft_base_init.create_command import SoftwareCommand
from functional_testing.helper_scripts.create_dir_struc import DirController
import sys
import shutil
from functional_testing.helper_scripts.switch_working_dir import SwitchDir
## Input: ldirectory to files, database

class CreateOutputs:
    def __init__(self, command = None):
        """
        :param files_dir: directory to files
        :param command: Command for running the program, e.g. 'resfinder'. None is default and will create the command autoamtically. If 
        """

        self.species = ["cjejuni", "ccoli", "ecoli", "senterica"]
        self.DirController = DirController()
        self.CommandRun = SwitchDir(command)
    #    self.check_kma_blast()
        self.command = command
    
    @staticmethod
    def check_individual_command(command):
        """
        :param command: Splitted command as a list with argumetns and parameters
        """
        if type(command) == str:
            command = command.split(" ")
        else:
            pass
        if "-s" in command:
            index_s = command[command.index("-s") + 1]
            command.pop(index_s)
            command.pop(index_s - 1)
        elif "--species" in command:
            index_s = command[command.index("--species") + 1]
            command.pop(index_s)
            command.pop(index_s - 1)
        
        return command
        
    
    def get_command(self, basename, filepath = None, species = "ecoli", output_dir_json = None):
        assert species in self.species, "Species not in database. Please choose one of the following: {}".format(self.species)
        self.Command = SoftwareCommand(self.DirController)
        if self.command == None: # FOr default which is just a resfinder run
            self.Command.base_command(filepath)
            self.Command.add_species(species)
            self.Command.add_mut_setting()
        else:
            # for case where you input a command 
            self.command = self.command.split(" ") # reformat in list object

            self.command = self.check_individual_command(self.command)
            self.Command.add_list(self.command)
            
        self.Command.add_file_input(filepath)
        self.Command.add_output()
        output_json_file = self.Command.add_json_output(species = species,basename=basename,
                                                        output_dir = output_dir_json)
        self.command = self.Command
        return output_json_file

            
        
    @staticmethod
    def check_kma_blast():
        try:
            subprocess.run(['kma', '-h'])
        except: 
            FileNotFoundError 
            sys.exit("KMA is not found. Please install KMA and add it to your PATH variable.")
            
        try:
            subprocess.run(['blastn', '-h'])
        except: 
            FileNotFoundError 
            sys.exit("KMA is not found. Please install KMA and add it to your PATH variable.")
        
    def run_command(self):
        self.CommandRun.cmd = self.command
        self.CommandRun.switcher() # runs command under consideration that tool directory is somewhere else and if so switches to that and back to old working dir
        shutil.rmtree(self.DirController.create_test_temp())
        
    