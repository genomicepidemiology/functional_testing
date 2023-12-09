import glob
import os
import subprocess
from .create_command import SoftwareCommand
from ..helper_scripts.create_dir_struc import DirController
import sys
import shutil
## Input: ldirectory to files, database

class CreateOutputs:
    def __init__(self, command = None):
        """
        :param files_dir: directory to files
        :param command: Command for running the program, e.g. 'resfinder'. None is default and will create the command autoamtically. If 
        """

        self.species = ["cjejuni", "ccoli", "ecoli", "senterica"]
        self.DirController = DirController()
    #    self.check_kma_blast()
        self.command = command
    
    
    def get_command(self, basename, filepath = None, species = "ecoli", mut_setting = "acq", output_dir_json = None):
        assert species in self.species, "Species not in database. Please choose one of the following: {}".format(self.species)
        self.Command = SoftwareCommand(self.DirController)
        if self.command == None: # FOr default which is just a resfinder run
            self.Command.base_command(filepath, species)
            self.Command.add_mut_setting(mut_setting)
        else:
            self.command = self.command.split(" ") # reformat in list object
            self.Command.add_list(self.command)
        
        self.Command.add_output()
        output_json_file = self.Command.add_json_output(species = species,basename=basename, output_dir = output_dir_json)
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
        subprocess.run(self.command)
        shutil.rmtree(self.DirController.create_test_temp())
        
    