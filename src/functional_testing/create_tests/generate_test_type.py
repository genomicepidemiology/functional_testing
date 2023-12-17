from helper_scripts.read_json_dic import check_output
from ft_base_init.create_test_base import CreateOutputs
from pandas import DataFrame
import json 
import os

class GenerateTest:
    def __init__(self, YourInput, TestCreator, cmd, dir_files_dic) -> None:
        self.YourInput = YourInput
        
        self.TestCreator = TestCreator
        
        self.cmd = cmd
        
        self.dir_files_dic = dir_files_dic
        
        self.test_active_json_out = YourInput.test_active_json_output
    
    @staticmethod
    def read_json(output_json):
        assert os.path.isfile(output_json), "Json output base file does not exist."
        with open(output_json, "r") as json_file:
            json_dic = json.load(json_file)
        json_base_df = DataFrame(json_dic["seq_variations"])
        return json_base_df
    
    def create_tests(self, output_json_base, output_json_new, CreateTestCommand):
        self.TestCreator.output_json_file = output_json_new
        self.TestCreator.software_command = CreateTestCommand.Command
        self.TestCreator.add_initial_command(output_json_base)
        self.TestCreator.test_base_dic()

     
    def create_all_tests(self, ):
        for basename, files in self.dir_files_dic.items():
            for mut_set in ["-c", "-acq", "-u"]:
                for species in [ "ecoli","cjejuni","ccoli", "senterica" ]:
                    CreateOutput = CreateOutputs(self.cmd)
                    output_json_base = CreateOutput.get_command(basename,
                                            filepath = files,
                                            species = species,
                                            mut_setting = mut_set,
                                            )
                    
                    CreateOutput.run_command()
                    
                    CreateTestCommand = CreateOutputs(self.cmd)
            
                    output_json_new = CreateTestCommand.get_command(basename,
                                                filepath = files,
                                                species = species,
                                                mut_setting = mut_set,
                                                output_dir_json = self.test_active_json_out)
                    
                    
                    if check_output(output_json_base) > 0:
                        self.TestCreator.species = species
                        
                        self.TestCreator.switch = CreateTestCommand.CommandRun.jump_to_dir # very important: tells TestCreator if directory for calling commnad needs to be changed
                        
                        self.create_tests(output_json_base, output_json_new, CreateTestCommand)
                        
                        self.TestCreator.test_seq_region()
                        
                        json_base_df = self.read_json(output_json_base) # has to be the one from original 
                        
                        
                        # both depend on json having a seq_variations content
                        self.TestCreator.test_mutations(mut_type = mut_set)
                        
                        if "-ifa" in CreateOutput.Command and len(json_base_df) > 0:
                            fasta_file = files[0]
                            self.TestCreator.add_contig_test(filepath_fasta = fasta_file)

                        self.TestCreator.remove_file_temp()
        
    def species_test(self, species_list, argument_identifier = "s"):
        for basename, files in self.dir_files_dic.items():
            for species in species_list:
                CreateOutput = CreateOutputs(self.cmd)
                output_json_base = CreateOutput.get_command(basename,
                                                            filepath = files)
                CreateOutput.Command.add_arg(arg = argument_identifier,
                                              value = species)
                CreateOutput.run_command()
                
                CreateTestCommand = CreateOutputs(self.cmd)
                
                output_json_new = CreateTestCommand.get_command(basename,
                                                                filepath = files,
                                                                output_dir_json = self.test_active_json_out)
                
                CreateTestCommand.Command.add_arg(arg = argument_identifier,
                                                    value = species)
                
                if check_output(output_json_base) > 0:
                    self.TestCreator.species = species
                    self.create_tests(output_json_base, output_json_new, CreateTestCommand)
                    self.TestCreator.test_seq_region()
                    

    def mutation_test(self, mut_setting_list):
        for basename, files in self.dir_files_dic.items():
            for mut_set in mut_setting_list:
                CreateOutput = CreateOutputs(self.cmd)
                output_json_base = CreateOutput.get_command(basename,
                                                            filepath = files,
                                                            mut_setting = mut_set)

                CreateOutput.run_command()
                
                CreateTestCommand = CreateOutputs(self.cmd)
                
                output_json_new = CreateTestCommand.get_command(basename,
                                                                filepath = files,
                                                                mut_setting = mut_set,
                                                                output_dir_json = self.test_active_json_out)
                
                
                if check_output(output_json_base) > 0:
                    self.TestCreator.species = "all"
                    self.create_tests(output_json_base, output_json_new, CreateTestCommand)
                    self.TestCreator.test_mutations(mut_type = mut_set)

    def custom_test(self, possible_values_list, argument_identifier):
        for basename, files in self.dir_files_dic.items():
            for input_value in possible_values_list:
                CreateOutput = CreateOutputs(self.cmd)
                output_json_base = CreateOutput.get_command(basename,
                                                            filepath = files)
                CreateOutput.Command.add_arg(arg = argument_identifier,
                                              value = input_value)
                CreateOutput.run_command()
                
                CreateTestCommand = CreateOutputs(self.cmd)
                
                output_json_new = CreateTestCommand.get_command(basename,
                                                                filepath = files,
                                                                output_dir_json = self.test_active_json_out)
                
                CreateTestCommand.Command.add_arg(arg = argument_identifier,
                                                    value = input_value)
                
                if check_output(output_json_base) > 0:
                    self.TestCreator.species = "all"
                    self.create_tests(output_json_base, output_json_new, CreateTestCommand)
                    self.TestCreator
    
    def save_tests(self, filename, filedir):
        self.TestCreator.write_markdown(filedir, filename)
        
    