from ..helper_scripts.read_json_dic import check_output
from ..ft_base_init.create_test_base import CreateOutputs





class GenerateTest:
    def __init__(self, YourInput, TestCreator, cmd, dir_files_dic) -> None:
        self.YourInput = YourInput
        
        self.TestCreator = TestCreator
        
        self.cmd = cmd
        
        self.dir_files_dic = dir_files_dic
        
        self.test_active_json_out = YourInput.test_active_json_output
            
    
    def create_tests(self, output_json_base, output_json_new, CreateTestCommand):
        self.TestCreator.output_json_file = output_json_new
        self.TestCreator.software_command = CreateTestCommand.Command
        self.TestCreator.add_initial_command(output_json_base)
        self.TestCreator.test_base_dic()

     
    def create_all_tests(self, ):
        for basename, files in self.dir_files_dic.items():
            for mut_set in ["-c", "-acq", "-u"]:
                for species in ["cjejuni", "ccoli", "ecoli", "senterica"]:
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
                        self.create_tests(output_json_base, output_json_new, CreateTestCommand)
                        
                        self.TestCreator.species = species
                        
                        self.TestCreator.test_seq_region()
                        
                        self.TestCreator.test_mutations(mut_type = mut_set)
                        
                        self.TestCreator.add_contig_test()
     
        
    def species_test(self, species_list, argument_identifier = "-s"):
        for basename, files in self.dir_files_dic.items():
            print(basename)
            for species in species_list:
                CreateOutput = CreateOutputs(self.cmd)
                output_json_base = CreateOutput.get_command(basename,
                                                            filepath = files)
                CreateOutput.Command.add_args(arg = argument_identifier,
                                              value = species)
                CreateOutput.run_command()
                
                CreateTestCommand = CreateOutputs(self.cmd)
                
                output_json_new = CreateTestCommand.get_command(basename,
                                                                output_dir_json = self.test_active_json_out)
                
                CreateTestCommand.Command.add_args(arg = argument_identifier,
                                                    value = species)
                
                if check_output(output_json_base) > 0:
                    self.TestCreator.species = species
                    self.create_tests(output_json_base, output_json_new, CreateTestCommand)
                    self.TestCreator.test_seq_region()
                    self.TestCreator.add_contig_test()

    def mutation_test(self, mut_setting_list, argument_identifier = "-c"):
        for basename, files in self.dir_files_dic.items():
            for mut_set in mut_setting_list:
                CreateOutput = CreateOutputs(self.cmd)
                output_json_base = CreateOutput.get_command(basename)
                CreateOutput.Command.add_args(arg = argument_identifier,
                                              value = mut_set)
                CreateOutput.run_command()
                
                CreateTestCommand = CreateOutputs(self.cmd)
                
                output_json_new = CreateTestCommand.get_command(basename,
                                                                output_dir_json = self.test_active_json_out)
                
                CreateTestCommand.Command.add_args(arg = argument_identifier,
                                                    value = mut_set)
                
                if check_output(output_json_base) > 0:
                    self.TestCreator.species = species
                    self.create_tests(output_json_base, output_json_new, CreateTestCommand)
                    self.TestCreator.test_seq_region()
                    self.TestCreator.add_contig_test()
        
    def save_tests(self, filename, filedir):
        self.TestCreator.write_tests(filedir, filename)
        
    