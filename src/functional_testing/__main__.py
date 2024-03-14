from functional_testing.helper_scripts.input_control import Input
from functional_testing.create_tests.test_creator import MakeTests
from functional_testing.create_tests.generate_test_type import GenerateTest
from functional_testing.input_control.species_control import species_control
from functional_testing.input_control.mutation_control import mutations_control
# I could not run with pdm because of module import errors
# Create like a package to create tests for different params

# create a test for contigs - find the resistant gene and split it in two contigs, then test if the position is still the same - resfinder should be able to consider this



# you need to create a test for that the flags are the same for the testing software and the software to create the tests for
if __name__ == "__main__":
    YourInput = Input()
    
    YourInput.enter_args()
    
    args = YourInput.parse_args()
    
    dir_files = args.test_files_dir[0]
    cmd = args.command
    type_test = args.type_test
    filename = args.filename
    dirname = args.save_directory
    diff_arg = args.arg_identifier
    
    dir_files_dic = YourInput.check_receive_dir(dir_files)
    
    test_active_json_out = YourInput.test_active_json_output
    
    YourInput.check_cmd(cmd)
        
    TestCreator = MakeTests() # initialize test creator
    
    TestCreator.add_imports()
    
    TestGenerator = GenerateTest(YourInput, TestCreator, cmd, dir_files_dic)

    if type_test == "all":
        TestGenerator.create_all_tests()
        
    elif type_test == "species":
        species_list = species_control(args)
        if diff_arg != None:
            TestGenerator.species_test(species_list, argument_identifier = diff_arg)
        else:
            TestGenerator.species_test(species_list)    
    elif type_test == "mutation":
        mutation_list = mutations_control(args)
        TestGenerator.mutation_test(mut_setting_list=mutation_list)

    elif type_test == "contig":
        TestGenerator.contig_test()
    elif type_test == "custom":
        TestGenerator.custom_test(argument_identifier=diff_arg, 
                                  possible_values_list=args.values_list)
    else:
        raise ValueError("Type of test not recognized. Please choose one of the following: all, species, mutation, contig, custom.")
        
    TestGenerator.save_tests(filename = filename, filedir = dirname)
        
    
    
    