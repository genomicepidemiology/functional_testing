from .helper_scripts.input_control import Input
from .create_tests.test_creator import MakeTests
from .create_tests.generate_test_type import GenerateTest
import sys
# problems: I needed to change the imports in Resfinder to have no module import error
# I could not run with pdm because of module import errors
# Create like a package to create tests for different params

# create a test for contigs - find the resistant gene and split it in two contigs, then test if the position is still the same - resfinder should be able to consider this


if __name__ == "__main__":
    YourInput = Input()
    
    YourInput.enter_args()
    
    args = YourInput.parse_args()
    
    
    dir_files = args.test_files_dir[0]
    cmd = args.command
    type_test = args.type_test[0]
    filename = args.filename
    dirname = args.save_directory
    
    dir_files_dic = YourInput.check_receive_dir(dir_files)
    
    print(dir_files_dic)
    
    test_active_json_out = YourInput.test_active_json_output
    
    YourInput.check_cmd(cmd)
    
    TestCreator = MakeTests() # initialize test creator
    
    TestCreator.add_imports()
    
    TestGenerator = GenerateTest(YourInput, TestCreator, cmd, dir_files_dic)

    if type_test == "all":
        TestGenerator.create_all_tests()
        
    elif type_test == "species":
        species_list = args.species_list
        if species_list == None:
            print("please provide a list of species you want to write the tests for. ")
            sys.exit(1)
        TestGenerator.species_test(species_list)
                    
    elif type_test == "mutation":
        TestGenerator.mutation_test()
    else:
        TestGenerator.contig_test()
        
    TestGenerator.save_tests(filename = filename, filedir = dirname)
        
    
    
    