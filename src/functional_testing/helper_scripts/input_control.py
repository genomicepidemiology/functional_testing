
from argparse import ArgumentParser
import os
import glob
from ..ft_base_init.filenames_sorting import sort_filenames
import subprocess, shlex

class Input:
    def __init__(self) -> None:
        self.parser = ArgumentParser()
        self.test_active_json_output = self.check_tmp_dir()
        
    def enter_args(self):
        self.parser.add_argument("-d", "--test_files_dir",
                help="Directory with FASTA or FASTQ input files.",
                nargs="+",
                type = str,
                default = None
                )
        self.parser.add_argument("-t", "--type_test",
                help="Options are all, species, mutation, contig. Using this command, you can specify if you only want to create certain tests.",
                nargs="+",
                type = str,
                default = "all"
                )
        
        self.parser.add_argument("-sd", "--save_directory",
                help="Directory where the tests should be saved.",
                nargs="+",
                type = str,
                default = None
                )
        
        self.parser.add_argument("-f", "--filename",
                                 help="Filename for the test file.",
                nargs="+",
                type = str,
                default = None)
        
        self.parser.add_argument("-c", "--command",
                            help = "Command as string if you want to modify the command for testing your tool.",
                            nargs = "+",
                            type = str,
                            default = None
                            )
    
        self.parser.add_argument("-sl", "--species_list",
                                 help = "List of species you want to test.",
                            nargs = "+",
                            type = str,
                            default = None
                            )
    
    def parse_args(self):
        args = self.parser.parse_args()

        return args
    
    @staticmethod
    def check_python(python_interpreter):
        if not os.path.isfile(python_interpreter):  
            return False
        try:
        # Attempt to run a simple Python command
            subprocess.run([python_interpreter, '-c', 'print("Hello, Python!")'], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            return True
        except subprocess.CalledProcessError:
            return False
        
    @staticmethod
    def check_python_interpreter_location(module_path, python_path):
        file_path = module_path.replace('.', '/') + '.py'
        path = os.path.normpath(file_path)
        splitted_path = path.split(os.path.sep)
        src_index = splitted_path.index("src")
        module_name = splitted_path[src_index + 1]
        splitted_python = os.path.normpath(python_path).split(os.path.sep)
        assert module_name in splitted_python, "Please provide the python interpreter which is installed in the directory of the module you are calling." 
        assert splitted_path[0] == "src", "The module path you provide must start with src and is the relative path from the parent directory of the module you are calling."
    
    @staticmethod   
    def is_module(module_path):
        if not isinstance(module_path, str):
            return False

        # Split the input string into parts using dots
        parts = module_path.split('.')

        # Check if each part is a valid identifier
        for part in parts:
            if not part.isidentifier():
                return False

        # Check that there is at least one part
        if not parts:
            return False
        
        if module_path.endswith(".py"):
            return False

        return True
    

    
    def check_cmd(self, cmd):
        if cmd != None:
            args = shlex.split(cmd)
            python_interpreter = args[0]
            assert "-j" not in args, "Please do not provide the -j argument in the command. As this will be generated automatically to track the json file directory."
            assert "-json" not in args, "Please do not provide the -json argument in the command. As this will be generated automatically to track the json file directory."
            assert args[1] == "-m", "Please insert the command -m as a second argument."
            assert self.check_python(python_interpreter) == True, "Please provide a valid python interpreter."
            assert self.is_module(args[2]) == True, "Please provide a valid module name where the relative path from "
            self.check_python_interpreter_location(module_path = args[2], python_path=python_interpreter)
        else:
            pass   
            
    
    def check_receive_dir(self, dir_files):
        assert os.path.isdir(dir_files), "Directory with test files does not exist."
        filenames = glob.glob(os.path.join(dir_files, '*.fast*'))
        filenames = filenames + glob.glob(os.path.join(dir_files, '*.fa*'))
        filenames = filenames + glob.glob(os.path.join(dir_files, '*.fq*'))
        
        final_filenames = list(set(filenames))
        
        assert len(final_filenames) != 0, "No files found in directory."
        
        dir_files_dic = sort_filenames(final_filenames)
        return dir_files_dic
    
    
    
    
    @staticmethod
    def check_tmp_dir():
        test_active_json_out =  os.path.join(os.getcwd(), "functional_testing_temp")
        if not os.path.isdir(test_active_json_out):
            os.mkdir(os.path.join(os.getcwd(), "functional_testing_temp"))
        return test_active_json_out
        
        
        