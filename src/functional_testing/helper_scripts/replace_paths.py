import os
import re
import glob


class ChangeMarkdown:
    def __init__(self, path_markdown):
        self.markdown_content = self.read_markdown(path_markdown=path_markdown)
        self.path_markdown = path_markdown
    @staticmethod  
    def read_markdown(path_markdown):
        """
        Static method that reads the content of the specified Markdown file.
        
        Parameters:
            path_markdown (str): Path to the Markdown file.
            
        Returns:
            str: Content of the Markdown file.
        """
        
        assert os.path.isfile(path_markdown) == True, "Markdown File does not exist."
        with open(path_markdown, 'r', encoding='utf-8') as file:
            markdown_content = file.read()
        return markdown_content

    def change_resfinder_path(self, old_path, new_path):
        """
        Changes the path of a ResFinder Python script in the Markdown content.
        
        Parameters:
            new_path (str): New path for the ResFinder Python script.
        """
        if new_path.startswith("'") and new_path.endswith("'"):
            pass
        else:
            new_path = "'{}'".format(new_path)
        if old_path.startswith("'") and old_path.endswith("'"):
            pass
        else:
            old_path = "'{}'".format(old_path)
        print(old_path, new_path)
        if new_path.endswith(".py") or new_path.endswith(".py'"):
            pass

        else:
            print("Please provide a path to a python file.")
            return
        self.markdown_content = self.markdown_content.replace(old_path, new_path)
     
    @staticmethod   
    def extract_substrings(original_string, search_string, end_character):
        pattern = search_string + r"(.*?)" + re.escape(end_character)
        substrings = re.findall(pattern, original_string)
        return substrings
    @staticmethod
    def extract_list_from_substring(substring):
        # Define a regular expression pattern to match the list
        pattern = r"\[.*?\]"
        
        # Use regular expression to find the list within the substring
        match = re.search(pattern, substring)
        
        if match:
            # Extract the matched substring (which is the list)
            list_string = match.group(0)
            
            # Evaluate the list string to convert it into a Python list
            result_list = eval(list_string)
            return result_list
        else:
            return None
    
    def remove_arg(self, arg_to_remove, substring_find = "subprocess"):
        """
        This command will look for a subprocess command in the string and then removes a specific argument from that command. This can be useful if your new program does not contain a certain flag anymore.
        
        Parameters:
            arg_to_remove (str): This is the argument you want to replace. It must start with - or -- .
            substring_find (str): This is the substring which acts as an identifier for the string you are looking for in the markdown document. For example subprocess looks for all substrings whcih starts with subprocess
        """
        all_subprocess = self.extract_substrings(self.markdown_content, substring_find, ")")
        for subprocess in all_subprocess:
            subprocess_new = subprocess
            command = self.extract_list_from_substring(subprocess) # extracts the list in subprocess 
            assert type(command) == list
            command_new = command.copy()
            if command == None:
                pass
            else:
                if arg_to_remove in command: # checks if arg is in the list
                    command_new.remove(arg_to_remove) # removes it if it is there
                else:
                    pass
                command_new = str(command_new)
                subprocess_new = subprocess_new.replace(str(command),command_new)
                self.markdown_content = self.markdown_content.replace(subprocess, subprocess_new)
        

        
        

    def change_data_paths(self, path_dir):
        """
        Changes the paths of data files in the Markdown content based on the provided directory.
        
        Parameters:
            path_dir (str): Path to directory with data files. E.g., files with .fastq or .fasta extension,
                           which can be found under the ResFinder tests.
        """
        
        assert os.path.isdir(path_dir) == True, "Directory does not exist."
        test_files = glob.glob(path_dir)
        for file in test_files:
            basename_file = os.path.basename(file)
            path_pattern_base = re.compile(r'(/[^/\s]+)+/' + re.escape(basename_file))
            self.markdown_content = re.sub(path_pattern_base, file, self.markdown_content)
        
    def change_json_temp(self, path_dir):
        """
        Changes the paths of JSON files in the Markdown content, specifically targeting directories with 'temp' in their names.
        
        Parameters:
            path_dir (str): Path to a directory containing JSON files, where the directory name includes the word 'temp'.
        """
        
        assert "temp" in path_dir, "Please provide a path with the word temp in one of the subdirectories"
        if not os.path.isdir(path_dir):
            os.mkdir(path_dir)
        test_files = glob.glob(path_dir)
        for file in test_files:
            basename_file = os.path.basename(file)
            pattern = re.compile(r'(/(?:[^/\s]+/)*[^/\s]*temp[^/\s]*)/' + re.escape(basename_file))
            self.markdown_content = re.sub(pattern, file, self.markdown_content)
            
    def change_json_true(self, path_dir):
        """
        Changes the paths of JSON files in the Markdown content, excluding directories with 'temp' in their names.
        
        Parameters:
            path_dir (str): Path to the directory containing JSON files.
        """
        
        assert os.path.isdir(path_dir) == True, "Directory does not exist. Please provide the correct directory with the json files."
        test_files = glob.glob(path_dir)
        for file in test_files:
            basename_file = os.path.basename(file)
            pattern = re.compile(r'(/(?:[^/\s]+/(?!.*temp)[^/\s]+/)*[^/\s]+/)' + re.escape(basename_file))
            self.markdown_content = re.sub(pattern, file, self.markdown_content)
            
            
    def change_flag(self, old_flag, new_flag):
        """
        Changes the flag of each command in the Markdown content.
        
        Parameters:
            old_flag (str): Old flag of the command.
            new_flag (str): New flag of the command.
        """
        if new_flag.startswith("'") and new_flag.endswith("'"):
            pass
        else:
            new_flag = "'{}'".format(new_flag)
        if old_flag.startswith("'") and old_flag.endswith("'"):
            pass
        else:
            old_flag = "'{}'".format(old_flag)
        if new_flag.endswith(".py") or new_flag.endswith(".py'"):
            pass
        self.markdown_content = self.markdown_content.replace(old_flag, new_flag)
        
    def save(self, filepath_new = None):
        if filepath_new == None:
            filepath_new = os.path.join(os.path.dirname(self.path_markdown), "new_" + os.path.basename(self.path_markdown))

        with open(filepath_new, 'w') as file:
            file.write(self.markdown_content)

def change_path(path_markdown, path_resfinder = None, path_dir_data = None,  path_json_true = None, path_json_temp = None):
    """
    :param path_markdown: Path to markdown file which should be changed.
    :param path_resfinder: Path to resfinder.py file.
    :param path_dir_data: Path to directory with data files. E.g. files with .fastq or .fasta extension, which can be find under the resfinder tests
    :param path_json_true: Path to directory with json files which are the base for the tests.
    :param path_json_temp: Path to directory with json files which are the base for the tests.
    :return: The modified markdown will be saved under: new_{basename of markdown file}
    """

    assert os.path.isfile(path_markdown) == True, "Markdown File does not exist."
    changer = ChangeMarkdown(path_markdown)
    if path_resfinder != None:
        changer.change_resfinder_path(path_resfinder)
    
    if path_dir_data != None:
        changer.change_data_paths(path_dir_data)
        
    if path_json_true != None:
        changer.change_json_true(path_json_true)

    if path_json_temp != None:
        changer.change_json_temp(path_json_temp)
        
    filepath_new = os.path.join(os.path.dirname(path_markdown), "new_" + os.path.basename(path_markdown))
    
    with open(filepath_new, 'w') as file:
        file.write(ChangeMarkdown.markdown_content)
    
changer = ChangeMarkdown("/home/people/s220672/resfinder/test_temp/test_new_test_all_formatted.md")
changer.remove_arg_in_subprocess("--point")
print(changer.markdown_content)