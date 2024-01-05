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
            old_path = "'{}'".format(old_flag)
        print(old_path, new_flag)
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
    
