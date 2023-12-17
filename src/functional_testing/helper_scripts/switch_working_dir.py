import subprocess
import os



class SwitchDir:
    def __init__(self, cmd):
        """
        Initializes the SwitchDir instance with a given command and determines the directory to switch to.
        
        Parameters:
            cmd (str): The command to be executed.
        """
        self.cmd = cmd
        if self.cmd != None:
            self.jump_to_dir = self.find_jump_dir()
        else:
            self.jump_to_dir = None
            
    def switcher(self):
        """
        Executes the command after switching to the specified directory, if available; otherwise, executes the command in the current working directory.
        """
        if self.jump_to_dir != None:
            initial_wd = os.getcwd()
            #assert os.path.isdir(jump_to_dir), "jump_to_dir is not a directory. Something went wrong"
            os.chdir(self.jump_to_dir)
            subprocess.run(self.cmd)
            os.chdir(initial_wd)
        else:
            subprocess.run(self.cmd)
        
    def find_jump_dir(self,):
        """
        Determines the directory to switch to based on the provided command.
        
        Returns:
            str: The directory path.
        """
        tool_names = ["resfinder", "pointfinder", "virulencefinder", "plasmidfinder"]
        self.cmd = self.cmd.split(" ")
        base_path = [item for item in self.cmd for tool in tool_names if tool in item]
        if len(base_path) == 0:
            raise ValueError("No tool name found in command. Please provide a valid command.")
        else:
            base_path = base_path[0]
        base_path = os.path.normpath(base_path)
        base_path = os.path.dirname(base_path)
        jump_to_dir = None
        subdirs = base_path.split("/")
        item_path = ""
        for item in subdirs:
            item_path += item + "/"
            if item in tool_names:
                jump_to_dir = item_path
                break
            else: 
                pass
        return jump_to_dir

        