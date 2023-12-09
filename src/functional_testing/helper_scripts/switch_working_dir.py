import subprocess
import os

def switcher(jump_to_dir, cmd):
    initial_wd = os.getcwd()
    #assert os.path.isdir(jump_to_dir), "jump_to_dir is not a directory. Something went wrong"
    os.chdir(jump_to_dir)
    subprocess.run(cmd, shell = True)
    os.chdir(initial_wd)
    