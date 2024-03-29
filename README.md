# functional_testing

## Purpose of this repository

This software was developed to create automatic functional tests for the softwares offers by the Center for genomic epidemiology at DTU. Initially, its purpose was to create automatic tests for resfinder which is why the default settings for it are designed for resfinder. However, the user can easily use this software for other softwares developed according to the standards of the institute. 



## How to use this software

### Installation

```bash
git clone https://github.com/genomicepidemiology/functional_testing.git

cd functional_testing

pdm install
```

### Make your first tests for resfinder 

```bash
pdm run python src/functional_testing/__main__.py -d /home/people/s220672/resfinder/tests/data  
``` 


### Flags

- `--test_files_dir`: This flag is *required*. Directory where the fastq and fasta files are located. Based on these files the software will create the tests for.
- `--type_test`: Default is all. This flag can be specified if only certain test types should be created. For instance you can only create tests for certain species if you add: --type_test species. Other options are:
<br>
*mutations*: Creates tests for different mutation types
<br>
*contigs*: Only tests on whether the software is stale for having the found gene on two contigs. 
<br>
*custom*: This creates tests for a customized input. You need to specify the flags arg_identifier and values_list for that.
<br>

- `--values_list`: Default is None. Here you can specify which values for your flag the tests should be written for. An example are the species where you can input multiple strings consecutively. For instance: --values_list ecoli ccoli 

- `--arg_identifier`: Default is None. Here you can specify the flag for which you want to create tests for. An example is the species flag. For instance: --arg_identifier species 
<br>
**NOTE**: Sometimes the flags between different softwares differ but they have the same meaning or do not have a different functionality. Then you can use this flag to input another flag for which a certain test should be created for. 

- `--save_directory`: Default is None. Here you can specify the directory where the tests should be saved. If you do not specify it, the tests will be saved in the directory where the software is located.

- `--filename`: Default is None. Here you can specify the name of the test file. If you do not specify it, the tests will be saved in the directory where the software is located.

- `--command`: Default is None which means that tests for Resfinder will be created. However, you can create tests also for other softwares by specifying the command for that. Note ,that you need to define the python interpreter of the corresponding software. 
<br>

### Custom testing - Enter your own command

This is a little more complex because your software will not be part of the module functional_testing. So we need to parse a default command of your old software as a string after the flag -c. You should consider the following:

1. Maybe you have changed the dependencies or similar in your new software. So you must parse a specific python interpreter from a virtual environment which was used or contains the dependencies for your **old** version. This ensures that no errors occur during the testing
2. You can parse the path to your main script of your old software version as module or as a script.
3. Add additional flags which are necessary for your software to run.
4. The type of the test should be custom if you will not analyse either contigs, mutations or species.
    - If you parse "custom", then a basic comparison of the dictionary will happen.

*Example*:
<br>
```bash
pdm run python -m src.functional_testing.__main__ -d /home/people/s220672/functional_testing/tests/data -c "/home/people/s220672/virulencefinder/.venv/bin/python -m src.virulencefinder.__main__ -p /home/people/s220672/databases/virulencefinder_db" --type_test "species" --values_list "virulence_ecoli" --arg_identifier "d"
```
The command flag is given as a string and contains the python interpreter from the virtual environment of virulencefinder. Further it has the relative path to the module based on the parent directory of virulencefinder. Note: YOu should put the old version The flag -p is added because it needs to be defined as an internal command for virulencefinder to be able to run. Every single flag which follows after that is a flag for the software *functional_testing*
<br>

## Output

The software outputs a markdown file which contains a series of tests based on the input files in the directory. If you move the output file to another directory you need to make sure that you can handle the imports in the markdown file. These are listed directly at the top.

**NOTE:** If the key seq regions in the output json file is empty when creating the tests, no tests will be created for the corresponding input file in the markdown document.
<br>

## Test information

The software creates various tests. Currenlty the following tests can be conducted:
- comparison of the whole json file. Checking whether the keys in both are the same

### Seq regions tests
- Checking whether found phenotypes are the same
- Checking whether the alignment length is the same
- Checking whether the reference start position is the same
- Checking whether the reference end position is the same

### Seq variations tests - mutations
- checking whether gene name is the same 
- Checking whether the reference start position is the same
- Checking whether the reference end position is the same
- Checking whether the reference codon is the same
- Checking whether the var codon is the same

### Seq variations tests - Contigs
Here, the input for the new run is first edited by a Fasta Handler which splits the Contig in the example fasta file in two contigs. Then the software is run on the file with two different contigs

- Checking whether the query end pos is the same
- Checking whether the query start pos is the same
- Checking whether the reference end pos is the same
- Checking whether the reference start pos is the same
 


### Example of a functional test

A functional test tests whether a new software version outputs the same results as the previous version. Therefore, a comparison between the value of a certain key in the json output will be conducted. An example for a test is:

``` python
assert seq_regions_table_true.loc['alignment_length'].values.tolist() == seq_regions_table_new.loc['alignment_length'].values.tolist(), 'Alignment length for species_name is not equal'
```

Here, the seq_regions_table_true was initially saved when creating the tests while the seq_regions_table_new is automatically generated in the test markdown file as a result of the new software version.



## A brief overview of how the software works

1. The software executes the commands of the old software version which is for instance in the package resfinder or the software parsed by the command
2. Then the software saves the output. This **MUST** be the standardized json file, since the software writes tests automatically based on the keys and values in this json file
3. The software writes a command for the new software version which maybe needs to be changed if it is created in a separated folder or similar. The output of this command is compared to the output of the old software version which was saved under functional_testing/functional_testing_data.
4. The software writes a markdown file whcih contains the tests.


## Modifying the created markdown script

It appears that sometimes you change crucial parts and identifiers of your software such as flags or paths. Therefore the tool provides commands to change these. 

**However, I think it is easier to just go to the test markdown file and replace the flags or similar with ctrl + F.**

*Example: Change flag -ifa to -f in all commands in the markdown file*

```python

from src.functional_testing.helper_scripts.replace_paths import ChangeMarkdown

changer = ChangeMarkdown(path_to_markdown_file)

changer.change_flag(old_flag = '-ifa', new_flag = '-f')

```

*Example: Change path of resfinder_python_script in the commands*

This might be useful if you did a major refactoring

```python

changer.change_resfinder_path(old_path = "/home/people/s220672/resfinder/src/resfinder/run_resfinder.py", new_path = "/home/people/s220672/resfinder/src/resfinder/__main__.py")

```

*Example: Remove an argument from a command in your markdown script*

This removes in all subprocess commands the argument --point. This is useful, if you deleted an argument in a new version of your program.

```python 

changer.remove_arg("--point", "subprocess")
```


### Save changes

```python

filepath_new = "my_new_markdown.md"
changer.save(filepath_new)

``` 
