

## Initialize Input

```
>>> from src.functional_testing.helper_scripts.input_control import Input
>>> import os
>>> YourInput = Input()
>>> test_temp_dir = YourInput.test_active_json_output
>>> os.path.isdir(test_temp_dir)
True

```

# check_python
You need to give a path to python interpreter yourself

```
>>> python_path = "/home/people/s220672/functional_testing/.venv/bin/python"
>>> YourInput.check_python(python_path)
True
>>> python_no = "/asdfjklasfd"
>>> YourInput.check_python(python_no)
False

```

## is_module

```
>>> YourInput.is_module("src.virulencefinder.__main__")
True
>>> YourInput.is_module("src.virulencefinder.__main__.py")
False

```

## check_python_interpreter_location

```
>>> python_path = "/home/people/s220672/virulencefinder/.venv/bin/python"
>>> YourInput.check_python_interpreter_location(module_path = "src.virulencefinder.__main__", python_path = python_path)

```
