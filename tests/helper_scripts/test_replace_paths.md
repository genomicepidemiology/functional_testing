
## read_markdown

```
>>> from src.functional_testing.helper_scripts.replace_paths import ChangePath
>>> import os
>>> cwd = os.getcwd()
>>> md = os.path.join(cwd, "tests", "data", "functional_tests.md")
>>> markdownCP = ChangePath(md)
>>> markdown_content = markdownCP.read_markdown(md)
>>> markdown_content.startswith("# Imports")
True

```

## change_resfinder_path

```
>>> import random
>>> import string
>>> import re
>>> random_string = ''.join(random.choice(string.ascii_letters) for _ in range(12))
>>> markdownCP.change_resfinder_path(random_string)
Please provide a path to a python file.
>>> python_random = "asdf/" + random_string + ".py"
>>> pattern = r'([^/\s]+/[^/\s]+\.py)' 
>>> current_pattern = re.findall(pattern, markdown_content)[0]
>>> markdownCP.change_resfinder_path(python_random)
>>> new = re.findall(python_random, markdownCP.markdown_content)
>>> assert current_pattern != new[0], "The path was not changed."
>>> assert new[-1] == python_random, "The path was not changed to the correct path."

```

## change_data_paths

```
>>> 

```
