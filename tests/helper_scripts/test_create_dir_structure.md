


## Initialize DirController

```
>>> from src.functional_testing.helper_scripts.create_dir_struc import DirController
>>> import shutil
>>> import os
>>> dc = DirController()
>>> temp = os.path.join(dc.cwd, "test_temp")
>>> os.path.exists(temp)
True
>>> shutil.rmtree(temp)
>>> os.path.exists(dc.data_storage)
True

```

## create_org_dir

```python
>>> species_dir = dc.create_org_dir("ecoli", "blabla")
>>> os.path.isdir(os.path.normpath(species_dir))
True
>>> shutil.rmtree(species_dir)

```
