

## Test - Create only species tests

```
>>> import subprocess
>>> subprocess.run(["pdm", "run", "python", "src/functional_testing/__main__.py", "-d", "/home/people/s220672/resfinder/tests/data", "--type_test", "species", "--values_list", "ecoli"])
CompletedProcess(args=['pdm', 'run', 'python', 'src/functional_testing/__main__.py', '-d', '/home/people/s220672/resfinder/tests/data', '--type_test', 'species', '--values_list', 'ecoli'], returncode=0)

```