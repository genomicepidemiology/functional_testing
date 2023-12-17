
# imports


```
>>> import json
>>> import pandas as pd
>>> import subprocess

```

# Call subrocess

```
>>> subprocess.run(['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'cjejuni', '-acq', '-ifq', '/home/people/s220672/resfinder/tests/data/test_isolate_01_2.fq', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_01_2.json'])
CompletedProcess(args=['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'cjejuni', '-acq', '-ifq', '/home/people/s220672/resfinder/tests/data/test_isolate_01_2.fq', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_01_2.json'], returncode=0)

```


## Read json files from initial generated test

```
>>> json_file = '/home/people/s220672/functional_testing/functional_testing_data/cjejuni/test_isolate_01_2/test_isolate_01_2_-acq.json'
>>> with open(json_file, 'r') as f:
...		json_data_true = json.load(f)
>>> json_file_new = '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_01_2.json'
>>> with open(json_file_new, 'r') as f:
...		json_data_new = json.load(f)

```

### get keys command

``` 
>>> keys_true = list(json_data_true.keys())

```


### get seq_regions table
```
>>> seq_regions = json_data_true['seq_regions']
>>> seq_regions_table_true = pd.DataFrame(seq_regions)
>>> seq_region_new = json_data_new['seq_regions']
>>> seq_regions_table_new = pd.DataFrame(seq_region_new)

```

### test seq regions tests - phenotypes

``` 
>>> assert seq_regions_table_true.loc['phenotypes'].values.tolist() == seq_regions_table_new.loc['phenotypes'].values.tolist(), 'Phenotypes are not equal'

```

### test alignment length

```
>>> assert seq_regions_table_true.loc['alignment_length'].values.tolist() == seq_regions_table_new.loc['alignment_length'].values.tolist(), 'Alignment length for {self.species} is not equal'

```

### test reference start pos

```
>>> assert seq_regions_table_true.loc['ref_start_pos'].values.tolist() == seq_regions_table_new.loc['ref_start_pos'].values.tolist(), 'Query start for is not equal'

```

## Test seq variations

```
>>> seq_variations_true = json_data_true['seq_variations']
>>> seq_variations_table_true = pd.DataFrame(seq_variations_true)
>>> seq_variations_new = json_data_new['seq_variations']
>>> seq_variations_table_new = pd.DataFrame(seq_variations_new)

```

### test_seq variations - gene
this test will be done with a different file and different settings

```
>>> json_file = '/home/people/s220672/functional_testing/functional_testing_data/ecoli/test_isolate_05/test_isolate_05_-c.json'
>>> with open(json_file, 'r') as f:
...		json_data_new = json.load(f)
>>> seq_variations_new_file = pd.DataFrame(json_data_new['seq_variations'])
>>> seq_variations_new_file.loc['key'].values.tolist()[0].split(';')[0]
'gyrA'

``` 


