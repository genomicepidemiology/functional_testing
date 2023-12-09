# Imports

```
>>> import subprocess
>>> import json
>>> import pandas as pd

```

## Software Command for ecoli

```
>>> subprocess.run(['python', '/home/people/s220672/functional_testing/env_ft/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-ifq', '/home/people/s220672/resfinder/tests/data/test_isolate_09b_1.fq', '-s', 'ecoli', '-c', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_09b_1.json'])
>>> json_file = /home/people/s220672/functional_testing/functional_testing_data/ecoli/test_isolate_09b_1/test_isolate_09b_1_-c.json
>>> with open(json_file, 'r') as f:
...		json_data_true = json.load(f)
>>> json_file_new = /home/people/s220672/functional_testing/functional_testing_temp/test_isolate_09b_1.json
>>> with open(json_file_new, 'r') as f:
...		json_data_new = json.load(f)

```

### Test keys and base structure of json dic

```
>>> assert json_file != json_file_new, 'Json files are the same'
>>> keys_true = list(json_data_true.keys())
>>> keys_new = list(json_data_new.keys())
>>> assert keys_true == keys_new, 'Keys are not equal'

```

### Test seq_region for ecoli

```
>>> seq_region = json_data_true['seq_region']
>>> seq_regions_table_true = pd.DataFrame(seq_region)
>>> seq_region_new = json_data_new['seq_region']
>>> seq_regions_table_new = pd.DataFrame(seq_region_new)
>>> assert seq_regions_table_true.loc['phenotypes'].values.tolist() == seq_regions_table_new['phenotypes'].values.tolist(), 'Phenotypes are not equal'
>>> assert seq_regions_table_true.loc['alignment_length'].values.tolist() == seq_regions_table_new['alignment_length'].values.tolist(), 'Alignment length for {self.species} is not equal'
>>> assert seq_regions_table_true.loc['query_start_pos'].values.tolist() == seq_regions_table_new['query_start_pos'].values.tolist(), 'Query start for {self.species} is not equal'
>>> assert seq_regions_table_true.loc['query_end_pos'].values.tolist() == seq_regions_table_new['query_end_pos'].values.tolist(), 'Query end for {self.species} is not equal'
>>> assert seq_regions_table_true.loc['ref_start_pos'].values.tolist() == seq_regions_table_new['ref_start_pos'].values.tolist(), 'Ref start for {self.species} is not equal'
>>> assert seq_regions_table_true.loc['ref_end_pos'].values.tolist() == seq_regions_table_new['ref_end_pos'].values.tolist(), 'Ref end for {self.species} is not equal'

```

# Test mutations for ecoli and mut_type: chromosomal

```
>>> seq_variations_true = json_data_true['seq_variations']
>>> seq_variations_table_true = pd.DataFrame(seq_variations_true)
>>> seq_variations_new = json_data_new['seq_variations']
>>> seq_variations_table_new = pd.DataFrame(seq_variations_new)
>>> assert seq_variations_table_new.loc['key'].values.tolist()[0].split(';')[0] == seq_variations_table_true.loc['key'].values.tolist()[0].split(';')[0], 'First element of key is not equal'
>>> assert seq_variations_table_new.loc['ref_start_pos'].values.tolist() == seq_variations_table_true.loc['ref_start_pos'].values.tolist(), 'Mutation start position for chromosomal and ecoli is not equal'
>>> assert seq_variations_table_new.loc['ref_end_pos'].values.tolist() == seq_variations_table_true.loc['ref_end_pos'].values.tolist(), 'Mutation end position for chromosomal and ecoli is not equal'
>>> assert seq_variations_table_new.loc['var_codon'].values.tolist() == seq_variations_table_true.loc['var_codon'].values.tolist(), 'Mutation codon for chromosomal and ecoli is not equal'
>>> assert seq_variations_table_new.loc['ref_codon'].values.tolist() == seq_variations_table_true.loc['ref_codon'].values.tolist(), 'Reference codon for chromosomal and ecoli is not equal'

```

## Software Command for ecoli

```
>>> subprocess.run(['python', '/home/people/s220672/functional_testing/env_ft/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-ifa', '/home/people/s220672/resfinder/tests/data/test_isolate_11.fa', '-s', 'ecoli', '-c', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_11.json'])
>>> json_file = /home/people/s220672/functional_testing/functional_testing_data/ecoli/test_isolate_11/test_isolate_11_-c.json
>>> with open(json_file, 'r') as f:
...		json_data_true = json.load(f)
>>> json_file_new = /home/people/s220672/functional_testing/functional_testing_temp/test_isolate_11.json
>>> with open(json_file_new, 'r') as f:
...		json_data_new = json.load(f)

```

### Test keys and base structure of json dic

```
>>> assert json_file != json_file_new, 'Json files are the same'
>>> keys_true = list(json_data_true.keys())
>>> keys_new = list(json_data_new.keys())
>>> assert keys_true == keys_new, 'Keys are not equal'

```

### Test seq_region for ecoli

```
>>> seq_region = json_data_true['seq_region']
>>> seq_regions_table_true = pd.DataFrame(seq_region)
>>> seq_region_new = json_data_new['seq_region']
>>> seq_regions_table_new = pd.DataFrame(seq_region_new)
>>> assert seq_regions_table_true.loc['phenotypes'].values.tolist() == seq_regions_table_new['phenotypes'].values.tolist(), 'Phenotypes are not equal'
>>> assert seq_regions_table_true.loc['alignment_length'].values.tolist() == seq_regions_table_new['alignment_length'].values.tolist(), 'Alignment length for {self.species} is not equal'
>>> assert seq_regions_table_true.loc['query_start_pos'].values.tolist() == seq_regions_table_new['query_start_pos'].values.tolist(), 'Query start for {self.species} is not equal'
>>> assert seq_regions_table_true.loc['query_end_pos'].values.tolist() == seq_regions_table_new['query_end_pos'].values.tolist(), 'Query end for {self.species} is not equal'
>>> assert seq_regions_table_true.loc['ref_start_pos'].values.tolist() == seq_regions_table_new['ref_start_pos'].values.tolist(), 'Ref start for {self.species} is not equal'
>>> assert seq_regions_table_true.loc['ref_end_pos'].values.tolist() == seq_regions_table_new['ref_end_pos'].values.tolist(), 'Ref end for {self.species} is not equal'

```

# Test mutations for ecoli and mut_type: chromosomal

```
>>> seq_variations_true = json_data_true['seq_variations']
>>> seq_variations_table_true = pd.DataFrame(seq_variations_true)
>>> seq_variations_new = json_data_new['seq_variations']
>>> seq_variations_table_new = pd.DataFrame(seq_variations_new)
>>> assert seq_variations_table_new.loc['key'].values.tolist()[0].split(';')[0] == seq_variations_table_true.loc['key'].values.tolist()[0].split(';')[0], 'First element of key is not equal'
>>> assert seq_variations_table_new.loc['ref_start_pos'].values.tolist() == seq_variations_table_true.loc['ref_start_pos'].values.tolist(), 'Mutation start position for chromosomal and ecoli is not equal'
>>> assert seq_variations_table_new.loc['ref_end_pos'].values.tolist() == seq_variations_table_true.loc['ref_end_pos'].values.tolist(), 'Mutation end position for chromosomal and ecoli is not equal'
>>> assert seq_variations_table_new.loc['var_codon'].values.tolist() == seq_variations_table_true.loc['var_codon'].values.tolist(), 'Mutation codon for chromosomal and ecoli is not equal'
>>> assert seq_variations_table_new.loc['ref_codon'].values.tolist() == seq_variations_table_true.loc['ref_codon'].values.tolist(), 'Reference codon for chromosomal and ecoli is not equal'

```

## Software Command for senterica

```
>>> subprocess.run(['python', '/home/people/s220672/functional_testing/env_ft/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-ifa', '/home/people/s220672/resfinder/tests/data/test_isolate_11.fa', '-s', 'senterica', '-c', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_11.json'])
>>> json_file = /home/people/s220672/functional_testing/functional_testing_data/senterica/test_isolate_11/test_isolate_11_-c.json
>>> with open(json_file, 'r') as f:
...		json_data_true = json.load(f)
>>> json_file_new = /home/people/s220672/functional_testing/functional_testing_temp/test_isolate_11.json
>>> with open(json_file_new, 'r') as f:
...		json_data_new = json.load(f)

```

### Test keys and base structure of json dic

```
>>> assert json_file != json_file_new, 'Json files are the same'
>>> keys_true = list(json_data_true.keys())
>>> keys_new = list(json_data_new.keys())
>>> assert keys_true == keys_new, 'Keys are not equal'

```

### Test seq_region for senterica

```
>>> seq_region = json_data_true['seq_region']
>>> seq_regions_table_true = pd.DataFrame(seq_region)
>>> seq_region_new = json_data_new['seq_region']
>>> seq_regions_table_new = pd.DataFrame(seq_region_new)
>>> assert seq_regions_table_true.loc['phenotypes'].values.tolist() == seq_regions_table_new['phenotypes'].values.tolist(), 'Phenotypes are not equal'
>>> assert seq_regions_table_true.loc['alignment_length'].values.tolist() == seq_regions_table_new['alignment_length'].values.tolist(), 'Alignment length for {self.species} is not equal'
>>> assert seq_regions_table_true.loc['query_start_pos'].values.tolist() == seq_regions_table_new['query_start_pos'].values.tolist(), 'Query start for {self.species} is not equal'
>>> assert seq_regions_table_true.loc['query_end_pos'].values.tolist() == seq_regions_table_new['query_end_pos'].values.tolist(), 'Query end for {self.species} is not equal'
>>> assert seq_regions_table_true.loc['ref_start_pos'].values.tolist() == seq_regions_table_new['ref_start_pos'].values.tolist(), 'Ref start for {self.species} is not equal'
>>> assert seq_regions_table_true.loc['ref_end_pos'].values.tolist() == seq_regions_table_new['ref_end_pos'].values.tolist(), 'Ref end for {self.species} is not equal'

```

# Test mutations for senterica and mut_type: chromosomal

```
>>> seq_variations_true = json_data_true['seq_variations']
>>> seq_variations_table_true = pd.DataFrame(seq_variations_true)
>>> seq_variations_new = json_data_new['seq_variations']
>>> seq_variations_table_new = pd.DataFrame(seq_variations_new)
>>> assert seq_variations_table_new.loc['key'].values.tolist()[0].split(';')[0] == seq_variations_table_true.loc['key'].values.tolist()[0].split(';')[0], 'First element of key is not equal'
>>> assert seq_variations_table_new.loc['ref_start_pos'].values.tolist() == seq_variations_table_true.loc['ref_start_pos'].values.tolist(), 'Mutation start position for chromosomal and senterica is not equal'
>>> assert seq_variations_table_new.loc['ref_end_pos'].values.tolist() == seq_variations_table_true.loc['ref_end_pos'].values.tolist(), 'Mutation end position for chromosomal and senterica is not equal'
>>> assert seq_variations_table_new.loc['var_codon'].values.tolist() == seq_variations_table_true.loc['var_codon'].values.tolist(), 'Mutation codon for chromosomal and senterica is not equal'
>>> assert seq_variations_table_new.loc['ref_codon'].values.tolist() == seq_variations_table_true.loc['ref_codon'].values.tolist(), 'Reference codon for chromosomal and senterica is not equal'

```

