# Imports

```
>>> import subprocess
>>> import json
>>> import os
>>> import pandas as pd

```

### Fasta Modifier Import

```
>>> from Bio import SeqIO
>>> from Bio.SeqRecord import SeqRecord
>>> from Bio.Seq import Seq
>>> class FastaModifier:
...		@staticmethod
...		def read_fasta(file_path):
...			with open(file_path, 'r') as fasta_file:
...				for record in SeqIO.parse(fasta_file, 'fasta'):
...					header = '>' + record.id
...					sequence = str(record.seq)
...			return sequence, header
...		@staticmethod
...		def create_fasta(output_file,sequences, headers):
...			records = []
...			for header, sequence in zip(headers, sequences):
...				seq_record = SeqRecord(Seq(sequence), id=header[1:], description='')
...				records.append(seq_record)
...			with open(output_file + '.fa', 'w') as fasta_file:
...				SeqIO.write(records, fasta_file, 'fasta')
...		@staticmethod
...		def calc_positions(json_data_true):
...			seq_variations_true = json_data_true['seq_regions']
...			seq_variations_true = seq_variations_true[next(iter(seq_variations_true))]
...			start_pos = seq_variations_true['ref_start_pos']
...			end_pos = seq_variations_true['ref_end_pos']
...			half_pos = int((int(end_pos) - int(start_pos))/2)
...			start_pos = int(start_pos)
...			end_pos = int(end_pos)
...			return start_pos, end_pos, half_pos
...		def modify_seqs(self, json_data_true, fasta_file):
...			start_pos, end_pos, half_pos = self.calc_positions(json_data_true)
...			sequence, header = self.read_fasta(fasta_file)
...			first_contig = sequence[:half_pos]
...			second_contig = sequence[half_pos:end_pos]
...			header2 = header + '_split2'
...			both_seqs = [first_contig, second_contig]
...			both_headers = [header, header2]
...			return both_seqs, both_headers

```

## Software Command for all

```
>>> subprocess.run(['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'ecoli', '-c', '-ifq', '/home/people/s220672/resfinder/tests/data/test_isolate_09b_1.fq', '/home/people/s220672/resfinder/tests/data/test_isolate_09b_2.fq', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_09b_1.json'])
CompletedProcess(args=['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'ecoli', '-c', '-ifq', '/home/people/s220672/resfinder/tests/data/test_isolate_09b_1.fq', '/home/people/s220672/resfinder/tests/data/test_isolate_09b_2.fq', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_09b_1.json'], returncode=0)
>>> json_file = '/home/people/s220672/functional_testing/functional_testing_data/ecoli/test_isolate_09b_1/ofEqB/test_isolate_09b_1_-c.json'
>>> with open(json_file, 'r') as f:
...		json_data_true = json.load(f)
>>> json_file_new = '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_09b_1.json'
>>> with open(json_file_new, 'r') as f:
...		json_data_new = json.load(f)
>>> if not os.path.isdir('test_temp'):
...		os.mkdir('test_temp')

```

### Test keys and base structure of json dic

```
>>> assert json_file != json_file_new, 'Json files are the same'
>>> keys_true = list(json_data_true.keys())
>>> keys_new = list(json_data_new.keys())
>>> assert keys_true == keys_new, 'Keys are not equal'

```

# Test mutations for all and mut_type: ecoli

```
>>> seq_variations_true = json_data_true['seq_variations']
>>> seq_variations_table_true = pd.DataFrame(seq_variations_true)
>>> seq_variations_new = json_data_new['seq_variations']
>>> seq_variations_table_new = pd.DataFrame(seq_variations_new)
>>> assert seq_variations_table_true.shape == seq_variations_table_new.shape, 'Shape of seq_variations is not equal'
>>> if len(seq_variations_table_true) > 0:
...		assert seq_variations_table_new.loc['key'].values.tolist()[0].split(';')[0] == seq_variations_table_true.loc['key'].values.tolist()[0].split(';')[0], 'First element of key is not equal'
...		assert seq_variations_table_new.loc['ref_start_pos'].values.tolist() == seq_variations_table_true.loc['ref_start_pos'].values.tolist(), 'Mutation start position for ecoli and all is not equal'
...		assert seq_variations_table_new.loc['ref_end_pos'].values.tolist() == seq_variations_table_true.loc['ref_end_pos'].values.tolist(), 'Mutation end position for ecoli and all is not equal'
...		assert seq_variations_table_new.loc['var_codon'].values.tolist() == seq_variations_table_true.loc['var_codon'].values.tolist(), 'Mutation codon for ecoli and all is not equal'
...		assert seq_variations_table_new.loc['ref_codon'].values.tolist() == seq_variations_table_true.loc['ref_codon'].values.tolist(), 'Reference codon for ecoli and all is not equal'

```

```
>>> os.remove('/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_09b_1.json')

```

## Software Command for all

```
>>> subprocess.run(['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'ecoli', '-c', '-ifa', '/home/people/s220672/resfinder/tests/data/test_isolate_11.fa', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_11.json'])
CompletedProcess(args=['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'ecoli', '-c', '-ifa', '/home/people/s220672/resfinder/tests/data/test_isolate_11.fa', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_11.json'], returncode=0)
>>> json_file = '/home/people/s220672/functional_testing/functional_testing_data/ecoli/test_isolate_11/uxIxD/test_isolate_11_-c.json'
>>> with open(json_file, 'r') as f:
...		json_data_true = json.load(f)
>>> json_file_new = '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_11.json'
>>> with open(json_file_new, 'r') as f:
...		json_data_new = json.load(f)
>>> if not os.path.isdir('test_temp'):
...		os.mkdir('test_temp')

```

### Test keys and base structure of json dic

```
>>> assert json_file != json_file_new, 'Json files are the same'
>>> keys_true = list(json_data_true.keys())
>>> keys_new = list(json_data_new.keys())
>>> assert keys_true == keys_new, 'Keys are not equal'

```

# Test mutations for all and mut_type: ecoli

```
>>> seq_variations_true = json_data_true['seq_variations']
>>> seq_variations_table_true = pd.DataFrame(seq_variations_true)
>>> seq_variations_new = json_data_new['seq_variations']
>>> seq_variations_table_new = pd.DataFrame(seq_variations_new)
>>> assert seq_variations_table_true.shape == seq_variations_table_new.shape, 'Shape of seq_variations is not equal'
>>> if len(seq_variations_table_true) > 0:
...		assert seq_variations_table_new.loc['key'].values.tolist()[0].split(';')[0] == seq_variations_table_true.loc['key'].values.tolist()[0].split(';')[0], 'First element of key is not equal'
...		assert seq_variations_table_new.loc['ref_start_pos'].values.tolist() == seq_variations_table_true.loc['ref_start_pos'].values.tolist(), 'Mutation start position for ecoli and all is not equal'
...		assert seq_variations_table_new.loc['ref_end_pos'].values.tolist() == seq_variations_table_true.loc['ref_end_pos'].values.tolist(), 'Mutation end position for ecoli and all is not equal'
...		assert seq_variations_table_new.loc['var_codon'].values.tolist() == seq_variations_table_true.loc['var_codon'].values.tolist(), 'Mutation codon for ecoli and all is not equal'
...		assert seq_variations_table_new.loc['ref_codon'].values.tolist() == seq_variations_table_true.loc['ref_codon'].values.tolist(), 'Reference codon for ecoli and all is not equal'

```

```
>>> os.remove('/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_11.json')

```

## Software Command for all

```
>>> subprocess.run(['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'ecoli', '-c', '-ifq', '/home/people/s220672/resfinder/tests/data/test_isolate_05_2.fq', '/home/people/s220672/resfinder/tests/data/test_isolate_05_1.fq', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_05_2.json'])
CompletedProcess(args=['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'ecoli', '-c', '-ifq', '/home/people/s220672/resfinder/tests/data/test_isolate_05_2.fq', '/home/people/s220672/resfinder/tests/data/test_isolate_05_1.fq', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_05_2.json'], returncode=0)
>>> json_file = '/home/people/s220672/functional_testing/functional_testing_data/ecoli/test_isolate_05_2/imxpA/test_isolate_05_2_-c.json'
>>> with open(json_file, 'r') as f:
...		json_data_true = json.load(f)
>>> json_file_new = '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_05_2.json'
>>> with open(json_file_new, 'r') as f:
...		json_data_new = json.load(f)
>>> if not os.path.isdir('test_temp'):
...		os.mkdir('test_temp')

```

### Test keys and base structure of json dic

```
>>> assert json_file != json_file_new, 'Json files are the same'
>>> keys_true = list(json_data_true.keys())
>>> keys_new = list(json_data_new.keys())
>>> assert keys_true == keys_new, 'Keys are not equal'

```

# Test mutations for all and mut_type: ecoli

```
>>> seq_variations_true = json_data_true['seq_variations']
>>> seq_variations_table_true = pd.DataFrame(seq_variations_true)
>>> seq_variations_new = json_data_new['seq_variations']
>>> seq_variations_table_new = pd.DataFrame(seq_variations_new)
>>> assert seq_variations_table_true.shape == seq_variations_table_new.shape, 'Shape of seq_variations is not equal'
>>> if len(seq_variations_table_true) > 0:
...		assert seq_variations_table_new.loc['key'].values.tolist()[0].split(';')[0] == seq_variations_table_true.loc['key'].values.tolist()[0].split(';')[0], 'First element of key is not equal'
...		assert seq_variations_table_new.loc['ref_start_pos'].values.tolist() == seq_variations_table_true.loc['ref_start_pos'].values.tolist(), 'Mutation start position for ecoli and all is not equal'
...		assert seq_variations_table_new.loc['ref_end_pos'].values.tolist() == seq_variations_table_true.loc['ref_end_pos'].values.tolist(), 'Mutation end position for ecoli and all is not equal'
...		assert seq_variations_table_new.loc['var_codon'].values.tolist() == seq_variations_table_true.loc['var_codon'].values.tolist(), 'Mutation codon for ecoli and all is not equal'
...		assert seq_variations_table_new.loc['ref_codon'].values.tolist() == seq_variations_table_true.loc['ref_codon'].values.tolist(), 'Reference codon for ecoli and all is not equal'

```

```
>>> os.remove('/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_05_2.json')

```

## Software Command for all

```
>>> subprocess.run(['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'ecoli', '-c', '-ifq', '/home/people/s220672/resfinder/tests/data/test_isolate_01_1.fq', '/home/people/s220672/resfinder/tests/data/test_isolate_01_2.fq', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_01_1.json'])
CompletedProcess(args=['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'ecoli', '-c', '-ifq', '/home/people/s220672/resfinder/tests/data/test_isolate_01_1.fq', '/home/people/s220672/resfinder/tests/data/test_isolate_01_2.fq', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_01_1.json'], returncode=0)
>>> json_file = '/home/people/s220672/functional_testing/functional_testing_data/ecoli/test_isolate_01_1/lYyyN/test_isolate_01_1_-c.json'
>>> with open(json_file, 'r') as f:
...		json_data_true = json.load(f)
>>> json_file_new = '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_01_1.json'
>>> with open(json_file_new, 'r') as f:
...		json_data_new = json.load(f)
>>> if not os.path.isdir('test_temp'):
...		os.mkdir('test_temp')

```

### Test keys and base structure of json dic

```
>>> assert json_file != json_file_new, 'Json files are the same'
>>> keys_true = list(json_data_true.keys())
>>> keys_new = list(json_data_new.keys())
>>> assert keys_true == keys_new, 'Keys are not equal'

```

# Test mutations for all and mut_type: ecoli

```
>>> seq_variations_true = json_data_true['seq_variations']
>>> seq_variations_table_true = pd.DataFrame(seq_variations_true)
>>> seq_variations_new = json_data_new['seq_variations']
>>> seq_variations_table_new = pd.DataFrame(seq_variations_new)
>>> assert seq_variations_table_true.shape == seq_variations_table_new.shape, 'Shape of seq_variations is not equal'
>>> if len(seq_variations_table_true) > 0:
...		assert seq_variations_table_new.loc['key'].values.tolist()[0].split(';')[0] == seq_variations_table_true.loc['key'].values.tolist()[0].split(';')[0], 'First element of key is not equal'
...		assert seq_variations_table_new.loc['ref_start_pos'].values.tolist() == seq_variations_table_true.loc['ref_start_pos'].values.tolist(), 'Mutation start position for ecoli and all is not equal'
...		assert seq_variations_table_new.loc['ref_end_pos'].values.tolist() == seq_variations_table_true.loc['ref_end_pos'].values.tolist(), 'Mutation end position for ecoli and all is not equal'
...		assert seq_variations_table_new.loc['var_codon'].values.tolist() == seq_variations_table_true.loc['var_codon'].values.tolist(), 'Mutation codon for ecoli and all is not equal'
...		assert seq_variations_table_new.loc['ref_codon'].values.tolist() == seq_variations_table_true.loc['ref_codon'].values.tolist(), 'Reference codon for ecoli and all is not equal'

```

```
>>> os.remove('/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_01_1.json')

```

## Software Command for all

```
>>> subprocess.run(['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'ecoli', '-c', '-ifq', '/home/people/s220672/resfinder/tests/data/test_isolate_11_2.fq', '/home/people/s220672/resfinder/tests/data/test_isolate_11_1.fq', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_11_2.json'])
CompletedProcess(args=['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'ecoli', '-c', '-ifq', '/home/people/s220672/resfinder/tests/data/test_isolate_11_2.fq', '/home/people/s220672/resfinder/tests/data/test_isolate_11_1.fq', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_11_2.json'], returncode=0)
>>> json_file = '/home/people/s220672/functional_testing/functional_testing_data/ecoli/test_isolate_11_2/nYJYj/test_isolate_11_2_-c.json'
>>> with open(json_file, 'r') as f:
...		json_data_true = json.load(f)
>>> json_file_new = '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_11_2.json'
>>> with open(json_file_new, 'r') as f:
...		json_data_new = json.load(f)
>>> if not os.path.isdir('test_temp'):
...		os.mkdir('test_temp')

```

### Test keys and base structure of json dic

```
>>> assert json_file != json_file_new, 'Json files are the same'
>>> keys_true = list(json_data_true.keys())
>>> keys_new = list(json_data_new.keys())
>>> assert keys_true == keys_new, 'Keys are not equal'

```

# Test mutations for all and mut_type: ecoli

```
>>> seq_variations_true = json_data_true['seq_variations']
>>> seq_variations_table_true = pd.DataFrame(seq_variations_true)
>>> seq_variations_new = json_data_new['seq_variations']
>>> seq_variations_table_new = pd.DataFrame(seq_variations_new)
>>> assert seq_variations_table_true.shape == seq_variations_table_new.shape, 'Shape of seq_variations is not equal'
>>> if len(seq_variations_table_true) > 0:
...		assert seq_variations_table_new.loc['key'].values.tolist()[0].split(';')[0] == seq_variations_table_true.loc['key'].values.tolist()[0].split(';')[0], 'First element of key is not equal'
...		assert seq_variations_table_new.loc['ref_start_pos'].values.tolist() == seq_variations_table_true.loc['ref_start_pos'].values.tolist(), 'Mutation start position for ecoli and all is not equal'
...		assert seq_variations_table_new.loc['ref_end_pos'].values.tolist() == seq_variations_table_true.loc['ref_end_pos'].values.tolist(), 'Mutation end position for ecoli and all is not equal'
...		assert seq_variations_table_new.loc['var_codon'].values.tolist() == seq_variations_table_true.loc['var_codon'].values.tolist(), 'Mutation codon for ecoli and all is not equal'
...		assert seq_variations_table_new.loc['ref_codon'].values.tolist() == seq_variations_table_true.loc['ref_codon'].values.tolist(), 'Reference codon for ecoli and all is not equal'

```

```
>>> os.remove('/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_11_2.json')

```

## Software Command for all

```
>>> subprocess.run(['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'ecoli', '-c', '-ifa', '/home/people/s220672/resfinder/tests/data/test_isolate_09a.fa', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_09a.json'])
CompletedProcess(args=['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'ecoli', '-c', '-ifa', '/home/people/s220672/resfinder/tests/data/test_isolate_09a.fa', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_09a.json'], returncode=0)
>>> json_file = '/home/people/s220672/functional_testing/functional_testing_data/ecoli/test_isolate_09a/yAcPE/test_isolate_09a_-c.json'
>>> with open(json_file, 'r') as f:
...		json_data_true = json.load(f)
>>> json_file_new = '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_09a.json'
>>> with open(json_file_new, 'r') as f:
...		json_data_new = json.load(f)
>>> if not os.path.isdir('test_temp'):
...		os.mkdir('test_temp')

```

### Test keys and base structure of json dic

```
>>> assert json_file != json_file_new, 'Json files are the same'
>>> keys_true = list(json_data_true.keys())
>>> keys_new = list(json_data_new.keys())
>>> assert keys_true == keys_new, 'Keys are not equal'

```

# Test mutations for all and mut_type: ecoli

```
>>> seq_variations_true = json_data_true['seq_variations']
>>> seq_variations_table_true = pd.DataFrame(seq_variations_true)
>>> seq_variations_new = json_data_new['seq_variations']
>>> seq_variations_table_new = pd.DataFrame(seq_variations_new)
>>> assert seq_variations_table_true.shape == seq_variations_table_new.shape, 'Shape of seq_variations is not equal'
>>> if len(seq_variations_table_true) > 0:
...		assert seq_variations_table_new.loc['key'].values.tolist()[0].split(';')[0] == seq_variations_table_true.loc['key'].values.tolist()[0].split(';')[0], 'First element of key is not equal'
...		assert seq_variations_table_new.loc['ref_start_pos'].values.tolist() == seq_variations_table_true.loc['ref_start_pos'].values.tolist(), 'Mutation start position for ecoli and all is not equal'
...		assert seq_variations_table_new.loc['ref_end_pos'].values.tolist() == seq_variations_table_true.loc['ref_end_pos'].values.tolist(), 'Mutation end position for ecoli and all is not equal'
...		assert seq_variations_table_new.loc['var_codon'].values.tolist() == seq_variations_table_true.loc['var_codon'].values.tolist(), 'Mutation codon for ecoli and all is not equal'
...		assert seq_variations_table_new.loc['ref_codon'].values.tolist() == seq_variations_table_true.loc['ref_codon'].values.tolist(), 'Reference codon for ecoli and all is not equal'

```

```
>>> os.remove('/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_09a.json')

```

## Software Command for all

```
>>> subprocess.run(['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'ecoli', '-c', '-ifa', '/home/people/s220672/resfinder/tests/data/test_isolate_03.fa', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_03.json'])
CompletedProcess(args=['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'ecoli', '-c', '-ifa', '/home/people/s220672/resfinder/tests/data/test_isolate_03.fa', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_03.json'], returncode=0)
>>> json_file = '/home/people/s220672/functional_testing/functional_testing_data/ecoli/test_isolate_03/sBwLu/test_isolate_03_-c.json'
>>> with open(json_file, 'r') as f:
...		json_data_true = json.load(f)
>>> json_file_new = '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_03.json'
>>> with open(json_file_new, 'r') as f:
...		json_data_new = json.load(f)
>>> if not os.path.isdir('test_temp'):
...		os.mkdir('test_temp')

```

### Test keys and base structure of json dic

```
>>> assert json_file != json_file_new, 'Json files are the same'
>>> keys_true = list(json_data_true.keys())
>>> keys_new = list(json_data_new.keys())
>>> assert keys_true == keys_new, 'Keys are not equal'

```

# Test mutations for all and mut_type: ecoli

```
>>> seq_variations_true = json_data_true['seq_variations']
>>> seq_variations_table_true = pd.DataFrame(seq_variations_true)
>>> seq_variations_new = json_data_new['seq_variations']
>>> seq_variations_table_new = pd.DataFrame(seq_variations_new)
>>> assert seq_variations_table_true.shape == seq_variations_table_new.shape, 'Shape of seq_variations is not equal'
>>> if len(seq_variations_table_true) > 0:
...		assert seq_variations_table_new.loc['key'].values.tolist()[0].split(';')[0] == seq_variations_table_true.loc['key'].values.tolist()[0].split(';')[0], 'First element of key is not equal'
...		assert seq_variations_table_new.loc['ref_start_pos'].values.tolist() == seq_variations_table_true.loc['ref_start_pos'].values.tolist(), 'Mutation start position for ecoli and all is not equal'
...		assert seq_variations_table_new.loc['ref_end_pos'].values.tolist() == seq_variations_table_true.loc['ref_end_pos'].values.tolist(), 'Mutation end position for ecoli and all is not equal'
...		assert seq_variations_table_new.loc['var_codon'].values.tolist() == seq_variations_table_true.loc['var_codon'].values.tolist(), 'Mutation codon for ecoli and all is not equal'
...		assert seq_variations_table_new.loc['ref_codon'].values.tolist() == seq_variations_table_true.loc['ref_codon'].values.tolist(), 'Reference codon for ecoli and all is not equal'

```

```
>>> os.remove('/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_03.json')

```

## Software Command for all

```
>>> subprocess.run(['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'ecoli', '-c', '-ifa', '/home/people/s220672/resfinder/tests/data/test_isolate_09b.fa', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_09b.json'])
CompletedProcess(args=['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'ecoli', '-c', '-ifa', '/home/people/s220672/resfinder/tests/data/test_isolate_09b.fa', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_09b.json'], returncode=0)
>>> json_file = '/home/people/s220672/functional_testing/functional_testing_data/ecoli/test_isolate_09b/sdgde/test_isolate_09b_-c.json'
>>> with open(json_file, 'r') as f:
...		json_data_true = json.load(f)
>>> json_file_new = '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_09b.json'
>>> with open(json_file_new, 'r') as f:
...		json_data_new = json.load(f)
>>> if not os.path.isdir('test_temp'):
...		os.mkdir('test_temp')

```

### Test keys and base structure of json dic

```
>>> assert json_file != json_file_new, 'Json files are the same'
>>> keys_true = list(json_data_true.keys())
>>> keys_new = list(json_data_new.keys())
>>> assert keys_true == keys_new, 'Keys are not equal'

```

# Test mutations for all and mut_type: ecoli

```
>>> seq_variations_true = json_data_true['seq_variations']
>>> seq_variations_table_true = pd.DataFrame(seq_variations_true)
>>> seq_variations_new = json_data_new['seq_variations']
>>> seq_variations_table_new = pd.DataFrame(seq_variations_new)
>>> assert seq_variations_table_true.shape == seq_variations_table_new.shape, 'Shape of seq_variations is not equal'
>>> if len(seq_variations_table_true) > 0:
...		assert seq_variations_table_new.loc['key'].values.tolist()[0].split(';')[0] == seq_variations_table_true.loc['key'].values.tolist()[0].split(';')[0], 'First element of key is not equal'
...		assert seq_variations_table_new.loc['ref_start_pos'].values.tolist() == seq_variations_table_true.loc['ref_start_pos'].values.tolist(), 'Mutation start position for ecoli and all is not equal'
...		assert seq_variations_table_new.loc['ref_end_pos'].values.tolist() == seq_variations_table_true.loc['ref_end_pos'].values.tolist(), 'Mutation end position for ecoli and all is not equal'
...		assert seq_variations_table_new.loc['var_codon'].values.tolist() == seq_variations_table_true.loc['var_codon'].values.tolist(), 'Mutation codon for ecoli and all is not equal'
...		assert seq_variations_table_new.loc['ref_codon'].values.tolist() == seq_variations_table_true.loc['ref_codon'].values.tolist(), 'Reference codon for ecoli and all is not equal'

```

```
>>> os.remove('/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_09b.json')

```

## Software Command for all

```
>>> subprocess.run(['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'ecoli', '-c', '-ifa', '/home/people/s220672/resfinder/tests/data/test_isolate_01.fa', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_01.json'])
CompletedProcess(args=['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'ecoli', '-c', '-ifa', '/home/people/s220672/resfinder/tests/data/test_isolate_01.fa', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_01.json'], returncode=0)
>>> json_file = '/home/people/s220672/functional_testing/functional_testing_data/ecoli/test_isolate_01/aqfdf/test_isolate_01_-c.json'
>>> with open(json_file, 'r') as f:
...		json_data_true = json.load(f)
>>> json_file_new = '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_01.json'
>>> with open(json_file_new, 'r') as f:
...		json_data_new = json.load(f)
>>> if not os.path.isdir('test_temp'):
...		os.mkdir('test_temp')

```

### Test keys and base structure of json dic

```
>>> assert json_file != json_file_new, 'Json files are the same'
>>> keys_true = list(json_data_true.keys())
>>> keys_new = list(json_data_new.keys())
>>> assert keys_true == keys_new, 'Keys are not equal'

```

# Test mutations for all and mut_type: ecoli

```
>>> seq_variations_true = json_data_true['seq_variations']
>>> seq_variations_table_true = pd.DataFrame(seq_variations_true)
>>> seq_variations_new = json_data_new['seq_variations']
>>> seq_variations_table_new = pd.DataFrame(seq_variations_new)
>>> assert seq_variations_table_true.shape == seq_variations_table_new.shape, 'Shape of seq_variations is not equal'
>>> if len(seq_variations_table_true) > 0:
...		assert seq_variations_table_new.loc['key'].values.tolist()[0].split(';')[0] == seq_variations_table_true.loc['key'].values.tolist()[0].split(';')[0], 'First element of key is not equal'
...		assert seq_variations_table_new.loc['ref_start_pos'].values.tolist() == seq_variations_table_true.loc['ref_start_pos'].values.tolist(), 'Mutation start position for ecoli and all is not equal'
...		assert seq_variations_table_new.loc['ref_end_pos'].values.tolist() == seq_variations_table_true.loc['ref_end_pos'].values.tolist(), 'Mutation end position for ecoli and all is not equal'
...		assert seq_variations_table_new.loc['var_codon'].values.tolist() == seq_variations_table_true.loc['var_codon'].values.tolist(), 'Mutation codon for ecoli and all is not equal'
...		assert seq_variations_table_new.loc['ref_codon'].values.tolist() == seq_variations_table_true.loc['ref_codon'].values.tolist(), 'Reference codon for ecoli and all is not equal'

```

```
>>> os.remove('/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_01.json')

```

## Software Command for all

```
>>> subprocess.run(['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'ecoli', '-c', '-ifa', '/home/people/s220672/resfinder/tests/data/test_isolate_05.fa', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_05.json'])
CompletedProcess(args=['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'ecoli', '-c', '-ifa', '/home/people/s220672/resfinder/tests/data/test_isolate_05.fa', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_05.json'], returncode=0)
>>> json_file = '/home/people/s220672/functional_testing/functional_testing_data/ecoli/test_isolate_05/fuivT/test_isolate_05_-c.json'
>>> with open(json_file, 'r') as f:
...		json_data_true = json.load(f)
>>> json_file_new = '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_05.json'
>>> with open(json_file_new, 'r') as f:
...		json_data_new = json.load(f)
>>> if not os.path.isdir('test_temp'):
...		os.mkdir('test_temp')

```

### Test keys and base structure of json dic

```
>>> assert json_file != json_file_new, 'Json files are the same'
>>> keys_true = list(json_data_true.keys())
>>> keys_new = list(json_data_new.keys())
>>> assert keys_true == keys_new, 'Keys are not equal'

```

# Test mutations for all and mut_type: ecoli

```
>>> seq_variations_true = json_data_true['seq_variations']
>>> seq_variations_table_true = pd.DataFrame(seq_variations_true)
>>> seq_variations_new = json_data_new['seq_variations']
>>> seq_variations_table_new = pd.DataFrame(seq_variations_new)
>>> assert seq_variations_table_true.shape == seq_variations_table_new.shape, 'Shape of seq_variations is not equal'
>>> if len(seq_variations_table_true) > 0:
...		assert seq_variations_table_new.loc['key'].values.tolist()[0].split(';')[0] == seq_variations_table_true.loc['key'].values.tolist()[0].split(';')[0], 'First element of key is not equal'
...		assert seq_variations_table_new.loc['ref_start_pos'].values.tolist() == seq_variations_table_true.loc['ref_start_pos'].values.tolist(), 'Mutation start position for ecoli and all is not equal'
...		assert seq_variations_table_new.loc['ref_end_pos'].values.tolist() == seq_variations_table_true.loc['ref_end_pos'].values.tolist(), 'Mutation end position for ecoli and all is not equal'
...		assert seq_variations_table_new.loc['var_codon'].values.tolist() == seq_variations_table_true.loc['var_codon'].values.tolist(), 'Mutation codon for ecoli and all is not equal'
...		assert seq_variations_table_new.loc['ref_codon'].values.tolist() == seq_variations_table_true.loc['ref_codon'].values.tolist(), 'Reference codon for ecoli and all is not equal'

```

```
>>> os.remove('/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_05.json')

```

## Software Command for all

```
>>> subprocess.run(['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'ecoli', '-c', '-ifa', '/home/people/s220672/resfinder/tests/data/test_isolate_10.fa', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_10.json'])
CompletedProcess(args=['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'ecoli', '-c', '-ifa', '/home/people/s220672/resfinder/tests/data/test_isolate_10.fa', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_10.json'], returncode=0)
>>> json_file = '/home/people/s220672/functional_testing/functional_testing_data/ecoli/test_isolate_10/wJhTa/test_isolate_10_-c.json'
>>> with open(json_file, 'r') as f:
...		json_data_true = json.load(f)
>>> json_file_new = '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_10.json'
>>> with open(json_file_new, 'r') as f:
...		json_data_new = json.load(f)
>>> if not os.path.isdir('test_temp'):
...		os.mkdir('test_temp')

```

### Test keys and base structure of json dic

```
>>> assert json_file != json_file_new, 'Json files are the same'
>>> keys_true = list(json_data_true.keys())
>>> keys_new = list(json_data_new.keys())
>>> assert keys_true == keys_new, 'Keys are not equal'

```

# Test mutations for all and mut_type: ecoli

```
>>> seq_variations_true = json_data_true['seq_variations']
>>> seq_variations_table_true = pd.DataFrame(seq_variations_true)
>>> seq_variations_new = json_data_new['seq_variations']
>>> seq_variations_table_new = pd.DataFrame(seq_variations_new)
>>> assert seq_variations_table_true.shape == seq_variations_table_new.shape, 'Shape of seq_variations is not equal'
>>> if len(seq_variations_table_true) > 0:
...		assert seq_variations_table_new.loc['key'].values.tolist()[0].split(';')[0] == seq_variations_table_true.loc['key'].values.tolist()[0].split(';')[0], 'First element of key is not equal'
...		assert seq_variations_table_new.loc['ref_start_pos'].values.tolist() == seq_variations_table_true.loc['ref_start_pos'].values.tolist(), 'Mutation start position for ecoli and all is not equal'
...		assert seq_variations_table_new.loc['ref_end_pos'].values.tolist() == seq_variations_table_true.loc['ref_end_pos'].values.tolist(), 'Mutation end position for ecoli and all is not equal'
...		assert seq_variations_table_new.loc['var_codon'].values.tolist() == seq_variations_table_true.loc['var_codon'].values.tolist(), 'Mutation codon for ecoli and all is not equal'
...		assert seq_variations_table_new.loc['ref_codon'].values.tolist() == seq_variations_table_true.loc['ref_codon'].values.tolist(), 'Reference codon for ecoli and all is not equal'

```

```
>>> os.remove('/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_10.json')

```

