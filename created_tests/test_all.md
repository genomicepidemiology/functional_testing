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

## Software Command for ecoli

```
>>> subprocess.run(['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'ecoli', '-c', '-ifa', '/home/people/s220672/resfinder/tests/data/test_isolate_09a.fa', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_09a.json'])
CompletedProcess(args=['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'ecoli', '-c', '-ifa', '/home/people/s220672/resfinder/tests/data/test_isolate_09a.fa', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_09a.json'], returncode=0)
>>> json_file = '/home/people/s220672/functional_testing/functional_testing_data/ecoli/test_isolate_09a/vzDbp/test_isolate_09a_-c.json'
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

### Test seq_regions for ecoli

```
>>> seq_regions = json_data_true['seq_regions']
>>> seq_regions_table_true = pd.DataFrame(seq_regions)
>>> seq_region_new = json_data_new['seq_regions']
>>> seq_regions_table_new = pd.DataFrame(seq_region_new)
>>> assert seq_regions_table_true.loc['phenotypes'].values.tolist() == seq_regions_table_new.loc['phenotypes'].values.tolist(), 'Phenotypes are not equal'
>>> assert seq_regions_table_true.loc['alignment_length'].values.tolist() == seq_regions_table_new.loc['alignment_length'].values.tolist(), f'Alignment length for ecoli is not equal'
>>> assert seq_regions_table_true.loc['ref_start_pos'].values.tolist() == seq_regions_table_new.loc['ref_start_pos'].values.tolist(), f'Ref start for ecoli is not equal'
>>> assert seq_regions_table_true.loc['ref_end_pos'].values.tolist() == seq_regions_table_new.loc['ref_end_pos'].values.tolist(), f'Ref end for ecoli is not equal'

```

# Test mutations for ecoli and mut_type: chromosomal

```
>>> seq_variations_true = json_data_true['seq_variations']
>>> seq_variations_table_true = pd.DataFrame(seq_variations_true)
>>> seq_variations_new = json_data_new['seq_variations']
>>> seq_variations_table_new = pd.DataFrame(seq_variations_new)
>>> assert seq_variations_table_true.shape == seq_variations_table_new.shape, 'Shape of seq_variations is not equal'
>>> if len(seq_variations_table_true) > 0:
...		assert seq_variations_table_new.loc['key'].values.tolist()[0].split(';')[0] == seq_variations_table_true.loc['key'].values.tolist()[0].split(';')[0], 'First element of key is not equal'
...		assert seq_variations_table_new.loc['ref_start_pos'].values.tolist() == seq_variations_table_true.loc['ref_start_pos'].values.tolist(), 'Mutation start position for chromosomal and ecoli is not equal'
...		assert seq_variations_table_new.loc['ref_end_pos'].values.tolist() == seq_variations_table_true.loc['ref_end_pos'].values.tolist(), 'Mutation end position for chromosomal and ecoli is not equal'
...		assert seq_variations_table_new.loc['var_codon'].values.tolist() == seq_variations_table_true.loc['var_codon'].values.tolist(), 'Mutation codon for chromosomal and ecoli is not equal'
...		assert seq_variations_table_new.loc['ref_codon'].values.tolist() == seq_variations_table_true.loc['ref_codon'].values.tolist(), 'Reference codon for chromosomal and ecoli is not equal'

```

### Test if resistant gene is same after splitting in two contigs.

```
>>> with open('/home/people/s220672/resfinder/tests/data/test_isolate_09a.fa', 'r') as f:
...		fasta = f.read()
>>> ModFasta = FastaModifier()
>>> both_seqs, both_headers = ModFasta.modify_seqs(json_data_true, '/home/people/s220672/resfinder/tests/data/test_isolate_09a.fa')
>>> ModFasta.create_fasta('test_temp/new_fasta', both_seqs, both_headers)
>>> subprocess.run(['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'ecoli', '-c', '-ifa', 'test_temp/new_fasta.fa', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_09a.json'])
CompletedProcess(args=['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'ecoli', '-c', '-ifa', 'test_temp/new_fasta.fa', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_09a.json'], returncode=0)
>>> json_file_new = '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_09a.json'
>>> with open(json_file_new, 'r') as f:
...		json_data_new = json.load(f)
>>> seq_variations = json_data_true['seq_variations']
>>> seq_variations_table_true = pd.DataFrame(seq_variations)
>>> if len(seq_variations_table_true) > 0:
...		seq_variations_new = json_data_new['seq_variations']
...		seq_variations_table_new = pd.DataFrame(seq_variations_new)
...		assert seq_variations_table_true.loc['seq_var'].values.tolist() == seq_variations_table_new.loc['seq_var'].values.tolist(), f'Reference sequence length for ecoli is not equal'
...		assert seq_variations_table_true.loc['ref_end_pos'].values.tolist() == seq_variations_table_new.loc['ref_end_pos'].values.tolist(), f'Ref id for ecoli is not equal'

```

```
>>> os.remove('/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_09a.json')

```

## Software Command for senterica

```
>>> subprocess.run(['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'senterica', '-c', '-ifa', '/home/people/s220672/resfinder/tests/data/test_isolate_09a.fa', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_09a.json'])
CompletedProcess(args=['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'senterica', '-c', '-ifa', '/home/people/s220672/resfinder/tests/data/test_isolate_09a.fa', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_09a.json'], returncode=0)
>>> json_file = '/home/people/s220672/functional_testing/functional_testing_data/senterica/test_isolate_09a/EkOHf/test_isolate_09a_-c.json'
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

### Test seq_regions for senterica

```
>>> seq_regions = json_data_true['seq_regions']
>>> seq_regions_table_true = pd.DataFrame(seq_regions)
>>> seq_region_new = json_data_new['seq_regions']
>>> seq_regions_table_new = pd.DataFrame(seq_region_new)
>>> assert seq_regions_table_true.loc['phenotypes'].values.tolist() == seq_regions_table_new.loc['phenotypes'].values.tolist(), 'Phenotypes are not equal'
>>> assert seq_regions_table_true.loc['alignment_length'].values.tolist() == seq_regions_table_new.loc['alignment_length'].values.tolist(), f'Alignment length for senterica is not equal'
>>> assert seq_regions_table_true.loc['ref_start_pos'].values.tolist() == seq_regions_table_new.loc['ref_start_pos'].values.tolist(), f'Ref start for senterica is not equal'
>>> assert seq_regions_table_true.loc['ref_end_pos'].values.tolist() == seq_regions_table_new.loc['ref_end_pos'].values.tolist(), f'Ref end for senterica is not equal'

```

# Test mutations for senterica and mut_type: chromosomal

```
>>> seq_variations_true = json_data_true['seq_variations']
>>> seq_variations_table_true = pd.DataFrame(seq_variations_true)
>>> seq_variations_new = json_data_new['seq_variations']
>>> seq_variations_table_new = pd.DataFrame(seq_variations_new)
>>> assert seq_variations_table_true.shape == seq_variations_table_new.shape, 'Shape of seq_variations is not equal'
>>> if len(seq_variations_table_true) > 0:
...		assert seq_variations_table_new.loc['key'].values.tolist()[0].split(';')[0] == seq_variations_table_true.loc['key'].values.tolist()[0].split(';')[0], 'First element of key is not equal'
...		assert seq_variations_table_new.loc['ref_start_pos'].values.tolist() == seq_variations_table_true.loc['ref_start_pos'].values.tolist(), 'Mutation start position for chromosomal and senterica is not equal'
...		assert seq_variations_table_new.loc['ref_end_pos'].values.tolist() == seq_variations_table_true.loc['ref_end_pos'].values.tolist(), 'Mutation end position for chromosomal and senterica is not equal'
...		assert seq_variations_table_new.loc['var_codon'].values.tolist() == seq_variations_table_true.loc['var_codon'].values.tolist(), 'Mutation codon for chromosomal and senterica is not equal'
...		assert seq_variations_table_new.loc['ref_codon'].values.tolist() == seq_variations_table_true.loc['ref_codon'].values.tolist(), 'Reference codon for chromosomal and senterica is not equal'

```

```
>>> os.remove('/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_09a.json')

```

## Software Command for ecoli

```
>>> subprocess.run(['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'ecoli', '-c', '-ifq', '/home/people/s220672/resfinder/tests/data/test_isolate_09b_2.fq', '/home/people/s220672/resfinder/tests/data/test_isolate_09b_1.fq', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_09b_2.json'])
CompletedProcess(args=['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'ecoli', '-c', '-ifq', '/home/people/s220672/resfinder/tests/data/test_isolate_09b_2.fq', '/home/people/s220672/resfinder/tests/data/test_isolate_09b_1.fq', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_09b_2.json'], returncode=0)
>>> json_file = '/home/people/s220672/functional_testing/functional_testing_data/ecoli/test_isolate_09b_2/lwIgW/test_isolate_09b_2_-c.json'
>>> with open(json_file, 'r') as f:
...		json_data_true = json.load(f)
>>> json_file_new = '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_09b_2.json'
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

### Test seq_regions for ecoli

```
>>> seq_regions = json_data_true['seq_regions']
>>> seq_regions_table_true = pd.DataFrame(seq_regions)
>>> seq_region_new = json_data_new['seq_regions']
>>> seq_regions_table_new = pd.DataFrame(seq_region_new)
>>> assert seq_regions_table_true.loc['phenotypes'].values.tolist() == seq_regions_table_new.loc['phenotypes'].values.tolist(), 'Phenotypes are not equal'
>>> assert seq_regions_table_true.loc['alignment_length'].values.tolist() == seq_regions_table_new.loc['alignment_length'].values.tolist(), f'Alignment length for ecoli is not equal'
>>> assert seq_regions_table_true.loc['ref_start_pos'].values.tolist() == seq_regions_table_new.loc['ref_start_pos'].values.tolist(), f'Ref start for ecoli is not equal'
>>> assert seq_regions_table_true.loc['ref_end_pos'].values.tolist() == seq_regions_table_new.loc['ref_end_pos'].values.tolist(), f'Ref end for ecoli is not equal'

```

# Test mutations for ecoli and mut_type: chromosomal

```
>>> seq_variations_true = json_data_true['seq_variations']
>>> seq_variations_table_true = pd.DataFrame(seq_variations_true)
>>> seq_variations_new = json_data_new['seq_variations']
>>> seq_variations_table_new = pd.DataFrame(seq_variations_new)
>>> assert seq_variations_table_true.shape == seq_variations_table_new.shape, 'Shape of seq_variations is not equal'
>>> if len(seq_variations_table_true) > 0:
...		assert seq_variations_table_new.loc['key'].values.tolist()[0].split(';')[0] == seq_variations_table_true.loc['key'].values.tolist()[0].split(';')[0], 'First element of key is not equal'
...		assert seq_variations_table_new.loc['ref_start_pos'].values.tolist() == seq_variations_table_true.loc['ref_start_pos'].values.tolist(), 'Mutation start position for chromosomal and ecoli is not equal'
...		assert seq_variations_table_new.loc['ref_end_pos'].values.tolist() == seq_variations_table_true.loc['ref_end_pos'].values.tolist(), 'Mutation end position for chromosomal and ecoli is not equal'
...		assert seq_variations_table_new.loc['var_codon'].values.tolist() == seq_variations_table_true.loc['var_codon'].values.tolist(), 'Mutation codon for chromosomal and ecoli is not equal'
...		assert seq_variations_table_new.loc['ref_codon'].values.tolist() == seq_variations_table_true.loc['ref_codon'].values.tolist(), 'Reference codon for chromosomal and ecoli is not equal'

```

```
>>> os.remove('/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_09b_2.json')

```

## Software Command for ecoli

```
>>> subprocess.run(['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'ecoli', '-c', '-ifq', '/home/people/s220672/resfinder/tests/data/test_isolate_11_2.fq', '/home/people/s220672/resfinder/tests/data/test_isolate_11_1.fq', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_11_2.json'])
CompletedProcess(args=['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'ecoli', '-c', '-ifq', '/home/people/s220672/resfinder/tests/data/test_isolate_11_2.fq', '/home/people/s220672/resfinder/tests/data/test_isolate_11_1.fq', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_11_2.json'], returncode=0)
>>> json_file = '/home/people/s220672/functional_testing/functional_testing_data/ecoli/test_isolate_11_2/KPfsc/test_isolate_11_2_-c.json'
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

### Test seq_regions for ecoli

```
>>> seq_regions = json_data_true['seq_regions']
>>> seq_regions_table_true = pd.DataFrame(seq_regions)
>>> seq_region_new = json_data_new['seq_regions']
>>> seq_regions_table_new = pd.DataFrame(seq_region_new)
>>> assert seq_regions_table_true.loc['phenotypes'].values.tolist() == seq_regions_table_new.loc['phenotypes'].values.tolist(), 'Phenotypes are not equal'
>>> assert seq_regions_table_true.loc['alignment_length'].values.tolist() == seq_regions_table_new.loc['alignment_length'].values.tolist(), f'Alignment length for ecoli is not equal'
>>> assert seq_regions_table_true.loc['ref_start_pos'].values.tolist() == seq_regions_table_new.loc['ref_start_pos'].values.tolist(), f'Ref start for ecoli is not equal'
>>> assert seq_regions_table_true.loc['ref_end_pos'].values.tolist() == seq_regions_table_new.loc['ref_end_pos'].values.tolist(), f'Ref end for ecoli is not equal'

```

# Test mutations for ecoli and mut_type: chromosomal

```
>>> seq_variations_true = json_data_true['seq_variations']
>>> seq_variations_table_true = pd.DataFrame(seq_variations_true)
>>> seq_variations_new = json_data_new['seq_variations']
>>> seq_variations_table_new = pd.DataFrame(seq_variations_new)
>>> assert seq_variations_table_true.shape == seq_variations_table_new.shape, 'Shape of seq_variations is not equal'
>>> if len(seq_variations_table_true) > 0:
...		assert seq_variations_table_new.loc['key'].values.tolist()[0].split(';')[0] == seq_variations_table_true.loc['key'].values.tolist()[0].split(';')[0], 'First element of key is not equal'
...		assert seq_variations_table_new.loc['ref_start_pos'].values.tolist() == seq_variations_table_true.loc['ref_start_pos'].values.tolist(), 'Mutation start position for chromosomal and ecoli is not equal'
...		assert seq_variations_table_new.loc['ref_end_pos'].values.tolist() == seq_variations_table_true.loc['ref_end_pos'].values.tolist(), 'Mutation end position for chromosomal and ecoli is not equal'
...		assert seq_variations_table_new.loc['var_codon'].values.tolist() == seq_variations_table_true.loc['var_codon'].values.tolist(), 'Mutation codon for chromosomal and ecoli is not equal'
...		assert seq_variations_table_new.loc['ref_codon'].values.tolist() == seq_variations_table_true.loc['ref_codon'].values.tolist(), 'Reference codon for chromosomal and ecoli is not equal'

```

```
>>> os.remove('/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_11_2.json')

```

## Software Command for senterica

```
>>> subprocess.run(['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'senterica', '-c', '-ifq', '/home/people/s220672/resfinder/tests/data/test_isolate_11_2.fq', '/home/people/s220672/resfinder/tests/data/test_isolate_11_1.fq', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_11_2.json'])
CompletedProcess(args=['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'senterica', '-c', '-ifq', '/home/people/s220672/resfinder/tests/data/test_isolate_11_2.fq', '/home/people/s220672/resfinder/tests/data/test_isolate_11_1.fq', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_11_2.json'], returncode=0)
>>> json_file = '/home/people/s220672/functional_testing/functional_testing_data/senterica/test_isolate_11_2/SItDr/test_isolate_11_2_-c.json'
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

### Test seq_regions for senterica

```
>>> seq_regions = json_data_true['seq_regions']
>>> seq_regions_table_true = pd.DataFrame(seq_regions)
>>> seq_region_new = json_data_new['seq_regions']
>>> seq_regions_table_new = pd.DataFrame(seq_region_new)
>>> assert seq_regions_table_true.loc['phenotypes'].values.tolist() == seq_regions_table_new.loc['phenotypes'].values.tolist(), 'Phenotypes are not equal'
>>> assert seq_regions_table_true.loc['alignment_length'].values.tolist() == seq_regions_table_new.loc['alignment_length'].values.tolist(), f'Alignment length for senterica is not equal'
>>> assert seq_regions_table_true.loc['ref_start_pos'].values.tolist() == seq_regions_table_new.loc['ref_start_pos'].values.tolist(), f'Ref start for senterica is not equal'
>>> assert seq_regions_table_true.loc['ref_end_pos'].values.tolist() == seq_regions_table_new.loc['ref_end_pos'].values.tolist(), f'Ref end for senterica is not equal'

```

# Test mutations for senterica and mut_type: chromosomal

```
>>> seq_variations_true = json_data_true['seq_variations']
>>> seq_variations_table_true = pd.DataFrame(seq_variations_true)
>>> seq_variations_new = json_data_new['seq_variations']
>>> seq_variations_table_new = pd.DataFrame(seq_variations_new)
>>> assert seq_variations_table_true.shape == seq_variations_table_new.shape, 'Shape of seq_variations is not equal'
>>> if len(seq_variations_table_true) > 0:
...		assert seq_variations_table_new.loc['key'].values.tolist()[0].split(';')[0] == seq_variations_table_true.loc['key'].values.tolist()[0].split(';')[0], 'First element of key is not equal'
...		assert seq_variations_table_new.loc['ref_start_pos'].values.tolist() == seq_variations_table_true.loc['ref_start_pos'].values.tolist(), 'Mutation start position for chromosomal and senterica is not equal'
...		assert seq_variations_table_new.loc['ref_end_pos'].values.tolist() == seq_variations_table_true.loc['ref_end_pos'].values.tolist(), 'Mutation end position for chromosomal and senterica is not equal'
...		assert seq_variations_table_new.loc['var_codon'].values.tolist() == seq_variations_table_true.loc['var_codon'].values.tolist(), 'Mutation codon for chromosomal and senterica is not equal'
...		assert seq_variations_table_new.loc['ref_codon'].values.tolist() == seq_variations_table_true.loc['ref_codon'].values.tolist(), 'Reference codon for chromosomal and senterica is not equal'

```

```
>>> os.remove('/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_11_2.json')

```

## Software Command for ecoli

```
>>> subprocess.run(['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'ecoli', '-c', '-ifa', '/home/people/s220672/resfinder/tests/data/test_isolate_10.fa', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_10.json'])
CompletedProcess(args=['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'ecoli', '-c', '-ifa', '/home/people/s220672/resfinder/tests/data/test_isolate_10.fa', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_10.json'], returncode=0)
>>> json_file = '/home/people/s220672/functional_testing/functional_testing_data/ecoli/test_isolate_10/Rzjre/test_isolate_10_-c.json'
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

### Test seq_regions for ecoli

```
>>> seq_regions = json_data_true['seq_regions']
>>> seq_regions_table_true = pd.DataFrame(seq_regions)
>>> seq_region_new = json_data_new['seq_regions']
>>> seq_regions_table_new = pd.DataFrame(seq_region_new)
>>> assert seq_regions_table_true.loc['phenotypes'].values.tolist() == seq_regions_table_new.loc['phenotypes'].values.tolist(), 'Phenotypes are not equal'
>>> assert seq_regions_table_true.loc['alignment_length'].values.tolist() == seq_regions_table_new.loc['alignment_length'].values.tolist(), f'Alignment length for ecoli is not equal'
>>> assert seq_regions_table_true.loc['ref_start_pos'].values.tolist() == seq_regions_table_new.loc['ref_start_pos'].values.tolist(), f'Ref start for ecoli is not equal'
>>> assert seq_regions_table_true.loc['ref_end_pos'].values.tolist() == seq_regions_table_new.loc['ref_end_pos'].values.tolist(), f'Ref end for ecoli is not equal'

```

# Test mutations for ecoli and mut_type: chromosomal

```
>>> seq_variations_true = json_data_true['seq_variations']
>>> seq_variations_table_true = pd.DataFrame(seq_variations_true)
>>> seq_variations_new = json_data_new['seq_variations']
>>> seq_variations_table_new = pd.DataFrame(seq_variations_new)
>>> assert seq_variations_table_true.shape == seq_variations_table_new.shape, 'Shape of seq_variations is not equal'
>>> if len(seq_variations_table_true) > 0:
...		assert seq_variations_table_new.loc['key'].values.tolist()[0].split(';')[0] == seq_variations_table_true.loc['key'].values.tolist()[0].split(';')[0], 'First element of key is not equal'
...		assert seq_variations_table_new.loc['ref_start_pos'].values.tolist() == seq_variations_table_true.loc['ref_start_pos'].values.tolist(), 'Mutation start position for chromosomal and ecoli is not equal'
...		assert seq_variations_table_new.loc['ref_end_pos'].values.tolist() == seq_variations_table_true.loc['ref_end_pos'].values.tolist(), 'Mutation end position for chromosomal and ecoli is not equal'
...		assert seq_variations_table_new.loc['var_codon'].values.tolist() == seq_variations_table_true.loc['var_codon'].values.tolist(), 'Mutation codon for chromosomal and ecoli is not equal'
...		assert seq_variations_table_new.loc['ref_codon'].values.tolist() == seq_variations_table_true.loc['ref_codon'].values.tolist(), 'Reference codon for chromosomal and ecoli is not equal'

```

### Test if resistant gene is same after splitting in two contigs.

```
>>> with open('/home/people/s220672/resfinder/tests/data/test_isolate_10.fa', 'r') as f:
...		fasta = f.read()
>>> ModFasta = FastaModifier()
>>> both_seqs, both_headers = ModFasta.modify_seqs(json_data_true, '/home/people/s220672/resfinder/tests/data/test_isolate_10.fa')
>>> ModFasta.create_fasta('test_temp/new_fasta', both_seqs, both_headers)
>>> subprocess.run(['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'ecoli', '-c', '-ifa', 'test_temp/new_fasta.fa', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_10.json'])
CompletedProcess(args=['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'ecoli', '-c', '-ifa', 'test_temp/new_fasta.fa', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_10.json'], returncode=0)
>>> json_file_new = '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_10.json'
>>> with open(json_file_new, 'r') as f:
...		json_data_new = json.load(f)
>>> seq_variations = json_data_true['seq_variations']
>>> seq_variations_table_true = pd.DataFrame(seq_variations)
>>> if len(seq_variations_table_true) > 0:
...		seq_variations_new = json_data_new['seq_variations']
...		seq_variations_table_new = pd.DataFrame(seq_variations_new)
...		assert seq_variations_table_true.loc['seq_var'].values.tolist() == seq_variations_table_new.loc['seq_var'].values.tolist(), f'Reference sequence length for ecoli is not equal'
...		assert seq_variations_table_true.loc['ref_end_pos'].values.tolist() == seq_variations_table_new.loc['ref_end_pos'].values.tolist(), f'Ref id for ecoli is not equal'

```

```
>>> os.remove('/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_10.json')

```

## Software Command for senterica

```
>>> subprocess.run(['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'senterica', '-c', '-ifa', '/home/people/s220672/resfinder/tests/data/test_isolate_10.fa', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_10.json'])
CompletedProcess(args=['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'senterica', '-c', '-ifa', '/home/people/s220672/resfinder/tests/data/test_isolate_10.fa', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_10.json'], returncode=0)
>>> json_file = '/home/people/s220672/functional_testing/functional_testing_data/senterica/test_isolate_10/kLmwP/test_isolate_10_-c.json'
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

### Test seq_regions for senterica

```
>>> seq_regions = json_data_true['seq_regions']
>>> seq_regions_table_true = pd.DataFrame(seq_regions)
>>> seq_region_new = json_data_new['seq_regions']
>>> seq_regions_table_new = pd.DataFrame(seq_region_new)
>>> assert seq_regions_table_true.loc['phenotypes'].values.tolist() == seq_regions_table_new.loc['phenotypes'].values.tolist(), 'Phenotypes are not equal'
>>> assert seq_regions_table_true.loc['alignment_length'].values.tolist() == seq_regions_table_new.loc['alignment_length'].values.tolist(), f'Alignment length for senterica is not equal'
>>> assert seq_regions_table_true.loc['ref_start_pos'].values.tolist() == seq_regions_table_new.loc['ref_start_pos'].values.tolist(), f'Ref start for senterica is not equal'
>>> assert seq_regions_table_true.loc['ref_end_pos'].values.tolist() == seq_regions_table_new.loc['ref_end_pos'].values.tolist(), f'Ref end for senterica is not equal'

```

# Test mutations for senterica and mut_type: chromosomal

```
>>> seq_variations_true = json_data_true['seq_variations']
>>> seq_variations_table_true = pd.DataFrame(seq_variations_true)
>>> seq_variations_new = json_data_new['seq_variations']
>>> seq_variations_table_new = pd.DataFrame(seq_variations_new)
>>> assert seq_variations_table_true.shape == seq_variations_table_new.shape, 'Shape of seq_variations is not equal'
>>> if len(seq_variations_table_true) > 0:
...		assert seq_variations_table_new.loc['key'].values.tolist()[0].split(';')[0] == seq_variations_table_true.loc['key'].values.tolist()[0].split(';')[0], 'First element of key is not equal'
...		assert seq_variations_table_new.loc['ref_start_pos'].values.tolist() == seq_variations_table_true.loc['ref_start_pos'].values.tolist(), 'Mutation start position for chromosomal and senterica is not equal'
...		assert seq_variations_table_new.loc['ref_end_pos'].values.tolist() == seq_variations_table_true.loc['ref_end_pos'].values.tolist(), 'Mutation end position for chromosomal and senterica is not equal'
...		assert seq_variations_table_new.loc['var_codon'].values.tolist() == seq_variations_table_true.loc['var_codon'].values.tolist(), 'Mutation codon for chromosomal and senterica is not equal'
...		assert seq_variations_table_new.loc['ref_codon'].values.tolist() == seq_variations_table_true.loc['ref_codon'].values.tolist(), 'Reference codon for chromosomal and senterica is not equal'

```

### Test if resistant gene is same after splitting in two contigs.

```
>>> with open('/home/people/s220672/resfinder/tests/data/test_isolate_10.fa', 'r') as f:
...		fasta = f.read()
>>> ModFasta = FastaModifier()
>>> both_seqs, both_headers = ModFasta.modify_seqs(json_data_true, '/home/people/s220672/resfinder/tests/data/test_isolate_10.fa')
>>> ModFasta.create_fasta('test_temp/new_fasta', both_seqs, both_headers)
>>> subprocess.run(['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'senterica', '-c', '-ifa', 'test_temp/new_fasta.fa', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_10.json'])
CompletedProcess(args=['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'senterica', '-c', '-ifa', 'test_temp/new_fasta.fa', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_10.json'], returncode=0)
>>> json_file_new = '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_10.json'
>>> with open(json_file_new, 'r') as f:
...		json_data_new = json.load(f)
>>> seq_variations = json_data_true['seq_variations']
>>> seq_variations_table_true = pd.DataFrame(seq_variations)
>>> if len(seq_variations_table_true) > 0:
...		seq_variations_new = json_data_new['seq_variations']
...		seq_variations_table_new = pd.DataFrame(seq_variations_new)
...		assert seq_variations_table_true.loc['seq_var'].values.tolist() == seq_variations_table_new.loc['seq_var'].values.tolist(), f'Reference sequence length for senterica is not equal'
...		assert seq_variations_table_true.loc['ref_end_pos'].values.tolist() == seq_variations_table_new.loc['ref_end_pos'].values.tolist(), f'Ref id for senterica is not equal'

```

```
>>> os.remove('/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_10.json')

```

## Software Command for ecoli

```
>>> subprocess.run(['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'ecoli', '-c', '-ifa', '/home/people/s220672/resfinder/tests/data/test_isolate_11.fa', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_11.json'])
CompletedProcess(args=['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'ecoli', '-c', '-ifa', '/home/people/s220672/resfinder/tests/data/test_isolate_11.fa', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_11.json'], returncode=0)
>>> json_file = '/home/people/s220672/functional_testing/functional_testing_data/ecoli/test_isolate_11/PEOug/test_isolate_11_-c.json'
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

### Test seq_regions for ecoli

```
>>> seq_regions = json_data_true['seq_regions']
>>> seq_regions_table_true = pd.DataFrame(seq_regions)
>>> seq_region_new = json_data_new['seq_regions']
>>> seq_regions_table_new = pd.DataFrame(seq_region_new)
>>> assert seq_regions_table_true.loc['phenotypes'].values.tolist() == seq_regions_table_new.loc['phenotypes'].values.tolist(), 'Phenotypes are not equal'
>>> assert seq_regions_table_true.loc['alignment_length'].values.tolist() == seq_regions_table_new.loc['alignment_length'].values.tolist(), f'Alignment length for ecoli is not equal'
>>> assert seq_regions_table_true.loc['ref_start_pos'].values.tolist() == seq_regions_table_new.loc['ref_start_pos'].values.tolist(), f'Ref start for ecoli is not equal'
>>> assert seq_regions_table_true.loc['ref_end_pos'].values.tolist() == seq_regions_table_new.loc['ref_end_pos'].values.tolist(), f'Ref end for ecoli is not equal'

```

# Test mutations for ecoli and mut_type: chromosomal

```
>>> seq_variations_true = json_data_true['seq_variations']
>>> seq_variations_table_true = pd.DataFrame(seq_variations_true)
>>> seq_variations_new = json_data_new['seq_variations']
>>> seq_variations_table_new = pd.DataFrame(seq_variations_new)
>>> assert seq_variations_table_true.shape == seq_variations_table_new.shape, 'Shape of seq_variations is not equal'
>>> if len(seq_variations_table_true) > 0:
...		assert seq_variations_table_new.loc['key'].values.tolist()[0].split(';')[0] == seq_variations_table_true.loc['key'].values.tolist()[0].split(';')[0], 'First element of key is not equal'
...		assert seq_variations_table_new.loc['ref_start_pos'].values.tolist() == seq_variations_table_true.loc['ref_start_pos'].values.tolist(), 'Mutation start position for chromosomal and ecoli is not equal'
...		assert seq_variations_table_new.loc['ref_end_pos'].values.tolist() == seq_variations_table_true.loc['ref_end_pos'].values.tolist(), 'Mutation end position for chromosomal and ecoli is not equal'
...		assert seq_variations_table_new.loc['var_codon'].values.tolist() == seq_variations_table_true.loc['var_codon'].values.tolist(), 'Mutation codon for chromosomal and ecoli is not equal'
...		assert seq_variations_table_new.loc['ref_codon'].values.tolist() == seq_variations_table_true.loc['ref_codon'].values.tolist(), 'Reference codon for chromosomal and ecoli is not equal'

```

### Test if resistant gene is same after splitting in two contigs.

```
>>> with open('/home/people/s220672/resfinder/tests/data/test_isolate_11.fa', 'r') as f:
...		fasta = f.read()
>>> ModFasta = FastaModifier()
>>> both_seqs, both_headers = ModFasta.modify_seqs(json_data_true, '/home/people/s220672/resfinder/tests/data/test_isolate_11.fa')
>>> ModFasta.create_fasta('test_temp/new_fasta', both_seqs, both_headers)
>>> subprocess.run(['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'ecoli', '-c', '-ifa', 'test_temp/new_fasta.fa', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_11.json'])
CompletedProcess(args=['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'ecoli', '-c', '-ifa', 'test_temp/new_fasta.fa', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_11.json'], returncode=0)
>>> json_file_new = '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_11.json'
>>> with open(json_file_new, 'r') as f:
...		json_data_new = json.load(f)
>>> seq_variations = json_data_true['seq_variations']
>>> seq_variations_table_true = pd.DataFrame(seq_variations)
>>> if len(seq_variations_table_true) > 0:
...		seq_variations_new = json_data_new['seq_variations']
...		seq_variations_table_new = pd.DataFrame(seq_variations_new)
...		assert seq_variations_table_true.loc['seq_var'].values.tolist() == seq_variations_table_new.loc['seq_var'].values.tolist(), f'Reference sequence length for ecoli is not equal'
...		assert seq_variations_table_true.loc['ref_end_pos'].values.tolist() == seq_variations_table_new.loc['ref_end_pos'].values.tolist(), f'Ref id for ecoli is not equal'

```

```
>>> os.remove('/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_11.json')

```

## Software Command for senterica

```
>>> subprocess.run(['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'senterica', '-c', '-ifa', '/home/people/s220672/resfinder/tests/data/test_isolate_11.fa', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_11.json'])
CompletedProcess(args=['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'senterica', '-c', '-ifa', '/home/people/s220672/resfinder/tests/data/test_isolate_11.fa', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_11.json'], returncode=0)
>>> json_file = '/home/people/s220672/functional_testing/functional_testing_data/senterica/test_isolate_11/iGMRE/test_isolate_11_-c.json'
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

### Test seq_regions for senterica

```
>>> seq_regions = json_data_true['seq_regions']
>>> seq_regions_table_true = pd.DataFrame(seq_regions)
>>> seq_region_new = json_data_new['seq_regions']
>>> seq_regions_table_new = pd.DataFrame(seq_region_new)
>>> assert seq_regions_table_true.loc['phenotypes'].values.tolist() == seq_regions_table_new.loc['phenotypes'].values.tolist(), 'Phenotypes are not equal'
>>> assert seq_regions_table_true.loc['alignment_length'].values.tolist() == seq_regions_table_new.loc['alignment_length'].values.tolist(), f'Alignment length for senterica is not equal'
>>> assert seq_regions_table_true.loc['ref_start_pos'].values.tolist() == seq_regions_table_new.loc['ref_start_pos'].values.tolist(), f'Ref start for senterica is not equal'
>>> assert seq_regions_table_true.loc['ref_end_pos'].values.tolist() == seq_regions_table_new.loc['ref_end_pos'].values.tolist(), f'Ref end for senterica is not equal'

```

# Test mutations for senterica and mut_type: chromosomal

```
>>> seq_variations_true = json_data_true['seq_variations']
>>> seq_variations_table_true = pd.DataFrame(seq_variations_true)
>>> seq_variations_new = json_data_new['seq_variations']
>>> seq_variations_table_new = pd.DataFrame(seq_variations_new)
>>> assert seq_variations_table_true.shape == seq_variations_table_new.shape, 'Shape of seq_variations is not equal'
>>> if len(seq_variations_table_true) > 0:
...		assert seq_variations_table_new.loc['key'].values.tolist()[0].split(';')[0] == seq_variations_table_true.loc['key'].values.tolist()[0].split(';')[0], 'First element of key is not equal'
...		assert seq_variations_table_new.loc['ref_start_pos'].values.tolist() == seq_variations_table_true.loc['ref_start_pos'].values.tolist(), 'Mutation start position for chromosomal and senterica is not equal'
...		assert seq_variations_table_new.loc['ref_end_pos'].values.tolist() == seq_variations_table_true.loc['ref_end_pos'].values.tolist(), 'Mutation end position for chromosomal and senterica is not equal'
...		assert seq_variations_table_new.loc['var_codon'].values.tolist() == seq_variations_table_true.loc['var_codon'].values.tolist(), 'Mutation codon for chromosomal and senterica is not equal'
...		assert seq_variations_table_new.loc['ref_codon'].values.tolist() == seq_variations_table_true.loc['ref_codon'].values.tolist(), 'Reference codon for chromosomal and senterica is not equal'

```

```
>>> os.remove('/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_11.json')

```

## Software Command for ecoli

```
>>> subprocess.run(['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'ecoli', '-c', '-ifa', '/home/people/s220672/resfinder/tests/data/test_isolate_03.fa', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_03.json'])
CompletedProcess(args=['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'ecoli', '-c', '-ifa', '/home/people/s220672/resfinder/tests/data/test_isolate_03.fa', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_03.json'], returncode=0)
>>> json_file = '/home/people/s220672/functional_testing/functional_testing_data/ecoli/test_isolate_03/ZvxHw/test_isolate_03_-c.json'
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

### Test seq_regions for ecoli

```
>>> seq_regions = json_data_true['seq_regions']
>>> seq_regions_table_true = pd.DataFrame(seq_regions)
>>> seq_region_new = json_data_new['seq_regions']
>>> seq_regions_table_new = pd.DataFrame(seq_region_new)
>>> assert seq_regions_table_true.loc['phenotypes'].values.tolist() == seq_regions_table_new.loc['phenotypes'].values.tolist(), 'Phenotypes are not equal'
>>> assert seq_regions_table_true.loc['alignment_length'].values.tolist() == seq_regions_table_new.loc['alignment_length'].values.tolist(), f'Alignment length for ecoli is not equal'
>>> assert seq_regions_table_true.loc['ref_start_pos'].values.tolist() == seq_regions_table_new.loc['ref_start_pos'].values.tolist(), f'Ref start for ecoli is not equal'
>>> assert seq_regions_table_true.loc['ref_end_pos'].values.tolist() == seq_regions_table_new.loc['ref_end_pos'].values.tolist(), f'Ref end for ecoli is not equal'

```

# Test mutations for ecoli and mut_type: chromosomal

```
>>> seq_variations_true = json_data_true['seq_variations']
>>> seq_variations_table_true = pd.DataFrame(seq_variations_true)
>>> seq_variations_new = json_data_new['seq_variations']
>>> seq_variations_table_new = pd.DataFrame(seq_variations_new)
>>> assert seq_variations_table_true.shape == seq_variations_table_new.shape, 'Shape of seq_variations is not equal'
>>> if len(seq_variations_table_true) > 0:
...		assert seq_variations_table_new.loc['key'].values.tolist()[0].split(';')[0] == seq_variations_table_true.loc['key'].values.tolist()[0].split(';')[0], 'First element of key is not equal'
...		assert seq_variations_table_new.loc['ref_start_pos'].values.tolist() == seq_variations_table_true.loc['ref_start_pos'].values.tolist(), 'Mutation start position for chromosomal and ecoli is not equal'
...		assert seq_variations_table_new.loc['ref_end_pos'].values.tolist() == seq_variations_table_true.loc['ref_end_pos'].values.tolist(), 'Mutation end position for chromosomal and ecoli is not equal'
...		assert seq_variations_table_new.loc['var_codon'].values.tolist() == seq_variations_table_true.loc['var_codon'].values.tolist(), 'Mutation codon for chromosomal and ecoli is not equal'
...		assert seq_variations_table_new.loc['ref_codon'].values.tolist() == seq_variations_table_true.loc['ref_codon'].values.tolist(), 'Reference codon for chromosomal and ecoli is not equal'

```

```
>>> os.remove('/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_03.json')

```

## Software Command for senterica

```
>>> subprocess.run(['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'senterica', '-c', '-ifa', '/home/people/s220672/resfinder/tests/data/test_isolate_03.fa', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_03.json'])
CompletedProcess(args=['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'senterica', '-c', '-ifa', '/home/people/s220672/resfinder/tests/data/test_isolate_03.fa', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_03.json'], returncode=0)
>>> json_file = '/home/people/s220672/functional_testing/functional_testing_data/senterica/test_isolate_03/VPOub/test_isolate_03_-c.json'
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

### Test seq_regions for senterica

```
>>> seq_regions = json_data_true['seq_regions']
>>> seq_regions_table_true = pd.DataFrame(seq_regions)
>>> seq_region_new = json_data_new['seq_regions']
>>> seq_regions_table_new = pd.DataFrame(seq_region_new)
>>> assert seq_regions_table_true.loc['phenotypes'].values.tolist() == seq_regions_table_new.loc['phenotypes'].values.tolist(), 'Phenotypes are not equal'
>>> assert seq_regions_table_true.loc['alignment_length'].values.tolist() == seq_regions_table_new.loc['alignment_length'].values.tolist(), f'Alignment length for senterica is not equal'
>>> assert seq_regions_table_true.loc['ref_start_pos'].values.tolist() == seq_regions_table_new.loc['ref_start_pos'].values.tolist(), f'Ref start for senterica is not equal'
>>> assert seq_regions_table_true.loc['ref_end_pos'].values.tolist() == seq_regions_table_new.loc['ref_end_pos'].values.tolist(), f'Ref end for senterica is not equal'

```

# Test mutations for senterica and mut_type: chromosomal

```
>>> seq_variations_true = json_data_true['seq_variations']
>>> seq_variations_table_true = pd.DataFrame(seq_variations_true)
>>> seq_variations_new = json_data_new['seq_variations']
>>> seq_variations_table_new = pd.DataFrame(seq_variations_new)
>>> assert seq_variations_table_true.shape == seq_variations_table_new.shape, 'Shape of seq_variations is not equal'
>>> if len(seq_variations_table_true) > 0:
...		assert seq_variations_table_new.loc['key'].values.tolist()[0].split(';')[0] == seq_variations_table_true.loc['key'].values.tolist()[0].split(';')[0], 'First element of key is not equal'
...		assert seq_variations_table_new.loc['ref_start_pos'].values.tolist() == seq_variations_table_true.loc['ref_start_pos'].values.tolist(), 'Mutation start position for chromosomal and senterica is not equal'
...		assert seq_variations_table_new.loc['ref_end_pos'].values.tolist() == seq_variations_table_true.loc['ref_end_pos'].values.tolist(), 'Mutation end position for chromosomal and senterica is not equal'
...		assert seq_variations_table_new.loc['var_codon'].values.tolist() == seq_variations_table_true.loc['var_codon'].values.tolist(), 'Mutation codon for chromosomal and senterica is not equal'
...		assert seq_variations_table_new.loc['ref_codon'].values.tolist() == seq_variations_table_true.loc['ref_codon'].values.tolist(), 'Reference codon for chromosomal and senterica is not equal'

```

```
>>> os.remove('/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_03.json')

```

## Software Command for ecoli

```
>>> subprocess.run(['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'ecoli', '-c', '-ifq', '/home/people/s220672/resfinder/tests/data/test_isolate_05_2.fq', '/home/people/s220672/resfinder/tests/data/test_isolate_05_1.fq', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_05_2.json'])
CompletedProcess(args=['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'ecoli', '-c', '-ifq', '/home/people/s220672/resfinder/tests/data/test_isolate_05_2.fq', '/home/people/s220672/resfinder/tests/data/test_isolate_05_1.fq', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_05_2.json'], returncode=0)
>>> json_file = '/home/people/s220672/functional_testing/functional_testing_data/ecoli/test_isolate_05_2/eMtPA/test_isolate_05_2_-c.json'
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

### Test seq_regions for ecoli

```
>>> seq_regions = json_data_true['seq_regions']
>>> seq_regions_table_true = pd.DataFrame(seq_regions)
>>> seq_region_new = json_data_new['seq_regions']
>>> seq_regions_table_new = pd.DataFrame(seq_region_new)
>>> assert seq_regions_table_true.loc['phenotypes'].values.tolist() == seq_regions_table_new.loc['phenotypes'].values.tolist(), 'Phenotypes are not equal'
>>> assert seq_regions_table_true.loc['alignment_length'].values.tolist() == seq_regions_table_new.loc['alignment_length'].values.tolist(), f'Alignment length for ecoli is not equal'
>>> assert seq_regions_table_true.loc['ref_start_pos'].values.tolist() == seq_regions_table_new.loc['ref_start_pos'].values.tolist(), f'Ref start for ecoli is not equal'
>>> assert seq_regions_table_true.loc['ref_end_pos'].values.tolist() == seq_regions_table_new.loc['ref_end_pos'].values.tolist(), f'Ref end for ecoli is not equal'

```

# Test mutations for ecoli and mut_type: chromosomal

```
>>> seq_variations_true = json_data_true['seq_variations']
>>> seq_variations_table_true = pd.DataFrame(seq_variations_true)
>>> seq_variations_new = json_data_new['seq_variations']
>>> seq_variations_table_new = pd.DataFrame(seq_variations_new)
>>> assert seq_variations_table_true.shape == seq_variations_table_new.shape, 'Shape of seq_variations is not equal'
>>> if len(seq_variations_table_true) > 0:
...		assert seq_variations_table_new.loc['key'].values.tolist()[0].split(';')[0] == seq_variations_table_true.loc['key'].values.tolist()[0].split(';')[0], 'First element of key is not equal'
...		assert seq_variations_table_new.loc['ref_start_pos'].values.tolist() == seq_variations_table_true.loc['ref_start_pos'].values.tolist(), 'Mutation start position for chromosomal and ecoli is not equal'
...		assert seq_variations_table_new.loc['ref_end_pos'].values.tolist() == seq_variations_table_true.loc['ref_end_pos'].values.tolist(), 'Mutation end position for chromosomal and ecoli is not equal'
...		assert seq_variations_table_new.loc['var_codon'].values.tolist() == seq_variations_table_true.loc['var_codon'].values.tolist(), 'Mutation codon for chromosomal and ecoli is not equal'
...		assert seq_variations_table_new.loc['ref_codon'].values.tolist() == seq_variations_table_true.loc['ref_codon'].values.tolist(), 'Reference codon for chromosomal and ecoli is not equal'

```

```
>>> os.remove('/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_05_2.json')

```

## Software Command for senterica

```
>>> subprocess.run(['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'senterica', '-c', '-ifq', '/home/people/s220672/resfinder/tests/data/test_isolate_05_2.fq', '/home/people/s220672/resfinder/tests/data/test_isolate_05_1.fq', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_05_2.json'])
CompletedProcess(args=['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'senterica', '-c', '-ifq', '/home/people/s220672/resfinder/tests/data/test_isolate_05_2.fq', '/home/people/s220672/resfinder/tests/data/test_isolate_05_1.fq', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_05_2.json'], returncode=0)
>>> json_file = '/home/people/s220672/functional_testing/functional_testing_data/senterica/test_isolate_05_2/gAnSx/test_isolate_05_2_-c.json'
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

### Test seq_regions for senterica

```
>>> seq_regions = json_data_true['seq_regions']
>>> seq_regions_table_true = pd.DataFrame(seq_regions)
>>> seq_region_new = json_data_new['seq_regions']
>>> seq_regions_table_new = pd.DataFrame(seq_region_new)
>>> assert seq_regions_table_true.loc['phenotypes'].values.tolist() == seq_regions_table_new.loc['phenotypes'].values.tolist(), 'Phenotypes are not equal'
>>> assert seq_regions_table_true.loc['alignment_length'].values.tolist() == seq_regions_table_new.loc['alignment_length'].values.tolist(), f'Alignment length for senterica is not equal'
>>> assert seq_regions_table_true.loc['ref_start_pos'].values.tolist() == seq_regions_table_new.loc['ref_start_pos'].values.tolist(), f'Ref start for senterica is not equal'
>>> assert seq_regions_table_true.loc['ref_end_pos'].values.tolist() == seq_regions_table_new.loc['ref_end_pos'].values.tolist(), f'Ref end for senterica is not equal'

```

# Test mutations for senterica and mut_type: chromosomal

```
>>> seq_variations_true = json_data_true['seq_variations']
>>> seq_variations_table_true = pd.DataFrame(seq_variations_true)
>>> seq_variations_new = json_data_new['seq_variations']
>>> seq_variations_table_new = pd.DataFrame(seq_variations_new)
>>> assert seq_variations_table_true.shape == seq_variations_table_new.shape, 'Shape of seq_variations is not equal'
>>> if len(seq_variations_table_true) > 0:
...		assert seq_variations_table_new.loc['key'].values.tolist()[0].split(';')[0] == seq_variations_table_true.loc['key'].values.tolist()[0].split(';')[0], 'First element of key is not equal'
...		assert seq_variations_table_new.loc['ref_start_pos'].values.tolist() == seq_variations_table_true.loc['ref_start_pos'].values.tolist(), 'Mutation start position for chromosomal and senterica is not equal'
...		assert seq_variations_table_new.loc['ref_end_pos'].values.tolist() == seq_variations_table_true.loc['ref_end_pos'].values.tolist(), 'Mutation end position for chromosomal and senterica is not equal'
...		assert seq_variations_table_new.loc['var_codon'].values.tolist() == seq_variations_table_true.loc['var_codon'].values.tolist(), 'Mutation codon for chromosomal and senterica is not equal'
...		assert seq_variations_table_new.loc['ref_codon'].values.tolist() == seq_variations_table_true.loc['ref_codon'].values.tolist(), 'Reference codon for chromosomal and senterica is not equal'

```

```
>>> os.remove('/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_05_2.json')

```

## Software Command for ecoli

```
>>> subprocess.run(['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'ecoli', '-c', '-ifa', '/home/people/s220672/resfinder/tests/data/test_isolate_01.fa', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_01.json'])
CompletedProcess(args=['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'ecoli', '-c', '-ifa', '/home/people/s220672/resfinder/tests/data/test_isolate_01.fa', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_01.json'], returncode=0)
>>> json_file = '/home/people/s220672/functional_testing/functional_testing_data/ecoli/test_isolate_01/BUiYz/test_isolate_01_-c.json'
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

### Test seq_regions for ecoli

```
>>> seq_regions = json_data_true['seq_regions']
>>> seq_regions_table_true = pd.DataFrame(seq_regions)
>>> seq_region_new = json_data_new['seq_regions']
>>> seq_regions_table_new = pd.DataFrame(seq_region_new)
>>> assert seq_regions_table_true.loc['phenotypes'].values.tolist() == seq_regions_table_new.loc['phenotypes'].values.tolist(), 'Phenotypes are not equal'
>>> assert seq_regions_table_true.loc['alignment_length'].values.tolist() == seq_regions_table_new.loc['alignment_length'].values.tolist(), f'Alignment length for ecoli is not equal'
>>> assert seq_regions_table_true.loc['ref_start_pos'].values.tolist() == seq_regions_table_new.loc['ref_start_pos'].values.tolist(), f'Ref start for ecoli is not equal'
>>> assert seq_regions_table_true.loc['ref_end_pos'].values.tolist() == seq_regions_table_new.loc['ref_end_pos'].values.tolist(), f'Ref end for ecoli is not equal'

```

# Test mutations for ecoli and mut_type: chromosomal

```
>>> seq_variations_true = json_data_true['seq_variations']
>>> seq_variations_table_true = pd.DataFrame(seq_variations_true)
>>> seq_variations_new = json_data_new['seq_variations']
>>> seq_variations_table_new = pd.DataFrame(seq_variations_new)
>>> assert seq_variations_table_true.shape == seq_variations_table_new.shape, 'Shape of seq_variations is not equal'
>>> if len(seq_variations_table_true) > 0:
...		assert seq_variations_table_new.loc['key'].values.tolist()[0].split(';')[0] == seq_variations_table_true.loc['key'].values.tolist()[0].split(';')[0], 'First element of key is not equal'
...		assert seq_variations_table_new.loc['ref_start_pos'].values.tolist() == seq_variations_table_true.loc['ref_start_pos'].values.tolist(), 'Mutation start position for chromosomal and ecoli is not equal'
...		assert seq_variations_table_new.loc['ref_end_pos'].values.tolist() == seq_variations_table_true.loc['ref_end_pos'].values.tolist(), 'Mutation end position for chromosomal and ecoli is not equal'
...		assert seq_variations_table_new.loc['var_codon'].values.tolist() == seq_variations_table_true.loc['var_codon'].values.tolist(), 'Mutation codon for chromosomal and ecoli is not equal'
...		assert seq_variations_table_new.loc['ref_codon'].values.tolist() == seq_variations_table_true.loc['ref_codon'].values.tolist(), 'Reference codon for chromosomal and ecoli is not equal'

```

```
>>> os.remove('/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_01.json')

```

## Software Command for senterica

```
>>> subprocess.run(['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'senterica', '-c', '-ifa', '/home/people/s220672/resfinder/tests/data/test_isolate_01.fa', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_01.json'])
CompletedProcess(args=['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'senterica', '-c', '-ifa', '/home/people/s220672/resfinder/tests/data/test_isolate_01.fa', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_01.json'], returncode=0)
>>> json_file = '/home/people/s220672/functional_testing/functional_testing_data/senterica/test_isolate_01/LWNxW/test_isolate_01_-c.json'
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

### Test seq_regions for senterica

```
>>> seq_regions = json_data_true['seq_regions']
>>> seq_regions_table_true = pd.DataFrame(seq_regions)
>>> seq_region_new = json_data_new['seq_regions']
>>> seq_regions_table_new = pd.DataFrame(seq_region_new)
>>> assert seq_regions_table_true.loc['phenotypes'].values.tolist() == seq_regions_table_new.loc['phenotypes'].values.tolist(), 'Phenotypes are not equal'
>>> assert seq_regions_table_true.loc['alignment_length'].values.tolist() == seq_regions_table_new.loc['alignment_length'].values.tolist(), f'Alignment length for senterica is not equal'
>>> assert seq_regions_table_true.loc['ref_start_pos'].values.tolist() == seq_regions_table_new.loc['ref_start_pos'].values.tolist(), f'Ref start for senterica is not equal'
>>> assert seq_regions_table_true.loc['ref_end_pos'].values.tolist() == seq_regions_table_new.loc['ref_end_pos'].values.tolist(), f'Ref end for senterica is not equal'

```

# Test mutations for senterica and mut_type: chromosomal

```
>>> seq_variations_true = json_data_true['seq_variations']
>>> seq_variations_table_true = pd.DataFrame(seq_variations_true)
>>> seq_variations_new = json_data_new['seq_variations']
>>> seq_variations_table_new = pd.DataFrame(seq_variations_new)
>>> assert seq_variations_table_true.shape == seq_variations_table_new.shape, 'Shape of seq_variations is not equal'
>>> if len(seq_variations_table_true) > 0:
...		assert seq_variations_table_new.loc['key'].values.tolist()[0].split(';')[0] == seq_variations_table_true.loc['key'].values.tolist()[0].split(';')[0], 'First element of key is not equal'
...		assert seq_variations_table_new.loc['ref_start_pos'].values.tolist() == seq_variations_table_true.loc['ref_start_pos'].values.tolist(), 'Mutation start position for chromosomal and senterica is not equal'
...		assert seq_variations_table_new.loc['ref_end_pos'].values.tolist() == seq_variations_table_true.loc['ref_end_pos'].values.tolist(), 'Mutation end position for chromosomal and senterica is not equal'
...		assert seq_variations_table_new.loc['var_codon'].values.tolist() == seq_variations_table_true.loc['var_codon'].values.tolist(), 'Mutation codon for chromosomal and senterica is not equal'
...		assert seq_variations_table_new.loc['ref_codon'].values.tolist() == seq_variations_table_true.loc['ref_codon'].values.tolist(), 'Reference codon for chromosomal and senterica is not equal'

```

```
>>> os.remove('/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_01.json')

```

## Software Command for ecoli

```
>>> subprocess.run(['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'ecoli', '-acq', '-ifa', '/home/people/s220672/resfinder/tests/data/test_isolate_01.fa', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_01.json'])
CompletedProcess(args=['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'ecoli', '-acq', '-ifa', '/home/people/s220672/resfinder/tests/data/test_isolate_01.fa', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_01.json'], returncode=0)
>>> json_file = '/home/people/s220672/functional_testing/functional_testing_data/ecoli/test_isolate_01/sZxzt/test_isolate_01_-acq.json'
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

### Test seq_regions for ecoli

```
>>> seq_regions = json_data_true['seq_regions']
>>> seq_regions_table_true = pd.DataFrame(seq_regions)
>>> seq_region_new = json_data_new['seq_regions']
>>> seq_regions_table_new = pd.DataFrame(seq_region_new)
>>> assert seq_regions_table_true.loc['phenotypes'].values.tolist() == seq_regions_table_new.loc['phenotypes'].values.tolist(), 'Phenotypes are not equal'
>>> assert seq_regions_table_true.loc['alignment_length'].values.tolist() == seq_regions_table_new.loc['alignment_length'].values.tolist(), f'Alignment length for ecoli is not equal'
>>> assert seq_regions_table_true.loc['ref_start_pos'].values.tolist() == seq_regions_table_new.loc['ref_start_pos'].values.tolist(), f'Ref start for ecoli is not equal'
>>> assert seq_regions_table_true.loc['ref_end_pos'].values.tolist() == seq_regions_table_new.loc['ref_end_pos'].values.tolist(), f'Ref end for ecoli is not equal'

```

# Test mutations for ecoli and mut_type: acquired

```
>>> seq_variations_true = json_data_true['seq_variations']
>>> seq_variations_table_true = pd.DataFrame(seq_variations_true)
>>> seq_variations_new = json_data_new['seq_variations']
>>> seq_variations_table_new = pd.DataFrame(seq_variations_new)
>>> assert seq_variations_table_true.shape == seq_variations_table_new.shape, 'Shape of seq_variations is not equal'
>>> if len(seq_variations_table_true) > 0:
...		assert seq_variations_table_new.loc['key'].values.tolist()[0].split(';')[0] == seq_variations_table_true.loc['key'].values.tolist()[0].split(';')[0], 'First element of key is not equal'
...		assert seq_variations_table_new.loc['ref_start_pos'].values.tolist() == seq_variations_table_true.loc['ref_start_pos'].values.tolist(), 'Mutation start position for acquired and ecoli is not equal'
...		assert seq_variations_table_new.loc['ref_end_pos'].values.tolist() == seq_variations_table_true.loc['ref_end_pos'].values.tolist(), 'Mutation end position for acquired and ecoli is not equal'
...		assert seq_variations_table_new.loc['var_codon'].values.tolist() == seq_variations_table_true.loc['var_codon'].values.tolist(), 'Mutation codon for acquired and ecoli is not equal'
...		assert seq_variations_table_new.loc['ref_codon'].values.tolist() == seq_variations_table_true.loc['ref_codon'].values.tolist(), 'Reference codon for acquired and ecoli is not equal'

```

```
>>> os.remove('/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_01.json')

```

## Software Command for cjejuni

```
>>> subprocess.run(['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'cjejuni', '-acq', '-ifa', '/home/people/s220672/resfinder/tests/data/test_isolate_01.fa', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_01.json'])
CompletedProcess(args=['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'cjejuni', '-acq', '-ifa', '/home/people/s220672/resfinder/tests/data/test_isolate_01.fa', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_01.json'], returncode=0)
>>> json_file = '/home/people/s220672/functional_testing/functional_testing_data/cjejuni/test_isolate_01/GaEjd/test_isolate_01_-acq.json'
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

### Test seq_regions for cjejuni

```
>>> seq_regions = json_data_true['seq_regions']
>>> seq_regions_table_true = pd.DataFrame(seq_regions)
>>> seq_region_new = json_data_new['seq_regions']
>>> seq_regions_table_new = pd.DataFrame(seq_region_new)
>>> assert seq_regions_table_true.loc['phenotypes'].values.tolist() == seq_regions_table_new.loc['phenotypes'].values.tolist(), 'Phenotypes are not equal'
>>> assert seq_regions_table_true.loc['alignment_length'].values.tolist() == seq_regions_table_new.loc['alignment_length'].values.tolist(), f'Alignment length for cjejuni is not equal'
>>> assert seq_regions_table_true.loc['ref_start_pos'].values.tolist() == seq_regions_table_new.loc['ref_start_pos'].values.tolist(), f'Ref start for cjejuni is not equal'
>>> assert seq_regions_table_true.loc['ref_end_pos'].values.tolist() == seq_regions_table_new.loc['ref_end_pos'].values.tolist(), f'Ref end for cjejuni is not equal'

```

# Test mutations for cjejuni and mut_type: acquired

```
>>> seq_variations_true = json_data_true['seq_variations']
>>> seq_variations_table_true = pd.DataFrame(seq_variations_true)
>>> seq_variations_new = json_data_new['seq_variations']
>>> seq_variations_table_new = pd.DataFrame(seq_variations_new)
>>> assert seq_variations_table_true.shape == seq_variations_table_new.shape, 'Shape of seq_variations is not equal'
>>> if len(seq_variations_table_true) > 0:
...		assert seq_variations_table_new.loc['key'].values.tolist()[0].split(';')[0] == seq_variations_table_true.loc['key'].values.tolist()[0].split(';')[0], 'First element of key is not equal'
...		assert seq_variations_table_new.loc['ref_start_pos'].values.tolist() == seq_variations_table_true.loc['ref_start_pos'].values.tolist(), 'Mutation start position for acquired and cjejuni is not equal'
...		assert seq_variations_table_new.loc['ref_end_pos'].values.tolist() == seq_variations_table_true.loc['ref_end_pos'].values.tolist(), 'Mutation end position for acquired and cjejuni is not equal'
...		assert seq_variations_table_new.loc['var_codon'].values.tolist() == seq_variations_table_true.loc['var_codon'].values.tolist(), 'Mutation codon for acquired and cjejuni is not equal'
...		assert seq_variations_table_new.loc['ref_codon'].values.tolist() == seq_variations_table_true.loc['ref_codon'].values.tolist(), 'Reference codon for acquired and cjejuni is not equal'

```

```
>>> os.remove('/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_01.json')

```

## Software Command for ccoli

```
>>> subprocess.run(['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'ccoli', '-acq', '-ifa', '/home/people/s220672/resfinder/tests/data/test_isolate_01.fa', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_01.json'])
CompletedProcess(args=['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'ccoli', '-acq', '-ifa', '/home/people/s220672/resfinder/tests/data/test_isolate_01.fa', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_01.json'], returncode=0)
>>> json_file = '/home/people/s220672/functional_testing/functional_testing_data/ccoli/test_isolate_01/abSeU/test_isolate_01_-acq.json'
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

### Test seq_regions for ccoli

```
>>> seq_regions = json_data_true['seq_regions']
>>> seq_regions_table_true = pd.DataFrame(seq_regions)
>>> seq_region_new = json_data_new['seq_regions']
>>> seq_regions_table_new = pd.DataFrame(seq_region_new)
>>> assert seq_regions_table_true.loc['phenotypes'].values.tolist() == seq_regions_table_new.loc['phenotypes'].values.tolist(), 'Phenotypes are not equal'
>>> assert seq_regions_table_true.loc['alignment_length'].values.tolist() == seq_regions_table_new.loc['alignment_length'].values.tolist(), f'Alignment length for ccoli is not equal'
>>> assert seq_regions_table_true.loc['ref_start_pos'].values.tolist() == seq_regions_table_new.loc['ref_start_pos'].values.tolist(), f'Ref start for ccoli is not equal'
>>> assert seq_regions_table_true.loc['ref_end_pos'].values.tolist() == seq_regions_table_new.loc['ref_end_pos'].values.tolist(), f'Ref end for ccoli is not equal'

```

# Test mutations for ccoli and mut_type: acquired

```
>>> seq_variations_true = json_data_true['seq_variations']
>>> seq_variations_table_true = pd.DataFrame(seq_variations_true)
>>> seq_variations_new = json_data_new['seq_variations']
>>> seq_variations_table_new = pd.DataFrame(seq_variations_new)
>>> assert seq_variations_table_true.shape == seq_variations_table_new.shape, 'Shape of seq_variations is not equal'
>>> if len(seq_variations_table_true) > 0:
...		assert seq_variations_table_new.loc['key'].values.tolist()[0].split(';')[0] == seq_variations_table_true.loc['key'].values.tolist()[0].split(';')[0], 'First element of key is not equal'
...		assert seq_variations_table_new.loc['ref_start_pos'].values.tolist() == seq_variations_table_true.loc['ref_start_pos'].values.tolist(), 'Mutation start position for acquired and ccoli is not equal'
...		assert seq_variations_table_new.loc['ref_end_pos'].values.tolist() == seq_variations_table_true.loc['ref_end_pos'].values.tolist(), 'Mutation end position for acquired and ccoli is not equal'
...		assert seq_variations_table_new.loc['var_codon'].values.tolist() == seq_variations_table_true.loc['var_codon'].values.tolist(), 'Mutation codon for acquired and ccoli is not equal'
...		assert seq_variations_table_new.loc['ref_codon'].values.tolist() == seq_variations_table_true.loc['ref_codon'].values.tolist(), 'Reference codon for acquired and ccoli is not equal'

```

```
>>> os.remove('/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_01.json')

```

## Software Command for senterica

```
>>> subprocess.run(['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'senterica', '-acq', '-ifa', '/home/people/s220672/resfinder/tests/data/test_isolate_01.fa', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_01.json'])
CompletedProcess(args=['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'senterica', '-acq', '-ifa', '/home/people/s220672/resfinder/tests/data/test_isolate_01.fa', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_01.json'], returncode=0)
>>> json_file = '/home/people/s220672/functional_testing/functional_testing_data/senterica/test_isolate_01/ztFMu/test_isolate_01_-acq.json'
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

### Test seq_regions for senterica

```
>>> seq_regions = json_data_true['seq_regions']
>>> seq_regions_table_true = pd.DataFrame(seq_regions)
>>> seq_region_new = json_data_new['seq_regions']
>>> seq_regions_table_new = pd.DataFrame(seq_region_new)
>>> assert seq_regions_table_true.loc['phenotypes'].values.tolist() == seq_regions_table_new.loc['phenotypes'].values.tolist(), 'Phenotypes are not equal'
>>> assert seq_regions_table_true.loc['alignment_length'].values.tolist() == seq_regions_table_new.loc['alignment_length'].values.tolist(), f'Alignment length for senterica is not equal'
>>> assert seq_regions_table_true.loc['ref_start_pos'].values.tolist() == seq_regions_table_new.loc['ref_start_pos'].values.tolist(), f'Ref start for senterica is not equal'
>>> assert seq_regions_table_true.loc['ref_end_pos'].values.tolist() == seq_regions_table_new.loc['ref_end_pos'].values.tolist(), f'Ref end for senterica is not equal'

```

# Test mutations for senterica and mut_type: acquired

```
>>> seq_variations_true = json_data_true['seq_variations']
>>> seq_variations_table_true = pd.DataFrame(seq_variations_true)
>>> seq_variations_new = json_data_new['seq_variations']
>>> seq_variations_table_new = pd.DataFrame(seq_variations_new)
>>> assert seq_variations_table_true.shape == seq_variations_table_new.shape, 'Shape of seq_variations is not equal'
>>> if len(seq_variations_table_true) > 0:
...		assert seq_variations_table_new.loc['key'].values.tolist()[0].split(';')[0] == seq_variations_table_true.loc['key'].values.tolist()[0].split(';')[0], 'First element of key is not equal'
...		assert seq_variations_table_new.loc['ref_start_pos'].values.tolist() == seq_variations_table_true.loc['ref_start_pos'].values.tolist(), 'Mutation start position for acquired and senterica is not equal'
...		assert seq_variations_table_new.loc['ref_end_pos'].values.tolist() == seq_variations_table_true.loc['ref_end_pos'].values.tolist(), 'Mutation end position for acquired and senterica is not equal'
...		assert seq_variations_table_new.loc['var_codon'].values.tolist() == seq_variations_table_true.loc['var_codon'].values.tolist(), 'Mutation codon for acquired and senterica is not equal'
...		assert seq_variations_table_new.loc['ref_codon'].values.tolist() == seq_variations_table_true.loc['ref_codon'].values.tolist(), 'Reference codon for acquired and senterica is not equal'

```

```
>>> os.remove('/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_01.json')

```

## Software Command for ecoli

```
>>> subprocess.run(['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'ecoli', '-u', '-acq', '-ifa', '/home/people/s220672/resfinder/tests/data/test_isolate_01.fa', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_01.json'])
CompletedProcess(args=['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'ecoli', '-u', '-acq', '-ifa', '/home/people/s220672/resfinder/tests/data/test_isolate_01.fa', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_01.json'], returncode=0)
>>> json_file = '/home/people/s220672/functional_testing/functional_testing_data/ecoli/test_isolate_01/fUTri/test_isolate_01_-acq-u.json'
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

### Test seq_regions for ecoli

```
>>> seq_regions = json_data_true['seq_regions']
>>> seq_regions_table_true = pd.DataFrame(seq_regions)
>>> seq_region_new = json_data_new['seq_regions']
>>> seq_regions_table_new = pd.DataFrame(seq_region_new)
>>> assert seq_regions_table_true.loc['phenotypes'].values.tolist() == seq_regions_table_new.loc['phenotypes'].values.tolist(), 'Phenotypes are not equal'
>>> assert seq_regions_table_true.loc['alignment_length'].values.tolist() == seq_regions_table_new.loc['alignment_length'].values.tolist(), f'Alignment length for ecoli is not equal'
>>> assert seq_regions_table_true.loc['ref_start_pos'].values.tolist() == seq_regions_table_new.loc['ref_start_pos'].values.tolist(), f'Ref start for ecoli is not equal'
>>> assert seq_regions_table_true.loc['ref_end_pos'].values.tolist() == seq_regions_table_new.loc['ref_end_pos'].values.tolist(), f'Ref end for ecoli is not equal'

```

# Test mutations for ecoli and mut_type: unknown

```
>>> seq_variations_true = json_data_true['seq_variations']
>>> seq_variations_table_true = pd.DataFrame(seq_variations_true)
>>> seq_variations_new = json_data_new['seq_variations']
>>> seq_variations_table_new = pd.DataFrame(seq_variations_new)
>>> assert seq_variations_table_true.shape == seq_variations_table_new.shape, 'Shape of seq_variations is not equal'
>>> if len(seq_variations_table_true) > 0:
...		assert seq_variations_table_new.loc['key'].values.tolist()[0].split(';')[0] == seq_variations_table_true.loc['key'].values.tolist()[0].split(';')[0], 'First element of key is not equal'
...		assert seq_variations_table_new.loc['ref_start_pos'].values.tolist() == seq_variations_table_true.loc['ref_start_pos'].values.tolist(), 'Mutation start position for unknown and ecoli is not equal'
...		assert seq_variations_table_new.loc['ref_end_pos'].values.tolist() == seq_variations_table_true.loc['ref_end_pos'].values.tolist(), 'Mutation end position for unknown and ecoli is not equal'
...		assert seq_variations_table_new.loc['var_codon'].values.tolist() == seq_variations_table_true.loc['var_codon'].values.tolist(), 'Mutation codon for unknown and ecoli is not equal'
...		assert seq_variations_table_new.loc['ref_codon'].values.tolist() == seq_variations_table_true.loc['ref_codon'].values.tolist(), 'Reference codon for unknown and ecoli is not equal'

```

```
>>> os.remove('/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_01.json')

```

## Software Command for cjejuni

```
>>> subprocess.run(['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'cjejuni', '-u', '-acq', '-ifa', '/home/people/s220672/resfinder/tests/data/test_isolate_01.fa', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_01.json'])
CompletedProcess(args=['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'cjejuni', '-u', '-acq', '-ifa', '/home/people/s220672/resfinder/tests/data/test_isolate_01.fa', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_01.json'], returncode=0)
>>> json_file = '/home/people/s220672/functional_testing/functional_testing_data/cjejuni/test_isolate_01/knygJ/test_isolate_01_-acq-u.json'
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

### Test seq_regions for cjejuni

```
>>> seq_regions = json_data_true['seq_regions']
>>> seq_regions_table_true = pd.DataFrame(seq_regions)
>>> seq_region_new = json_data_new['seq_regions']
>>> seq_regions_table_new = pd.DataFrame(seq_region_new)
>>> assert seq_regions_table_true.loc['phenotypes'].values.tolist() == seq_regions_table_new.loc['phenotypes'].values.tolist(), 'Phenotypes are not equal'
>>> assert seq_regions_table_true.loc['alignment_length'].values.tolist() == seq_regions_table_new.loc['alignment_length'].values.tolist(), f'Alignment length for cjejuni is not equal'
>>> assert seq_regions_table_true.loc['ref_start_pos'].values.tolist() == seq_regions_table_new.loc['ref_start_pos'].values.tolist(), f'Ref start for cjejuni is not equal'
>>> assert seq_regions_table_true.loc['ref_end_pos'].values.tolist() == seq_regions_table_new.loc['ref_end_pos'].values.tolist(), f'Ref end for cjejuni is not equal'

```

# Test mutations for cjejuni and mut_type: unknown

```
>>> seq_variations_true = json_data_true['seq_variations']
>>> seq_variations_table_true = pd.DataFrame(seq_variations_true)
>>> seq_variations_new = json_data_new['seq_variations']
>>> seq_variations_table_new = pd.DataFrame(seq_variations_new)
>>> assert seq_variations_table_true.shape == seq_variations_table_new.shape, 'Shape of seq_variations is not equal'
>>> if len(seq_variations_table_true) > 0:
...		assert seq_variations_table_new.loc['key'].values.tolist()[0].split(';')[0] == seq_variations_table_true.loc['key'].values.tolist()[0].split(';')[0], 'First element of key is not equal'
...		assert seq_variations_table_new.loc['ref_start_pos'].values.tolist() == seq_variations_table_true.loc['ref_start_pos'].values.tolist(), 'Mutation start position for unknown and cjejuni is not equal'
...		assert seq_variations_table_new.loc['ref_end_pos'].values.tolist() == seq_variations_table_true.loc['ref_end_pos'].values.tolist(), 'Mutation end position for unknown and cjejuni is not equal'
...		assert seq_variations_table_new.loc['var_codon'].values.tolist() == seq_variations_table_true.loc['var_codon'].values.tolist(), 'Mutation codon for unknown and cjejuni is not equal'
...		assert seq_variations_table_new.loc['ref_codon'].values.tolist() == seq_variations_table_true.loc['ref_codon'].values.tolist(), 'Reference codon for unknown and cjejuni is not equal'

```

```
>>> os.remove('/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_01.json')

```

## Software Command for ccoli

```
>>> subprocess.run(['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'ccoli', '-u', '-acq', '-ifa', '/home/people/s220672/resfinder/tests/data/test_isolate_01.fa', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_01.json'])
CompletedProcess(args=['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'ccoli', '-u', '-acq', '-ifa', '/home/people/s220672/resfinder/tests/data/test_isolate_01.fa', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_01.json'], returncode=0)
>>> json_file = '/home/people/s220672/functional_testing/functional_testing_data/ccoli/test_isolate_01/DahyA/test_isolate_01_-acq-u.json'
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

### Test seq_regions for ccoli

```
>>> seq_regions = json_data_true['seq_regions']
>>> seq_regions_table_true = pd.DataFrame(seq_regions)
>>> seq_region_new = json_data_new['seq_regions']
>>> seq_regions_table_new = pd.DataFrame(seq_region_new)
>>> assert seq_regions_table_true.loc['phenotypes'].values.tolist() == seq_regions_table_new.loc['phenotypes'].values.tolist(), 'Phenotypes are not equal'
>>> assert seq_regions_table_true.loc['alignment_length'].values.tolist() == seq_regions_table_new.loc['alignment_length'].values.tolist(), f'Alignment length for ccoli is not equal'
>>> assert seq_regions_table_true.loc['ref_start_pos'].values.tolist() == seq_regions_table_new.loc['ref_start_pos'].values.tolist(), f'Ref start for ccoli is not equal'
>>> assert seq_regions_table_true.loc['ref_end_pos'].values.tolist() == seq_regions_table_new.loc['ref_end_pos'].values.tolist(), f'Ref end for ccoli is not equal'

```

# Test mutations for ccoli and mut_type: unknown

```
>>> seq_variations_true = json_data_true['seq_variations']
>>> seq_variations_table_true = pd.DataFrame(seq_variations_true)
>>> seq_variations_new = json_data_new['seq_variations']
>>> seq_variations_table_new = pd.DataFrame(seq_variations_new)
>>> assert seq_variations_table_true.shape == seq_variations_table_new.shape, 'Shape of seq_variations is not equal'
>>> if len(seq_variations_table_true) > 0:
...		assert seq_variations_table_new.loc['key'].values.tolist()[0].split(';')[0] == seq_variations_table_true.loc['key'].values.tolist()[0].split(';')[0], 'First element of key is not equal'
...		assert seq_variations_table_new.loc['ref_start_pos'].values.tolist() == seq_variations_table_true.loc['ref_start_pos'].values.tolist(), 'Mutation start position for unknown and ccoli is not equal'
...		assert seq_variations_table_new.loc['ref_end_pos'].values.tolist() == seq_variations_table_true.loc['ref_end_pos'].values.tolist(), 'Mutation end position for unknown and ccoli is not equal'
...		assert seq_variations_table_new.loc['var_codon'].values.tolist() == seq_variations_table_true.loc['var_codon'].values.tolist(), 'Mutation codon for unknown and ccoli is not equal'
...		assert seq_variations_table_new.loc['ref_codon'].values.tolist() == seq_variations_table_true.loc['ref_codon'].values.tolist(), 'Reference codon for unknown and ccoli is not equal'

```

```
>>> os.remove('/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_01.json')

```

## Software Command for senterica

```
>>> subprocess.run(['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'senterica', '-u', '-acq', '-ifa', '/home/people/s220672/resfinder/tests/data/test_isolate_01.fa', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_01.json'])
CompletedProcess(args=['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'senterica', '-u', '-acq', '-ifa', '/home/people/s220672/resfinder/tests/data/test_isolate_01.fa', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_01.json'], returncode=0)
>>> json_file = '/home/people/s220672/functional_testing/functional_testing_data/senterica/test_isolate_01/sklSP/test_isolate_01_-acq-u.json'
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

### Test seq_regions for senterica

```
>>> seq_regions = json_data_true['seq_regions']
>>> seq_regions_table_true = pd.DataFrame(seq_regions)
>>> seq_region_new = json_data_new['seq_regions']
>>> seq_regions_table_new = pd.DataFrame(seq_region_new)
>>> assert seq_regions_table_true.loc['phenotypes'].values.tolist() == seq_regions_table_new.loc['phenotypes'].values.tolist(), 'Phenotypes are not equal'
>>> assert seq_regions_table_true.loc['alignment_length'].values.tolist() == seq_regions_table_new.loc['alignment_length'].values.tolist(), f'Alignment length for senterica is not equal'
>>> assert seq_regions_table_true.loc['ref_start_pos'].values.tolist() == seq_regions_table_new.loc['ref_start_pos'].values.tolist(), f'Ref start for senterica is not equal'
>>> assert seq_regions_table_true.loc['ref_end_pos'].values.tolist() == seq_regions_table_new.loc['ref_end_pos'].values.tolist(), f'Ref end for senterica is not equal'

```

# Test mutations for senterica and mut_type: unknown

```
>>> seq_variations_true = json_data_true['seq_variations']
>>> seq_variations_table_true = pd.DataFrame(seq_variations_true)
>>> seq_variations_new = json_data_new['seq_variations']
>>> seq_variations_table_new = pd.DataFrame(seq_variations_new)
>>> assert seq_variations_table_true.shape == seq_variations_table_new.shape, 'Shape of seq_variations is not equal'
>>> if len(seq_variations_table_true) > 0:
...		assert seq_variations_table_new.loc['key'].values.tolist()[0].split(';')[0] == seq_variations_table_true.loc['key'].values.tolist()[0].split(';')[0], 'First element of key is not equal'
...		assert seq_variations_table_new.loc['ref_start_pos'].values.tolist() == seq_variations_table_true.loc['ref_start_pos'].values.tolist(), 'Mutation start position for unknown and senterica is not equal'
...		assert seq_variations_table_new.loc['ref_end_pos'].values.tolist() == seq_variations_table_true.loc['ref_end_pos'].values.tolist(), 'Mutation end position for unknown and senterica is not equal'
...		assert seq_variations_table_new.loc['var_codon'].values.tolist() == seq_variations_table_true.loc['var_codon'].values.tolist(), 'Mutation codon for unknown and senterica is not equal'
...		assert seq_variations_table_new.loc['ref_codon'].values.tolist() == seq_variations_table_true.loc['ref_codon'].values.tolist(), 'Reference codon for unknown and senterica is not equal'

```

```
>>> os.remove('/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_01.json')

```

## Software Command for ecoli

```
>>> subprocess.run(['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'ecoli', '-c', '-ifq', '/home/people/s220672/resfinder/tests/data/test_isolate_01_2.fq', '/home/people/s220672/resfinder/tests/data/test_isolate_01_1.fq', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_01_2.json'])
CompletedProcess(args=['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'ecoli', '-c', '-ifq', '/home/people/s220672/resfinder/tests/data/test_isolate_01_2.fq', '/home/people/s220672/resfinder/tests/data/test_isolate_01_1.fq', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_01_2.json'], returncode=0)
>>> json_file = '/home/people/s220672/functional_testing/functional_testing_data/ecoli/test_isolate_01_2/zmBeE/test_isolate_01_2_-c.json'
>>> with open(json_file, 'r') as f:
...		json_data_true = json.load(f)
>>> json_file_new = '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_01_2.json'
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

### Test seq_regions for ecoli

```
>>> seq_regions = json_data_true['seq_regions']
>>> seq_regions_table_true = pd.DataFrame(seq_regions)
>>> seq_region_new = json_data_new['seq_regions']
>>> seq_regions_table_new = pd.DataFrame(seq_region_new)
>>> assert seq_regions_table_true.loc['phenotypes'].values.tolist() == seq_regions_table_new.loc['phenotypes'].values.tolist(), 'Phenotypes are not equal'
>>> assert seq_regions_table_true.loc['alignment_length'].values.tolist() == seq_regions_table_new.loc['alignment_length'].values.tolist(), f'Alignment length for ecoli is not equal'
>>> assert seq_regions_table_true.loc['ref_start_pos'].values.tolist() == seq_regions_table_new.loc['ref_start_pos'].values.tolist(), f'Ref start for ecoli is not equal'
>>> assert seq_regions_table_true.loc['ref_end_pos'].values.tolist() == seq_regions_table_new.loc['ref_end_pos'].values.tolist(), f'Ref end for ecoli is not equal'

```

# Test mutations for ecoli and mut_type: chromosomal

```
>>> seq_variations_true = json_data_true['seq_variations']
>>> seq_variations_table_true = pd.DataFrame(seq_variations_true)
>>> seq_variations_new = json_data_new['seq_variations']
>>> seq_variations_table_new = pd.DataFrame(seq_variations_new)
>>> assert seq_variations_table_true.shape == seq_variations_table_new.shape, 'Shape of seq_variations is not equal'
>>> if len(seq_variations_table_true) > 0:
...		assert seq_variations_table_new.loc['key'].values.tolist()[0].split(';')[0] == seq_variations_table_true.loc['key'].values.tolist()[0].split(';')[0], 'First element of key is not equal'
...		assert seq_variations_table_new.loc['ref_start_pos'].values.tolist() == seq_variations_table_true.loc['ref_start_pos'].values.tolist(), 'Mutation start position for chromosomal and ecoli is not equal'
...		assert seq_variations_table_new.loc['ref_end_pos'].values.tolist() == seq_variations_table_true.loc['ref_end_pos'].values.tolist(), 'Mutation end position for chromosomal and ecoli is not equal'
...		assert seq_variations_table_new.loc['var_codon'].values.tolist() == seq_variations_table_true.loc['var_codon'].values.tolist(), 'Mutation codon for chromosomal and ecoli is not equal'
...		assert seq_variations_table_new.loc['ref_codon'].values.tolist() == seq_variations_table_true.loc['ref_codon'].values.tolist(), 'Reference codon for chromosomal and ecoli is not equal'

```

```
>>> os.remove('/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_01_2.json')

```

## Software Command for senterica

```
>>> subprocess.run(['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'senterica', '-c', '-ifq', '/home/people/s220672/resfinder/tests/data/test_isolate_01_2.fq', '/home/people/s220672/resfinder/tests/data/test_isolate_01_1.fq', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_01_2.json'])
CompletedProcess(args=['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'senterica', '-c', '-ifq', '/home/people/s220672/resfinder/tests/data/test_isolate_01_2.fq', '/home/people/s220672/resfinder/tests/data/test_isolate_01_1.fq', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_01_2.json'], returncode=0)
>>> json_file = '/home/people/s220672/functional_testing/functional_testing_data/senterica/test_isolate_01_2/Gyozy/test_isolate_01_2_-c.json'
>>> with open(json_file, 'r') as f:
...		json_data_true = json.load(f)
>>> json_file_new = '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_01_2.json'
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

### Test seq_regions for senterica

```
>>> seq_regions = json_data_true['seq_regions']
>>> seq_regions_table_true = pd.DataFrame(seq_regions)
>>> seq_region_new = json_data_new['seq_regions']
>>> seq_regions_table_new = pd.DataFrame(seq_region_new)
>>> assert seq_regions_table_true.loc['phenotypes'].values.tolist() == seq_regions_table_new.loc['phenotypes'].values.tolist(), 'Phenotypes are not equal'
>>> assert seq_regions_table_true.loc['alignment_length'].values.tolist() == seq_regions_table_new.loc['alignment_length'].values.tolist(), f'Alignment length for senterica is not equal'
>>> assert seq_regions_table_true.loc['ref_start_pos'].values.tolist() == seq_regions_table_new.loc['ref_start_pos'].values.tolist(), f'Ref start for senterica is not equal'
>>> assert seq_regions_table_true.loc['ref_end_pos'].values.tolist() == seq_regions_table_new.loc['ref_end_pos'].values.tolist(), f'Ref end for senterica is not equal'

```

# Test mutations for senterica and mut_type: chromosomal

```
>>> seq_variations_true = json_data_true['seq_variations']
>>> seq_variations_table_true = pd.DataFrame(seq_variations_true)
>>> seq_variations_new = json_data_new['seq_variations']
>>> seq_variations_table_new = pd.DataFrame(seq_variations_new)
>>> assert seq_variations_table_true.shape == seq_variations_table_new.shape, 'Shape of seq_variations is not equal'
>>> if len(seq_variations_table_true) > 0:
...		assert seq_variations_table_new.loc['key'].values.tolist()[0].split(';')[0] == seq_variations_table_true.loc['key'].values.tolist()[0].split(';')[0], 'First element of key is not equal'
...		assert seq_variations_table_new.loc['ref_start_pos'].values.tolist() == seq_variations_table_true.loc['ref_start_pos'].values.tolist(), 'Mutation start position for chromosomal and senterica is not equal'
...		assert seq_variations_table_new.loc['ref_end_pos'].values.tolist() == seq_variations_table_true.loc['ref_end_pos'].values.tolist(), 'Mutation end position for chromosomal and senterica is not equal'
...		assert seq_variations_table_new.loc['var_codon'].values.tolist() == seq_variations_table_true.loc['var_codon'].values.tolist(), 'Mutation codon for chromosomal and senterica is not equal'
...		assert seq_variations_table_new.loc['ref_codon'].values.tolist() == seq_variations_table_true.loc['ref_codon'].values.tolist(), 'Reference codon for chromosomal and senterica is not equal'

```

```
>>> os.remove('/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_01_2.json')

```

## Software Command for ecoli

```
>>> subprocess.run(['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'ecoli', '-acq', '-ifq', '/home/people/s220672/resfinder/tests/data/test_isolate_01_2.fq', '/home/people/s220672/resfinder/tests/data/test_isolate_01_1.fq', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_01_2.json'])
CompletedProcess(args=['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'ecoli', '-acq', '-ifq', '/home/people/s220672/resfinder/tests/data/test_isolate_01_2.fq', '/home/people/s220672/resfinder/tests/data/test_isolate_01_1.fq', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_01_2.json'], returncode=0)
>>> json_file = '/home/people/s220672/functional_testing/functional_testing_data/ecoli/test_isolate_01_2/gjkQh/test_isolate_01_2_-acq.json'
>>> with open(json_file, 'r') as f:
...		json_data_true = json.load(f)
>>> json_file_new = '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_01_2.json'
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

### Test seq_regions for ecoli

```
>>> seq_regions = json_data_true['seq_regions']
>>> seq_regions_table_true = pd.DataFrame(seq_regions)
>>> seq_region_new = json_data_new['seq_regions']
>>> seq_regions_table_new = pd.DataFrame(seq_region_new)
>>> assert seq_regions_table_true.loc['phenotypes'].values.tolist() == seq_regions_table_new.loc['phenotypes'].values.tolist(), 'Phenotypes are not equal'
>>> assert seq_regions_table_true.loc['alignment_length'].values.tolist() == seq_regions_table_new.loc['alignment_length'].values.tolist(), f'Alignment length for ecoli is not equal'
>>> assert seq_regions_table_true.loc['ref_start_pos'].values.tolist() == seq_regions_table_new.loc['ref_start_pos'].values.tolist(), f'Ref start for ecoli is not equal'
>>> assert seq_regions_table_true.loc['ref_end_pos'].values.tolist() == seq_regions_table_new.loc['ref_end_pos'].values.tolist(), f'Ref end for ecoli is not equal'

```

# Test mutations for ecoli and mut_type: acquired

```
>>> seq_variations_true = json_data_true['seq_variations']
>>> seq_variations_table_true = pd.DataFrame(seq_variations_true)
>>> seq_variations_new = json_data_new['seq_variations']
>>> seq_variations_table_new = pd.DataFrame(seq_variations_new)
>>> assert seq_variations_table_true.shape == seq_variations_table_new.shape, 'Shape of seq_variations is not equal'
>>> if len(seq_variations_table_true) > 0:
...		assert seq_variations_table_new.loc['key'].values.tolist()[0].split(';')[0] == seq_variations_table_true.loc['key'].values.tolist()[0].split(';')[0], 'First element of key is not equal'
...		assert seq_variations_table_new.loc['ref_start_pos'].values.tolist() == seq_variations_table_true.loc['ref_start_pos'].values.tolist(), 'Mutation start position for acquired and ecoli is not equal'
...		assert seq_variations_table_new.loc['ref_end_pos'].values.tolist() == seq_variations_table_true.loc['ref_end_pos'].values.tolist(), 'Mutation end position for acquired and ecoli is not equal'
...		assert seq_variations_table_new.loc['var_codon'].values.tolist() == seq_variations_table_true.loc['var_codon'].values.tolist(), 'Mutation codon for acquired and ecoli is not equal'
...		assert seq_variations_table_new.loc['ref_codon'].values.tolist() == seq_variations_table_true.loc['ref_codon'].values.tolist(), 'Reference codon for acquired and ecoli is not equal'

```

```
>>> os.remove('/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_01_2.json')

```

## Software Command for cjejuni

```
>>> subprocess.run(['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'cjejuni', '-acq', '-ifq', '/home/people/s220672/resfinder/tests/data/test_isolate_01_2.fq', '/home/people/s220672/resfinder/tests/data/test_isolate_01_1.fq', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_01_2.json'])
CompletedProcess(args=['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'cjejuni', '-acq', '-ifq', '/home/people/s220672/resfinder/tests/data/test_isolate_01_2.fq', '/home/people/s220672/resfinder/tests/data/test_isolate_01_1.fq', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_01_2.json'], returncode=0)
>>> json_file = '/home/people/s220672/functional_testing/functional_testing_data/cjejuni/test_isolate_01_2/ouPSH/test_isolate_01_2_-acq.json'
>>> with open(json_file, 'r') as f:
...		json_data_true = json.load(f)
>>> json_file_new = '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_01_2.json'
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

### Test seq_regions for cjejuni

```
>>> seq_regions = json_data_true['seq_regions']
>>> seq_regions_table_true = pd.DataFrame(seq_regions)
>>> seq_region_new = json_data_new['seq_regions']
>>> seq_regions_table_new = pd.DataFrame(seq_region_new)
>>> assert seq_regions_table_true.loc['phenotypes'].values.tolist() == seq_regions_table_new.loc['phenotypes'].values.tolist(), 'Phenotypes are not equal'
>>> assert seq_regions_table_true.loc['alignment_length'].values.tolist() == seq_regions_table_new.loc['alignment_length'].values.tolist(), f'Alignment length for cjejuni is not equal'
>>> assert seq_regions_table_true.loc['ref_start_pos'].values.tolist() == seq_regions_table_new.loc['ref_start_pos'].values.tolist(), f'Ref start for cjejuni is not equal'
>>> assert seq_regions_table_true.loc['ref_end_pos'].values.tolist() == seq_regions_table_new.loc['ref_end_pos'].values.tolist(), f'Ref end for cjejuni is not equal'

```

# Test mutations for cjejuni and mut_type: acquired

```
>>> seq_variations_true = json_data_true['seq_variations']
>>> seq_variations_table_true = pd.DataFrame(seq_variations_true)
>>> seq_variations_new = json_data_new['seq_variations']
>>> seq_variations_table_new = pd.DataFrame(seq_variations_new)
>>> assert seq_variations_table_true.shape == seq_variations_table_new.shape, 'Shape of seq_variations is not equal'
>>> if len(seq_variations_table_true) > 0:
...		assert seq_variations_table_new.loc['key'].values.tolist()[0].split(';')[0] == seq_variations_table_true.loc['key'].values.tolist()[0].split(';')[0], 'First element of key is not equal'
...		assert seq_variations_table_new.loc['ref_start_pos'].values.tolist() == seq_variations_table_true.loc['ref_start_pos'].values.tolist(), 'Mutation start position for acquired and cjejuni is not equal'
...		assert seq_variations_table_new.loc['ref_end_pos'].values.tolist() == seq_variations_table_true.loc['ref_end_pos'].values.tolist(), 'Mutation end position for acquired and cjejuni is not equal'
...		assert seq_variations_table_new.loc['var_codon'].values.tolist() == seq_variations_table_true.loc['var_codon'].values.tolist(), 'Mutation codon for acquired and cjejuni is not equal'
...		assert seq_variations_table_new.loc['ref_codon'].values.tolist() == seq_variations_table_true.loc['ref_codon'].values.tolist(), 'Reference codon for acquired and cjejuni is not equal'

```

```
>>> os.remove('/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_01_2.json')

```

## Software Command for ccoli

```
>>> subprocess.run(['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'ccoli', '-acq', '-ifq', '/home/people/s220672/resfinder/tests/data/test_isolate_01_2.fq', '/home/people/s220672/resfinder/tests/data/test_isolate_01_1.fq', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_01_2.json'])
CompletedProcess(args=['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'ccoli', '-acq', '-ifq', '/home/people/s220672/resfinder/tests/data/test_isolate_01_2.fq', '/home/people/s220672/resfinder/tests/data/test_isolate_01_1.fq', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_01_2.json'], returncode=0)
>>> json_file = '/home/people/s220672/functional_testing/functional_testing_data/ccoli/test_isolate_01_2/svyTl/test_isolate_01_2_-acq.json'
>>> with open(json_file, 'r') as f:
...		json_data_true = json.load(f)
>>> json_file_new = '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_01_2.json'
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

### Test seq_regions for ccoli

```
>>> seq_regions = json_data_true['seq_regions']
>>> seq_regions_table_true = pd.DataFrame(seq_regions)
>>> seq_region_new = json_data_new['seq_regions']
>>> seq_regions_table_new = pd.DataFrame(seq_region_new)
>>> assert seq_regions_table_true.loc['phenotypes'].values.tolist() == seq_regions_table_new.loc['phenotypes'].values.tolist(), 'Phenotypes are not equal'
>>> assert seq_regions_table_true.loc['alignment_length'].values.tolist() == seq_regions_table_new.loc['alignment_length'].values.tolist(), f'Alignment length for ccoli is not equal'
>>> assert seq_regions_table_true.loc['ref_start_pos'].values.tolist() == seq_regions_table_new.loc['ref_start_pos'].values.tolist(), f'Ref start for ccoli is not equal'
>>> assert seq_regions_table_true.loc['ref_end_pos'].values.tolist() == seq_regions_table_new.loc['ref_end_pos'].values.tolist(), f'Ref end for ccoli is not equal'

```

# Test mutations for ccoli and mut_type: acquired

```
>>> seq_variations_true = json_data_true['seq_variations']
>>> seq_variations_table_true = pd.DataFrame(seq_variations_true)
>>> seq_variations_new = json_data_new['seq_variations']
>>> seq_variations_table_new = pd.DataFrame(seq_variations_new)
>>> assert seq_variations_table_true.shape == seq_variations_table_new.shape, 'Shape of seq_variations is not equal'
>>> if len(seq_variations_table_true) > 0:
...		assert seq_variations_table_new.loc['key'].values.tolist()[0].split(';')[0] == seq_variations_table_true.loc['key'].values.tolist()[0].split(';')[0], 'First element of key is not equal'
...		assert seq_variations_table_new.loc['ref_start_pos'].values.tolist() == seq_variations_table_true.loc['ref_start_pos'].values.tolist(), 'Mutation start position for acquired and ccoli is not equal'
...		assert seq_variations_table_new.loc['ref_end_pos'].values.tolist() == seq_variations_table_true.loc['ref_end_pos'].values.tolist(), 'Mutation end position for acquired and ccoli is not equal'
...		assert seq_variations_table_new.loc['var_codon'].values.tolist() == seq_variations_table_true.loc['var_codon'].values.tolist(), 'Mutation codon for acquired and ccoli is not equal'
...		assert seq_variations_table_new.loc['ref_codon'].values.tolist() == seq_variations_table_true.loc['ref_codon'].values.tolist(), 'Reference codon for acquired and ccoli is not equal'

```

```
>>> os.remove('/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_01_2.json')

```

## Software Command for senterica

```
>>> subprocess.run(['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'senterica', '-acq', '-ifq', '/home/people/s220672/resfinder/tests/data/test_isolate_01_2.fq', '/home/people/s220672/resfinder/tests/data/test_isolate_01_1.fq', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_01_2.json'])
CompletedProcess(args=['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'senterica', '-acq', '-ifq', '/home/people/s220672/resfinder/tests/data/test_isolate_01_2.fq', '/home/people/s220672/resfinder/tests/data/test_isolate_01_1.fq', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_01_2.json'], returncode=0)
>>> json_file = '/home/people/s220672/functional_testing/functional_testing_data/senterica/test_isolate_01_2/gOeXe/test_isolate_01_2_-acq.json'
>>> with open(json_file, 'r') as f:
...		json_data_true = json.load(f)
>>> json_file_new = '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_01_2.json'
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

### Test seq_regions for senterica

```
>>> seq_regions = json_data_true['seq_regions']
>>> seq_regions_table_true = pd.DataFrame(seq_regions)
>>> seq_region_new = json_data_new['seq_regions']
>>> seq_regions_table_new = pd.DataFrame(seq_region_new)
>>> assert seq_regions_table_true.loc['phenotypes'].values.tolist() == seq_regions_table_new.loc['phenotypes'].values.tolist(), 'Phenotypes are not equal'
>>> assert seq_regions_table_true.loc['alignment_length'].values.tolist() == seq_regions_table_new.loc['alignment_length'].values.tolist(), f'Alignment length for senterica is not equal'
>>> assert seq_regions_table_true.loc['ref_start_pos'].values.tolist() == seq_regions_table_new.loc['ref_start_pos'].values.tolist(), f'Ref start for senterica is not equal'
>>> assert seq_regions_table_true.loc['ref_end_pos'].values.tolist() == seq_regions_table_new.loc['ref_end_pos'].values.tolist(), f'Ref end for senterica is not equal'

```

# Test mutations for senterica and mut_type: acquired

```
>>> seq_variations_true = json_data_true['seq_variations']
>>> seq_variations_table_true = pd.DataFrame(seq_variations_true)
>>> seq_variations_new = json_data_new['seq_variations']
>>> seq_variations_table_new = pd.DataFrame(seq_variations_new)
>>> assert seq_variations_table_true.shape == seq_variations_table_new.shape, 'Shape of seq_variations is not equal'
>>> if len(seq_variations_table_true) > 0:
...		assert seq_variations_table_new.loc['key'].values.tolist()[0].split(';')[0] == seq_variations_table_true.loc['key'].values.tolist()[0].split(';')[0], 'First element of key is not equal'
...		assert seq_variations_table_new.loc['ref_start_pos'].values.tolist() == seq_variations_table_true.loc['ref_start_pos'].values.tolist(), 'Mutation start position for acquired and senterica is not equal'
...		assert seq_variations_table_new.loc['ref_end_pos'].values.tolist() == seq_variations_table_true.loc['ref_end_pos'].values.tolist(), 'Mutation end position for acquired and senterica is not equal'
...		assert seq_variations_table_new.loc['var_codon'].values.tolist() == seq_variations_table_true.loc['var_codon'].values.tolist(), 'Mutation codon for acquired and senterica is not equal'
...		assert seq_variations_table_new.loc['ref_codon'].values.tolist() == seq_variations_table_true.loc['ref_codon'].values.tolist(), 'Reference codon for acquired and senterica is not equal'

```

```
>>> os.remove('/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_01_2.json')

```

## Software Command for ecoli

```
>>> subprocess.run(['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'ecoli', '-u', '-acq', '-ifq', '/home/people/s220672/resfinder/tests/data/test_isolate_01_2.fq', '/home/people/s220672/resfinder/tests/data/test_isolate_01_1.fq', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_01_2.json'])
CompletedProcess(args=['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'ecoli', '-u', '-acq', '-ifq', '/home/people/s220672/resfinder/tests/data/test_isolate_01_2.fq', '/home/people/s220672/resfinder/tests/data/test_isolate_01_1.fq', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_01_2.json'], returncode=0)
>>> json_file = '/home/people/s220672/functional_testing/functional_testing_data/ecoli/test_isolate_01_2/tFgpZ/test_isolate_01_2_-acq-u.json'
>>> with open(json_file, 'r') as f:
...		json_data_true = json.load(f)
>>> json_file_new = '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_01_2.json'
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

### Test seq_regions for ecoli

```
>>> seq_regions = json_data_true['seq_regions']
>>> seq_regions_table_true = pd.DataFrame(seq_regions)
>>> seq_region_new = json_data_new['seq_regions']
>>> seq_regions_table_new = pd.DataFrame(seq_region_new)
>>> assert seq_regions_table_true.loc['phenotypes'].values.tolist() == seq_regions_table_new.loc['phenotypes'].values.tolist(), 'Phenotypes are not equal'
>>> assert seq_regions_table_true.loc['alignment_length'].values.tolist() == seq_regions_table_new.loc['alignment_length'].values.tolist(), f'Alignment length for ecoli is not equal'
>>> assert seq_regions_table_true.loc['ref_start_pos'].values.tolist() == seq_regions_table_new.loc['ref_start_pos'].values.tolist(), f'Ref start for ecoli is not equal'
>>> assert seq_regions_table_true.loc['ref_end_pos'].values.tolist() == seq_regions_table_new.loc['ref_end_pos'].values.tolist(), f'Ref end for ecoli is not equal'

```

# Test mutations for ecoli and mut_type: unknown

```
>>> seq_variations_true = json_data_true['seq_variations']
>>> seq_variations_table_true = pd.DataFrame(seq_variations_true)
>>> seq_variations_new = json_data_new['seq_variations']
>>> seq_variations_table_new = pd.DataFrame(seq_variations_new)
>>> assert seq_variations_table_true.shape == seq_variations_table_new.shape, 'Shape of seq_variations is not equal'
>>> if len(seq_variations_table_true) > 0:
...		assert seq_variations_table_new.loc['key'].values.tolist()[0].split(';')[0] == seq_variations_table_true.loc['key'].values.tolist()[0].split(';')[0], 'First element of key is not equal'
...		assert seq_variations_table_new.loc['ref_start_pos'].values.tolist() == seq_variations_table_true.loc['ref_start_pos'].values.tolist(), 'Mutation start position for unknown and ecoli is not equal'
...		assert seq_variations_table_new.loc['ref_end_pos'].values.tolist() == seq_variations_table_true.loc['ref_end_pos'].values.tolist(), 'Mutation end position for unknown and ecoli is not equal'
...		assert seq_variations_table_new.loc['var_codon'].values.tolist() == seq_variations_table_true.loc['var_codon'].values.tolist(), 'Mutation codon for unknown and ecoli is not equal'
...		assert seq_variations_table_new.loc['ref_codon'].values.tolist() == seq_variations_table_true.loc['ref_codon'].values.tolist(), 'Reference codon for unknown and ecoli is not equal'

```

```
>>> os.remove('/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_01_2.json')

```

## Software Command for cjejuni

```
>>> subprocess.run(['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'cjejuni', '-u', '-acq', '-ifq', '/home/people/s220672/resfinder/tests/data/test_isolate_01_2.fq', '/home/people/s220672/resfinder/tests/data/test_isolate_01_1.fq', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_01_2.json'])
CompletedProcess(args=['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'cjejuni', '-u', '-acq', '-ifq', '/home/people/s220672/resfinder/tests/data/test_isolate_01_2.fq', '/home/people/s220672/resfinder/tests/data/test_isolate_01_1.fq', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_01_2.json'], returncode=0)
>>> json_file = '/home/people/s220672/functional_testing/functional_testing_data/cjejuni/test_isolate_01_2/lIfif/test_isolate_01_2_-acq-u.json'
>>> with open(json_file, 'r') as f:
...		json_data_true = json.load(f)
>>> json_file_new = '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_01_2.json'
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

### Test seq_regions for cjejuni

```
>>> seq_regions = json_data_true['seq_regions']
>>> seq_regions_table_true = pd.DataFrame(seq_regions)
>>> seq_region_new = json_data_new['seq_regions']
>>> seq_regions_table_new = pd.DataFrame(seq_region_new)
>>> assert seq_regions_table_true.loc['phenotypes'].values.tolist() == seq_regions_table_new.loc['phenotypes'].values.tolist(), 'Phenotypes are not equal'
>>> assert seq_regions_table_true.loc['alignment_length'].values.tolist() == seq_regions_table_new.loc['alignment_length'].values.tolist(), f'Alignment length for cjejuni is not equal'
>>> assert seq_regions_table_true.loc['ref_start_pos'].values.tolist() == seq_regions_table_new.loc['ref_start_pos'].values.tolist(), f'Ref start for cjejuni is not equal'
>>> assert seq_regions_table_true.loc['ref_end_pos'].values.tolist() == seq_regions_table_new.loc['ref_end_pos'].values.tolist(), f'Ref end for cjejuni is not equal'

```

# Test mutations for cjejuni and mut_type: unknown

```
>>> seq_variations_true = json_data_true['seq_variations']
>>> seq_variations_table_true = pd.DataFrame(seq_variations_true)
>>> seq_variations_new = json_data_new['seq_variations']
>>> seq_variations_table_new = pd.DataFrame(seq_variations_new)
>>> assert seq_variations_table_true.shape == seq_variations_table_new.shape, 'Shape of seq_variations is not equal'
>>> if len(seq_variations_table_true) > 0:
...		assert seq_variations_table_new.loc['key'].values.tolist()[0].split(';')[0] == seq_variations_table_true.loc['key'].values.tolist()[0].split(';')[0], 'First element of key is not equal'
...		assert seq_variations_table_new.loc['ref_start_pos'].values.tolist() == seq_variations_table_true.loc['ref_start_pos'].values.tolist(), 'Mutation start position for unknown and cjejuni is not equal'
...		assert seq_variations_table_new.loc['ref_end_pos'].values.tolist() == seq_variations_table_true.loc['ref_end_pos'].values.tolist(), 'Mutation end position for unknown and cjejuni is not equal'
...		assert seq_variations_table_new.loc['var_codon'].values.tolist() == seq_variations_table_true.loc['var_codon'].values.tolist(), 'Mutation codon for unknown and cjejuni is not equal'
...		assert seq_variations_table_new.loc['ref_codon'].values.tolist() == seq_variations_table_true.loc['ref_codon'].values.tolist(), 'Reference codon for unknown and cjejuni is not equal'

```

```
>>> os.remove('/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_01_2.json')

```

## Software Command for ccoli

```
>>> subprocess.run(['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'ccoli', '-u', '-acq', '-ifq', '/home/people/s220672/resfinder/tests/data/test_isolate_01_2.fq', '/home/people/s220672/resfinder/tests/data/test_isolate_01_1.fq', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_01_2.json'])
CompletedProcess(args=['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'ccoli', '-u', '-acq', '-ifq', '/home/people/s220672/resfinder/tests/data/test_isolate_01_2.fq', '/home/people/s220672/resfinder/tests/data/test_isolate_01_1.fq', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_01_2.json'], returncode=0)
>>> json_file = '/home/people/s220672/functional_testing/functional_testing_data/ccoli/test_isolate_01_2/GrCpa/test_isolate_01_2_-acq-u.json'
>>> with open(json_file, 'r') as f:
...		json_data_true = json.load(f)
>>> json_file_new = '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_01_2.json'
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

### Test seq_regions for ccoli

```
>>> seq_regions = json_data_true['seq_regions']
>>> seq_regions_table_true = pd.DataFrame(seq_regions)
>>> seq_region_new = json_data_new['seq_regions']
>>> seq_regions_table_new = pd.DataFrame(seq_region_new)
>>> assert seq_regions_table_true.loc['phenotypes'].values.tolist() == seq_regions_table_new.loc['phenotypes'].values.tolist(), 'Phenotypes are not equal'
>>> assert seq_regions_table_true.loc['alignment_length'].values.tolist() == seq_regions_table_new.loc['alignment_length'].values.tolist(), f'Alignment length for ccoli is not equal'
>>> assert seq_regions_table_true.loc['ref_start_pos'].values.tolist() == seq_regions_table_new.loc['ref_start_pos'].values.tolist(), f'Ref start for ccoli is not equal'
>>> assert seq_regions_table_true.loc['ref_end_pos'].values.tolist() == seq_regions_table_new.loc['ref_end_pos'].values.tolist(), f'Ref end for ccoli is not equal'

```

# Test mutations for ccoli and mut_type: unknown

```
>>> seq_variations_true = json_data_true['seq_variations']
>>> seq_variations_table_true = pd.DataFrame(seq_variations_true)
>>> seq_variations_new = json_data_new['seq_variations']
>>> seq_variations_table_new = pd.DataFrame(seq_variations_new)
>>> assert seq_variations_table_true.shape == seq_variations_table_new.shape, 'Shape of seq_variations is not equal'
>>> if len(seq_variations_table_true) > 0:
...		assert seq_variations_table_new.loc['key'].values.tolist()[0].split(';')[0] == seq_variations_table_true.loc['key'].values.tolist()[0].split(';')[0], 'First element of key is not equal'
...		assert seq_variations_table_new.loc['ref_start_pos'].values.tolist() == seq_variations_table_true.loc['ref_start_pos'].values.tolist(), 'Mutation start position for unknown and ccoli is not equal'
...		assert seq_variations_table_new.loc['ref_end_pos'].values.tolist() == seq_variations_table_true.loc['ref_end_pos'].values.tolist(), 'Mutation end position for unknown and ccoli is not equal'
...		assert seq_variations_table_new.loc['var_codon'].values.tolist() == seq_variations_table_true.loc['var_codon'].values.tolist(), 'Mutation codon for unknown and ccoli is not equal'
...		assert seq_variations_table_new.loc['ref_codon'].values.tolist() == seq_variations_table_true.loc['ref_codon'].values.tolist(), 'Reference codon for unknown and ccoli is not equal'

```

```
>>> os.remove('/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_01_2.json')

```

## Software Command for senterica

```
>>> subprocess.run(['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'senterica', '-u', '-acq', '-ifq', '/home/people/s220672/resfinder/tests/data/test_isolate_01_2.fq', '/home/people/s220672/resfinder/tests/data/test_isolate_01_1.fq', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_01_2.json'])
CompletedProcess(args=['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'senterica', '-u', '-acq', '-ifq', '/home/people/s220672/resfinder/tests/data/test_isolate_01_2.fq', '/home/people/s220672/resfinder/tests/data/test_isolate_01_1.fq', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_01_2.json'], returncode=0)
>>> json_file = '/home/people/s220672/functional_testing/functional_testing_data/senterica/test_isolate_01_2/Dqube/test_isolate_01_2_-acq-u.json'
>>> with open(json_file, 'r') as f:
...		json_data_true = json.load(f)
>>> json_file_new = '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_01_2.json'
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

### Test seq_regions for senterica

```
>>> seq_regions = json_data_true['seq_regions']
>>> seq_regions_table_true = pd.DataFrame(seq_regions)
>>> seq_region_new = json_data_new['seq_regions']
>>> seq_regions_table_new = pd.DataFrame(seq_region_new)
>>> assert seq_regions_table_true.loc['phenotypes'].values.tolist() == seq_regions_table_new.loc['phenotypes'].values.tolist(), 'Phenotypes are not equal'
>>> assert seq_regions_table_true.loc['alignment_length'].values.tolist() == seq_regions_table_new.loc['alignment_length'].values.tolist(), f'Alignment length for senterica is not equal'
>>> assert seq_regions_table_true.loc['ref_start_pos'].values.tolist() == seq_regions_table_new.loc['ref_start_pos'].values.tolist(), f'Ref start for senterica is not equal'
>>> assert seq_regions_table_true.loc['ref_end_pos'].values.tolist() == seq_regions_table_new.loc['ref_end_pos'].values.tolist(), f'Ref end for senterica is not equal'

```

# Test mutations for senterica and mut_type: unknown

```
>>> seq_variations_true = json_data_true['seq_variations']
>>> seq_variations_table_true = pd.DataFrame(seq_variations_true)
>>> seq_variations_new = json_data_new['seq_variations']
>>> seq_variations_table_new = pd.DataFrame(seq_variations_new)
>>> assert seq_variations_table_true.shape == seq_variations_table_new.shape, 'Shape of seq_variations is not equal'
>>> if len(seq_variations_table_true) > 0:
...		assert seq_variations_table_new.loc['key'].values.tolist()[0].split(';')[0] == seq_variations_table_true.loc['key'].values.tolist()[0].split(';')[0], 'First element of key is not equal'
...		assert seq_variations_table_new.loc['ref_start_pos'].values.tolist() == seq_variations_table_true.loc['ref_start_pos'].values.tolist(), 'Mutation start position for unknown and senterica is not equal'
...		assert seq_variations_table_new.loc['ref_end_pos'].values.tolist() == seq_variations_table_true.loc['ref_end_pos'].values.tolist(), 'Mutation end position for unknown and senterica is not equal'
...		assert seq_variations_table_new.loc['var_codon'].values.tolist() == seq_variations_table_true.loc['var_codon'].values.tolist(), 'Mutation codon for unknown and senterica is not equal'
...		assert seq_variations_table_new.loc['ref_codon'].values.tolist() == seq_variations_table_true.loc['ref_codon'].values.tolist(), 'Reference codon for unknown and senterica is not equal'

```

```
>>> os.remove('/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_01_2.json')

```

## Software Command for ecoli

```
>>> subprocess.run(['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'ecoli', '-acq', '-ifa', '/home/people/s220672/resfinder/tests/data/test_isolate_02.fa', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_02.json'])
CompletedProcess(args=['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'ecoli', '-acq', '-ifa', '/home/people/s220672/resfinder/tests/data/test_isolate_02.fa', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_02.json'], returncode=0)
>>> json_file = '/home/people/s220672/functional_testing/functional_testing_data/ecoli/test_isolate_02/JXQsM/test_isolate_02_-acq.json'
>>> with open(json_file, 'r') as f:
...		json_data_true = json.load(f)
>>> json_file_new = '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_02.json'
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

### Test seq_regions for ecoli

```
>>> seq_regions = json_data_true['seq_regions']
>>> seq_regions_table_true = pd.DataFrame(seq_regions)
>>> seq_region_new = json_data_new['seq_regions']
>>> seq_regions_table_new = pd.DataFrame(seq_region_new)
>>> assert seq_regions_table_true.loc['phenotypes'].values.tolist() == seq_regions_table_new.loc['phenotypes'].values.tolist(), 'Phenotypes are not equal'
>>> assert seq_regions_table_true.loc['alignment_length'].values.tolist() == seq_regions_table_new.loc['alignment_length'].values.tolist(), f'Alignment length for ecoli is not equal'
>>> assert seq_regions_table_true.loc['ref_start_pos'].values.tolist() == seq_regions_table_new.loc['ref_start_pos'].values.tolist(), f'Ref start for ecoli is not equal'
>>> assert seq_regions_table_true.loc['ref_end_pos'].values.tolist() == seq_regions_table_new.loc['ref_end_pos'].values.tolist(), f'Ref end for ecoli is not equal'

```

# Test mutations for ecoli and mut_type: acquired

```
>>> seq_variations_true = json_data_true['seq_variations']
>>> seq_variations_table_true = pd.DataFrame(seq_variations_true)
>>> seq_variations_new = json_data_new['seq_variations']
>>> seq_variations_table_new = pd.DataFrame(seq_variations_new)
>>> assert seq_variations_table_true.shape == seq_variations_table_new.shape, 'Shape of seq_variations is not equal'
>>> if len(seq_variations_table_true) > 0:
...		assert seq_variations_table_new.loc['key'].values.tolist()[0].split(';')[0] == seq_variations_table_true.loc['key'].values.tolist()[0].split(';')[0], 'First element of key is not equal'
...		assert seq_variations_table_new.loc['ref_start_pos'].values.tolist() == seq_variations_table_true.loc['ref_start_pos'].values.tolist(), 'Mutation start position for acquired and ecoli is not equal'
...		assert seq_variations_table_new.loc['ref_end_pos'].values.tolist() == seq_variations_table_true.loc['ref_end_pos'].values.tolist(), 'Mutation end position for acquired and ecoli is not equal'
...		assert seq_variations_table_new.loc['var_codon'].values.tolist() == seq_variations_table_true.loc['var_codon'].values.tolist(), 'Mutation codon for acquired and ecoli is not equal'
...		assert seq_variations_table_new.loc['ref_codon'].values.tolist() == seq_variations_table_true.loc['ref_codon'].values.tolist(), 'Reference codon for acquired and ecoli is not equal'

```

```
>>> os.remove('/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_02.json')

```

## Software Command for cjejuni

```
>>> subprocess.run(['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'cjejuni', '-acq', '-ifa', '/home/people/s220672/resfinder/tests/data/test_isolate_02.fa', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_02.json'])
CompletedProcess(args=['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'cjejuni', '-acq', '-ifa', '/home/people/s220672/resfinder/tests/data/test_isolate_02.fa', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_02.json'], returncode=0)
>>> json_file = '/home/people/s220672/functional_testing/functional_testing_data/cjejuni/test_isolate_02/KCYoi/test_isolate_02_-acq.json'
>>> with open(json_file, 'r') as f:
...		json_data_true = json.load(f)
>>> json_file_new = '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_02.json'
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

### Test seq_regions for cjejuni

```
>>> seq_regions = json_data_true['seq_regions']
>>> seq_regions_table_true = pd.DataFrame(seq_regions)
>>> seq_region_new = json_data_new['seq_regions']
>>> seq_regions_table_new = pd.DataFrame(seq_region_new)
>>> assert seq_regions_table_true.loc['phenotypes'].values.tolist() == seq_regions_table_new.loc['phenotypes'].values.tolist(), 'Phenotypes are not equal'
>>> assert seq_regions_table_true.loc['alignment_length'].values.tolist() == seq_regions_table_new.loc['alignment_length'].values.tolist(), f'Alignment length for cjejuni is not equal'
>>> assert seq_regions_table_true.loc['ref_start_pos'].values.tolist() == seq_regions_table_new.loc['ref_start_pos'].values.tolist(), f'Ref start for cjejuni is not equal'
>>> assert seq_regions_table_true.loc['ref_end_pos'].values.tolist() == seq_regions_table_new.loc['ref_end_pos'].values.tolist(), f'Ref end for cjejuni is not equal'

```

# Test mutations for cjejuni and mut_type: acquired

```
>>> seq_variations_true = json_data_true['seq_variations']
>>> seq_variations_table_true = pd.DataFrame(seq_variations_true)
>>> seq_variations_new = json_data_new['seq_variations']
>>> seq_variations_table_new = pd.DataFrame(seq_variations_new)
>>> assert seq_variations_table_true.shape == seq_variations_table_new.shape, 'Shape of seq_variations is not equal'
>>> if len(seq_variations_table_true) > 0:
...		assert seq_variations_table_new.loc['key'].values.tolist()[0].split(';')[0] == seq_variations_table_true.loc['key'].values.tolist()[0].split(';')[0], 'First element of key is not equal'
...		assert seq_variations_table_new.loc['ref_start_pos'].values.tolist() == seq_variations_table_true.loc['ref_start_pos'].values.tolist(), 'Mutation start position for acquired and cjejuni is not equal'
...		assert seq_variations_table_new.loc['ref_end_pos'].values.tolist() == seq_variations_table_true.loc['ref_end_pos'].values.tolist(), 'Mutation end position for acquired and cjejuni is not equal'
...		assert seq_variations_table_new.loc['var_codon'].values.tolist() == seq_variations_table_true.loc['var_codon'].values.tolist(), 'Mutation codon for acquired and cjejuni is not equal'
...		assert seq_variations_table_new.loc['ref_codon'].values.tolist() == seq_variations_table_true.loc['ref_codon'].values.tolist(), 'Reference codon for acquired and cjejuni is not equal'

```

```
>>> os.remove('/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_02.json')

```

## Software Command for ccoli

```
>>> subprocess.run(['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'ccoli', '-acq', '-ifa', '/home/people/s220672/resfinder/tests/data/test_isolate_02.fa', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_02.json'])
CompletedProcess(args=['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'ccoli', '-acq', '-ifa', '/home/people/s220672/resfinder/tests/data/test_isolate_02.fa', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_02.json'], returncode=0)
>>> json_file = '/home/people/s220672/functional_testing/functional_testing_data/ccoli/test_isolate_02/yXUzi/test_isolate_02_-acq.json'
>>> with open(json_file, 'r') as f:
...		json_data_true = json.load(f)
>>> json_file_new = '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_02.json'
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

### Test seq_regions for ccoli

```
>>> seq_regions = json_data_true['seq_regions']
>>> seq_regions_table_true = pd.DataFrame(seq_regions)
>>> seq_region_new = json_data_new['seq_regions']
>>> seq_regions_table_new = pd.DataFrame(seq_region_new)
>>> assert seq_regions_table_true.loc['phenotypes'].values.tolist() == seq_regions_table_new.loc['phenotypes'].values.tolist(), 'Phenotypes are not equal'
>>> assert seq_regions_table_true.loc['alignment_length'].values.tolist() == seq_regions_table_new.loc['alignment_length'].values.tolist(), f'Alignment length for ccoli is not equal'
>>> assert seq_regions_table_true.loc['ref_start_pos'].values.tolist() == seq_regions_table_new.loc['ref_start_pos'].values.tolist(), f'Ref start for ccoli is not equal'
>>> assert seq_regions_table_true.loc['ref_end_pos'].values.tolist() == seq_regions_table_new.loc['ref_end_pos'].values.tolist(), f'Ref end for ccoli is not equal'

```

# Test mutations for ccoli and mut_type: acquired

```
>>> seq_variations_true = json_data_true['seq_variations']
>>> seq_variations_table_true = pd.DataFrame(seq_variations_true)
>>> seq_variations_new = json_data_new['seq_variations']
>>> seq_variations_table_new = pd.DataFrame(seq_variations_new)
>>> assert seq_variations_table_true.shape == seq_variations_table_new.shape, 'Shape of seq_variations is not equal'
>>> if len(seq_variations_table_true) > 0:
...		assert seq_variations_table_new.loc['key'].values.tolist()[0].split(';')[0] == seq_variations_table_true.loc['key'].values.tolist()[0].split(';')[0], 'First element of key is not equal'
...		assert seq_variations_table_new.loc['ref_start_pos'].values.tolist() == seq_variations_table_true.loc['ref_start_pos'].values.tolist(), 'Mutation start position for acquired and ccoli is not equal'
...		assert seq_variations_table_new.loc['ref_end_pos'].values.tolist() == seq_variations_table_true.loc['ref_end_pos'].values.tolist(), 'Mutation end position for acquired and ccoli is not equal'
...		assert seq_variations_table_new.loc['var_codon'].values.tolist() == seq_variations_table_true.loc['var_codon'].values.tolist(), 'Mutation codon for acquired and ccoli is not equal'
...		assert seq_variations_table_new.loc['ref_codon'].values.tolist() == seq_variations_table_true.loc['ref_codon'].values.tolist(), 'Reference codon for acquired and ccoli is not equal'

```

```
>>> os.remove('/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_02.json')

```

## Software Command for senterica

```
>>> subprocess.run(['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'senterica', '-acq', '-ifa', '/home/people/s220672/resfinder/tests/data/test_isolate_02.fa', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_02.json'])
CompletedProcess(args=['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'senterica', '-acq', '-ifa', '/home/people/s220672/resfinder/tests/data/test_isolate_02.fa', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_02.json'], returncode=0)
>>> json_file = '/home/people/s220672/functional_testing/functional_testing_data/senterica/test_isolate_02/OoWYL/test_isolate_02_-acq.json'
>>> with open(json_file, 'r') as f:
...		json_data_true = json.load(f)
>>> json_file_new = '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_02.json'
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

### Test seq_regions for senterica

```
>>> seq_regions = json_data_true['seq_regions']
>>> seq_regions_table_true = pd.DataFrame(seq_regions)
>>> seq_region_new = json_data_new['seq_regions']
>>> seq_regions_table_new = pd.DataFrame(seq_region_new)
>>> assert seq_regions_table_true.loc['phenotypes'].values.tolist() == seq_regions_table_new.loc['phenotypes'].values.tolist(), 'Phenotypes are not equal'
>>> assert seq_regions_table_true.loc['alignment_length'].values.tolist() == seq_regions_table_new.loc['alignment_length'].values.tolist(), f'Alignment length for senterica is not equal'
>>> assert seq_regions_table_true.loc['ref_start_pos'].values.tolist() == seq_regions_table_new.loc['ref_start_pos'].values.tolist(), f'Ref start for senterica is not equal'
>>> assert seq_regions_table_true.loc['ref_end_pos'].values.tolist() == seq_regions_table_new.loc['ref_end_pos'].values.tolist(), f'Ref end for senterica is not equal'

```

# Test mutations for senterica and mut_type: acquired

```
>>> seq_variations_true = json_data_true['seq_variations']
>>> seq_variations_table_true = pd.DataFrame(seq_variations_true)
>>> seq_variations_new = json_data_new['seq_variations']
>>> seq_variations_table_new = pd.DataFrame(seq_variations_new)
>>> assert seq_variations_table_true.shape == seq_variations_table_new.shape, 'Shape of seq_variations is not equal'
>>> if len(seq_variations_table_true) > 0:
...		assert seq_variations_table_new.loc['key'].values.tolist()[0].split(';')[0] == seq_variations_table_true.loc['key'].values.tolist()[0].split(';')[0], 'First element of key is not equal'
...		assert seq_variations_table_new.loc['ref_start_pos'].values.tolist() == seq_variations_table_true.loc['ref_start_pos'].values.tolist(), 'Mutation start position for acquired and senterica is not equal'
...		assert seq_variations_table_new.loc['ref_end_pos'].values.tolist() == seq_variations_table_true.loc['ref_end_pos'].values.tolist(), 'Mutation end position for acquired and senterica is not equal'
...		assert seq_variations_table_new.loc['var_codon'].values.tolist() == seq_variations_table_true.loc['var_codon'].values.tolist(), 'Mutation codon for acquired and senterica is not equal'
...		assert seq_variations_table_new.loc['ref_codon'].values.tolist() == seq_variations_table_true.loc['ref_codon'].values.tolist(), 'Reference codon for acquired and senterica is not equal'

```

```
>>> os.remove('/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_02.json')

```

## Software Command for ecoli

```
>>> subprocess.run(['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'ecoli', '-u', '-acq', '-ifa', '/home/people/s220672/resfinder/tests/data/test_isolate_02.fa', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_02.json'])
CompletedProcess(args=['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'ecoli', '-u', '-acq', '-ifa', '/home/people/s220672/resfinder/tests/data/test_isolate_02.fa', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_02.json'], returncode=0)
>>> json_file = '/home/people/s220672/functional_testing/functional_testing_data/ecoli/test_isolate_02/LpnNp/test_isolate_02_-acq-u.json'
>>> with open(json_file, 'r') as f:
...		json_data_true = json.load(f)
>>> json_file_new = '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_02.json'
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

### Test seq_regions for ecoli

```
>>> seq_regions = json_data_true['seq_regions']
>>> seq_regions_table_true = pd.DataFrame(seq_regions)
>>> seq_region_new = json_data_new['seq_regions']
>>> seq_regions_table_new = pd.DataFrame(seq_region_new)
>>> assert seq_regions_table_true.loc['phenotypes'].values.tolist() == seq_regions_table_new.loc['phenotypes'].values.tolist(), 'Phenotypes are not equal'
>>> assert seq_regions_table_true.loc['alignment_length'].values.tolist() == seq_regions_table_new.loc['alignment_length'].values.tolist(), f'Alignment length for ecoli is not equal'
>>> assert seq_regions_table_true.loc['ref_start_pos'].values.tolist() == seq_regions_table_new.loc['ref_start_pos'].values.tolist(), f'Ref start for ecoli is not equal'
>>> assert seq_regions_table_true.loc['ref_end_pos'].values.tolist() == seq_regions_table_new.loc['ref_end_pos'].values.tolist(), f'Ref end for ecoli is not equal'

```

# Test mutations for ecoli and mut_type: unknown

```
>>> seq_variations_true = json_data_true['seq_variations']
>>> seq_variations_table_true = pd.DataFrame(seq_variations_true)
>>> seq_variations_new = json_data_new['seq_variations']
>>> seq_variations_table_new = pd.DataFrame(seq_variations_new)
>>> assert seq_variations_table_true.shape == seq_variations_table_new.shape, 'Shape of seq_variations is not equal'
>>> if len(seq_variations_table_true) > 0:
...		assert seq_variations_table_new.loc['key'].values.tolist()[0].split(';')[0] == seq_variations_table_true.loc['key'].values.tolist()[0].split(';')[0], 'First element of key is not equal'
...		assert seq_variations_table_new.loc['ref_start_pos'].values.tolist() == seq_variations_table_true.loc['ref_start_pos'].values.tolist(), 'Mutation start position for unknown and ecoli is not equal'
...		assert seq_variations_table_new.loc['ref_end_pos'].values.tolist() == seq_variations_table_true.loc['ref_end_pos'].values.tolist(), 'Mutation end position for unknown and ecoli is not equal'
...		assert seq_variations_table_new.loc['var_codon'].values.tolist() == seq_variations_table_true.loc['var_codon'].values.tolist(), 'Mutation codon for unknown and ecoli is not equal'
...		assert seq_variations_table_new.loc['ref_codon'].values.tolist() == seq_variations_table_true.loc['ref_codon'].values.tolist(), 'Reference codon for unknown and ecoli is not equal'

```

```
>>> os.remove('/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_02.json')

```

## Software Command for cjejuni

```
>>> subprocess.run(['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'cjejuni', '-u', '-acq', '-ifa', '/home/people/s220672/resfinder/tests/data/test_isolate_02.fa', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_02.json'])
CompletedProcess(args=['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'cjejuni', '-u', '-acq', '-ifa', '/home/people/s220672/resfinder/tests/data/test_isolate_02.fa', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_02.json'], returncode=0)
>>> json_file = '/home/people/s220672/functional_testing/functional_testing_data/cjejuni/test_isolate_02/PMcNA/test_isolate_02_-acq-u.json'
>>> with open(json_file, 'r') as f:
...		json_data_true = json.load(f)
>>> json_file_new = '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_02.json'
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

### Test seq_regions for cjejuni

```
>>> seq_regions = json_data_true['seq_regions']
>>> seq_regions_table_true = pd.DataFrame(seq_regions)
>>> seq_region_new = json_data_new['seq_regions']
>>> seq_regions_table_new = pd.DataFrame(seq_region_new)
>>> assert seq_regions_table_true.loc['phenotypes'].values.tolist() == seq_regions_table_new.loc['phenotypes'].values.tolist(), 'Phenotypes are not equal'
>>> assert seq_regions_table_true.loc['alignment_length'].values.tolist() == seq_regions_table_new.loc['alignment_length'].values.tolist(), f'Alignment length for cjejuni is not equal'
>>> assert seq_regions_table_true.loc['ref_start_pos'].values.tolist() == seq_regions_table_new.loc['ref_start_pos'].values.tolist(), f'Ref start for cjejuni is not equal'
>>> assert seq_regions_table_true.loc['ref_end_pos'].values.tolist() == seq_regions_table_new.loc['ref_end_pos'].values.tolist(), f'Ref end for cjejuni is not equal'

```

# Test mutations for cjejuni and mut_type: unknown

```
>>> seq_variations_true = json_data_true['seq_variations']
>>> seq_variations_table_true = pd.DataFrame(seq_variations_true)
>>> seq_variations_new = json_data_new['seq_variations']
>>> seq_variations_table_new = pd.DataFrame(seq_variations_new)
>>> assert seq_variations_table_true.shape == seq_variations_table_new.shape, 'Shape of seq_variations is not equal'
>>> if len(seq_variations_table_true) > 0:
...		assert seq_variations_table_new.loc['key'].values.tolist()[0].split(';')[0] == seq_variations_table_true.loc['key'].values.tolist()[0].split(';')[0], 'First element of key is not equal'
...		assert seq_variations_table_new.loc['ref_start_pos'].values.tolist() == seq_variations_table_true.loc['ref_start_pos'].values.tolist(), 'Mutation start position for unknown and cjejuni is not equal'
...		assert seq_variations_table_new.loc['ref_end_pos'].values.tolist() == seq_variations_table_true.loc['ref_end_pos'].values.tolist(), 'Mutation end position for unknown and cjejuni is not equal'
...		assert seq_variations_table_new.loc['var_codon'].values.tolist() == seq_variations_table_true.loc['var_codon'].values.tolist(), 'Mutation codon for unknown and cjejuni is not equal'
...		assert seq_variations_table_new.loc['ref_codon'].values.tolist() == seq_variations_table_true.loc['ref_codon'].values.tolist(), 'Reference codon for unknown and cjejuni is not equal'

```

```
>>> os.remove('/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_02.json')

```

## Software Command for ccoli

```
>>> subprocess.run(['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'ccoli', '-u', '-acq', '-ifa', '/home/people/s220672/resfinder/tests/data/test_isolate_02.fa', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_02.json'])
CompletedProcess(args=['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'ccoli', '-u', '-acq', '-ifa', '/home/people/s220672/resfinder/tests/data/test_isolate_02.fa', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_02.json'], returncode=0)
>>> json_file = '/home/people/s220672/functional_testing/functional_testing_data/ccoli/test_isolate_02/JYVqr/test_isolate_02_-acq-u.json'
>>> with open(json_file, 'r') as f:
...		json_data_true = json.load(f)
>>> json_file_new = '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_02.json'
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

### Test seq_regions for ccoli

```
>>> seq_regions = json_data_true['seq_regions']
>>> seq_regions_table_true = pd.DataFrame(seq_regions)
>>> seq_region_new = json_data_new['seq_regions']
>>> seq_regions_table_new = pd.DataFrame(seq_region_new)
>>> assert seq_regions_table_true.loc['phenotypes'].values.tolist() == seq_regions_table_new.loc['phenotypes'].values.tolist(), 'Phenotypes are not equal'
>>> assert seq_regions_table_true.loc['alignment_length'].values.tolist() == seq_regions_table_new.loc['alignment_length'].values.tolist(), f'Alignment length for ccoli is not equal'
>>> assert seq_regions_table_true.loc['ref_start_pos'].values.tolist() == seq_regions_table_new.loc['ref_start_pos'].values.tolist(), f'Ref start for ccoli is not equal'
>>> assert seq_regions_table_true.loc['ref_end_pos'].values.tolist() == seq_regions_table_new.loc['ref_end_pos'].values.tolist(), f'Ref end for ccoli is not equal'

```

# Test mutations for ccoli and mut_type: unknown

```
>>> seq_variations_true = json_data_true['seq_variations']
>>> seq_variations_table_true = pd.DataFrame(seq_variations_true)
>>> seq_variations_new = json_data_new['seq_variations']
>>> seq_variations_table_new = pd.DataFrame(seq_variations_new)
>>> assert seq_variations_table_true.shape == seq_variations_table_new.shape, 'Shape of seq_variations is not equal'
>>> if len(seq_variations_table_true) > 0:
...		assert seq_variations_table_new.loc['key'].values.tolist()[0].split(';')[0] == seq_variations_table_true.loc['key'].values.tolist()[0].split(';')[0], 'First element of key is not equal'
...		assert seq_variations_table_new.loc['ref_start_pos'].values.tolist() == seq_variations_table_true.loc['ref_start_pos'].values.tolist(), 'Mutation start position for unknown and ccoli is not equal'
...		assert seq_variations_table_new.loc['ref_end_pos'].values.tolist() == seq_variations_table_true.loc['ref_end_pos'].values.tolist(), 'Mutation end position for unknown and ccoli is not equal'
...		assert seq_variations_table_new.loc['var_codon'].values.tolist() == seq_variations_table_true.loc['var_codon'].values.tolist(), 'Mutation codon for unknown and ccoli is not equal'
...		assert seq_variations_table_new.loc['ref_codon'].values.tolist() == seq_variations_table_true.loc['ref_codon'].values.tolist(), 'Reference codon for unknown and ccoli is not equal'

```

```
>>> os.remove('/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_02.json')

```

## Software Command for senterica

```
>>> subprocess.run(['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'senterica', '-u', '-acq', '-ifa', '/home/people/s220672/resfinder/tests/data/test_isolate_02.fa', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_02.json'])
CompletedProcess(args=['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'senterica', '-u', '-acq', '-ifa', '/home/people/s220672/resfinder/tests/data/test_isolate_02.fa', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_02.json'], returncode=0)
>>> json_file = '/home/people/s220672/functional_testing/functional_testing_data/senterica/test_isolate_02/whimv/test_isolate_02_-acq-u.json'
>>> with open(json_file, 'r') as f:
...		json_data_true = json.load(f)
>>> json_file_new = '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_02.json'
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

### Test seq_regions for senterica

```
>>> seq_regions = json_data_true['seq_regions']
>>> seq_regions_table_true = pd.DataFrame(seq_regions)
>>> seq_region_new = json_data_new['seq_regions']
>>> seq_regions_table_new = pd.DataFrame(seq_region_new)
>>> assert seq_regions_table_true.loc['phenotypes'].values.tolist() == seq_regions_table_new.loc['phenotypes'].values.tolist(), 'Phenotypes are not equal'
>>> assert seq_regions_table_true.loc['alignment_length'].values.tolist() == seq_regions_table_new.loc['alignment_length'].values.tolist(), f'Alignment length for senterica is not equal'
>>> assert seq_regions_table_true.loc['ref_start_pos'].values.tolist() == seq_regions_table_new.loc['ref_start_pos'].values.tolist(), f'Ref start for senterica is not equal'
>>> assert seq_regions_table_true.loc['ref_end_pos'].values.tolist() == seq_regions_table_new.loc['ref_end_pos'].values.tolist(), f'Ref end for senterica is not equal'

```

# Test mutations for senterica and mut_type: unknown

```
>>> seq_variations_true = json_data_true['seq_variations']
>>> seq_variations_table_true = pd.DataFrame(seq_variations_true)
>>> seq_variations_new = json_data_new['seq_variations']
>>> seq_variations_table_new = pd.DataFrame(seq_variations_new)
>>> assert seq_variations_table_true.shape == seq_variations_table_new.shape, 'Shape of seq_variations is not equal'
>>> if len(seq_variations_table_true) > 0:
...		assert seq_variations_table_new.loc['key'].values.tolist()[0].split(';')[0] == seq_variations_table_true.loc['key'].values.tolist()[0].split(';')[0], 'First element of key is not equal'
...		assert seq_variations_table_new.loc['ref_start_pos'].values.tolist() == seq_variations_table_true.loc['ref_start_pos'].values.tolist(), 'Mutation start position for unknown and senterica is not equal'
...		assert seq_variations_table_new.loc['ref_end_pos'].values.tolist() == seq_variations_table_true.loc['ref_end_pos'].values.tolist(), 'Mutation end position for unknown and senterica is not equal'
...		assert seq_variations_table_new.loc['var_codon'].values.tolist() == seq_variations_table_true.loc['var_codon'].values.tolist(), 'Mutation codon for unknown and senterica is not equal'
...		assert seq_variations_table_new.loc['ref_codon'].values.tolist() == seq_variations_table_true.loc['ref_codon'].values.tolist(), 'Reference codon for unknown and senterica is not equal'

```

```
>>> os.remove('/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_02.json')

```

## Software Command for ecoli

```
>>> subprocess.run(['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'ecoli', '-c', '-ifa', '/home/people/s220672/resfinder/tests/data/test_isolate_05.fa', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_05.json'])
CompletedProcess(args=['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'ecoli', '-c', '-ifa', '/home/people/s220672/resfinder/tests/data/test_isolate_05.fa', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_05.json'], returncode=0)
>>> json_file = '/home/people/s220672/functional_testing/functional_testing_data/ecoli/test_isolate_05/GOhqy/test_isolate_05_-c.json'
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

### Test seq_regions for ecoli

```
>>> seq_regions = json_data_true['seq_regions']
>>> seq_regions_table_true = pd.DataFrame(seq_regions)
>>> seq_region_new = json_data_new['seq_regions']
>>> seq_regions_table_new = pd.DataFrame(seq_region_new)
>>> assert seq_regions_table_true.loc['phenotypes'].values.tolist() == seq_regions_table_new.loc['phenotypes'].values.tolist(), 'Phenotypes are not equal'
>>> assert seq_regions_table_true.loc['alignment_length'].values.tolist() == seq_regions_table_new.loc['alignment_length'].values.tolist(), f'Alignment length for ecoli is not equal'
>>> assert seq_regions_table_true.loc['ref_start_pos'].values.tolist() == seq_regions_table_new.loc['ref_start_pos'].values.tolist(), f'Ref start for ecoli is not equal'
>>> assert seq_regions_table_true.loc['ref_end_pos'].values.tolist() == seq_regions_table_new.loc['ref_end_pos'].values.tolist(), f'Ref end for ecoli is not equal'

```

# Test mutations for ecoli and mut_type: chromosomal

```
>>> seq_variations_true = json_data_true['seq_variations']
>>> seq_variations_table_true = pd.DataFrame(seq_variations_true)
>>> seq_variations_new = json_data_new['seq_variations']
>>> seq_variations_table_new = pd.DataFrame(seq_variations_new)
>>> assert seq_variations_table_true.shape == seq_variations_table_new.shape, 'Shape of seq_variations is not equal'
>>> if len(seq_variations_table_true) > 0:
...		assert seq_variations_table_new.loc['key'].values.tolist()[0].split(';')[0] == seq_variations_table_true.loc['key'].values.tolist()[0].split(';')[0], 'First element of key is not equal'
...		assert seq_variations_table_new.loc['ref_start_pos'].values.tolist() == seq_variations_table_true.loc['ref_start_pos'].values.tolist(), 'Mutation start position for chromosomal and ecoli is not equal'
...		assert seq_variations_table_new.loc['ref_end_pos'].values.tolist() == seq_variations_table_true.loc['ref_end_pos'].values.tolist(), 'Mutation end position for chromosomal and ecoli is not equal'
...		assert seq_variations_table_new.loc['var_codon'].values.tolist() == seq_variations_table_true.loc['var_codon'].values.tolist(), 'Mutation codon for chromosomal and ecoli is not equal'
...		assert seq_variations_table_new.loc['ref_codon'].values.tolist() == seq_variations_table_true.loc['ref_codon'].values.tolist(), 'Reference codon for chromosomal and ecoli is not equal'

```

### Test if resistant gene is same after splitting in two contigs.

```
>>> with open('/home/people/s220672/resfinder/tests/data/test_isolate_05.fa', 'r') as f:
...		fasta = f.read()
>>> ModFasta = FastaModifier()
>>> both_seqs, both_headers = ModFasta.modify_seqs(json_data_true, '/home/people/s220672/resfinder/tests/data/test_isolate_05.fa')
>>> ModFasta.create_fasta('test_temp/new_fasta', both_seqs, both_headers)
>>> subprocess.run(['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'ecoli', '-c', '-ifa', 'test_temp/new_fasta.fa', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_05.json'])
CompletedProcess(args=['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'ecoli', '-c', '-ifa', 'test_temp/new_fasta.fa', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_05.json'], returncode=0)
>>> json_file_new = '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_05.json'
>>> with open(json_file_new, 'r') as f:
...		json_data_new = json.load(f)
>>> seq_variations = json_data_true['seq_variations']
>>> seq_variations_table_true = pd.DataFrame(seq_variations)
>>> if len(seq_variations_table_true) > 0:
...		seq_variations_new = json_data_new['seq_variations']
...		seq_variations_table_new = pd.DataFrame(seq_variations_new)
...		assert seq_variations_table_true.loc['seq_var'].values.tolist() == seq_variations_table_new.loc['seq_var'].values.tolist(), f'Reference sequence length for ecoli is not equal'
...		assert seq_variations_table_true.loc['ref_end_pos'].values.tolist() == seq_variations_table_new.loc['ref_end_pos'].values.tolist(), f'Ref id for ecoli is not equal'

```

```
>>> os.remove('/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_05.json')

```

## Software Command for senterica

```
>>> subprocess.run(['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'senterica', '-c', '-ifa', '/home/people/s220672/resfinder/tests/data/test_isolate_05.fa', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_05.json'])
CompletedProcess(args=['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'senterica', '-c', '-ifa', '/home/people/s220672/resfinder/tests/data/test_isolate_05.fa', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_05.json'], returncode=0)
>>> json_file = '/home/people/s220672/functional_testing/functional_testing_data/senterica/test_isolate_05/KrMge/test_isolate_05_-c.json'
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

### Test seq_regions for senterica

```
>>> seq_regions = json_data_true['seq_regions']
>>> seq_regions_table_true = pd.DataFrame(seq_regions)
>>> seq_region_new = json_data_new['seq_regions']
>>> seq_regions_table_new = pd.DataFrame(seq_region_new)
>>> assert seq_regions_table_true.loc['phenotypes'].values.tolist() == seq_regions_table_new.loc['phenotypes'].values.tolist(), 'Phenotypes are not equal'
>>> assert seq_regions_table_true.loc['alignment_length'].values.tolist() == seq_regions_table_new.loc['alignment_length'].values.tolist(), f'Alignment length for senterica is not equal'
>>> assert seq_regions_table_true.loc['ref_start_pos'].values.tolist() == seq_regions_table_new.loc['ref_start_pos'].values.tolist(), f'Ref start for senterica is not equal'
>>> assert seq_regions_table_true.loc['ref_end_pos'].values.tolist() == seq_regions_table_new.loc['ref_end_pos'].values.tolist(), f'Ref end for senterica is not equal'

```

# Test mutations for senterica and mut_type: chromosomal

```
>>> seq_variations_true = json_data_true['seq_variations']
>>> seq_variations_table_true = pd.DataFrame(seq_variations_true)
>>> seq_variations_new = json_data_new['seq_variations']
>>> seq_variations_table_new = pd.DataFrame(seq_variations_new)
>>> assert seq_variations_table_true.shape == seq_variations_table_new.shape, 'Shape of seq_variations is not equal'
>>> if len(seq_variations_table_true) > 0:
...		assert seq_variations_table_new.loc['key'].values.tolist()[0].split(';')[0] == seq_variations_table_true.loc['key'].values.tolist()[0].split(';')[0], 'First element of key is not equal'
...		assert seq_variations_table_new.loc['ref_start_pos'].values.tolist() == seq_variations_table_true.loc['ref_start_pos'].values.tolist(), 'Mutation start position for chromosomal and senterica is not equal'
...		assert seq_variations_table_new.loc['ref_end_pos'].values.tolist() == seq_variations_table_true.loc['ref_end_pos'].values.tolist(), 'Mutation end position for chromosomal and senterica is not equal'
...		assert seq_variations_table_new.loc['var_codon'].values.tolist() == seq_variations_table_true.loc['var_codon'].values.tolist(), 'Mutation codon for chromosomal and senterica is not equal'
...		assert seq_variations_table_new.loc['ref_codon'].values.tolist() == seq_variations_table_true.loc['ref_codon'].values.tolist(), 'Reference codon for chromosomal and senterica is not equal'

```

```
>>> os.remove('/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_05.json')

```

## Software Command for ecoli

```
>>> subprocess.run(['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'ecoli', '-c', '-ifa', '/home/people/s220672/resfinder/tests/data/test_isolate_09b.fa', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_09b.json'])
CompletedProcess(args=['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'ecoli', '-c', '-ifa', '/home/people/s220672/resfinder/tests/data/test_isolate_09b.fa', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_09b.json'], returncode=0)
>>> json_file = '/home/people/s220672/functional_testing/functional_testing_data/ecoli/test_isolate_09b/gioop/test_isolate_09b_-c.json'
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

### Test seq_regions for ecoli

```
>>> seq_regions = json_data_true['seq_regions']
>>> seq_regions_table_true = pd.DataFrame(seq_regions)
>>> seq_region_new = json_data_new['seq_regions']
>>> seq_regions_table_new = pd.DataFrame(seq_region_new)
>>> assert seq_regions_table_true.loc['phenotypes'].values.tolist() == seq_regions_table_new.loc['phenotypes'].values.tolist(), 'Phenotypes are not equal'
>>> assert seq_regions_table_true.loc['alignment_length'].values.tolist() == seq_regions_table_new.loc['alignment_length'].values.tolist(), f'Alignment length for ecoli is not equal'
>>> assert seq_regions_table_true.loc['ref_start_pos'].values.tolist() == seq_regions_table_new.loc['ref_start_pos'].values.tolist(), f'Ref start for ecoli is not equal'
>>> assert seq_regions_table_true.loc['ref_end_pos'].values.tolist() == seq_regions_table_new.loc['ref_end_pos'].values.tolist(), f'Ref end for ecoli is not equal'

```

# Test mutations for ecoli and mut_type: chromosomal

```
>>> seq_variations_true = json_data_true['seq_variations']
>>> seq_variations_table_true = pd.DataFrame(seq_variations_true)
>>> seq_variations_new = json_data_new['seq_variations']
>>> seq_variations_table_new = pd.DataFrame(seq_variations_new)
>>> assert seq_variations_table_true.shape == seq_variations_table_new.shape, 'Shape of seq_variations is not equal'
>>> if len(seq_variations_table_true) > 0:
...		assert seq_variations_table_new.loc['key'].values.tolist()[0].split(';')[0] == seq_variations_table_true.loc['key'].values.tolist()[0].split(';')[0], 'First element of key is not equal'
...		assert seq_variations_table_new.loc['ref_start_pos'].values.tolist() == seq_variations_table_true.loc['ref_start_pos'].values.tolist(), 'Mutation start position for chromosomal and ecoli is not equal'
...		assert seq_variations_table_new.loc['ref_end_pos'].values.tolist() == seq_variations_table_true.loc['ref_end_pos'].values.tolist(), 'Mutation end position for chromosomal and ecoli is not equal'
...		assert seq_variations_table_new.loc['var_codon'].values.tolist() == seq_variations_table_true.loc['var_codon'].values.tolist(), 'Mutation codon for chromosomal and ecoli is not equal'
...		assert seq_variations_table_new.loc['ref_codon'].values.tolist() == seq_variations_table_true.loc['ref_codon'].values.tolist(), 'Reference codon for chromosomal and ecoli is not equal'

```

### Test if resistant gene is same after splitting in two contigs.

```
>>> with open('/home/people/s220672/resfinder/tests/data/test_isolate_09b.fa', 'r') as f:
...		fasta = f.read()
>>> ModFasta = FastaModifier()
>>> both_seqs, both_headers = ModFasta.modify_seqs(json_data_true, '/home/people/s220672/resfinder/tests/data/test_isolate_09b.fa')
>>> ModFasta.create_fasta('test_temp/new_fasta', both_seqs, both_headers)
>>> subprocess.run(['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'ecoli', '-c', '-ifa', 'test_temp/new_fasta.fa', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_09b.json'])
CompletedProcess(args=['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'ecoli', '-c', '-ifa', 'test_temp/new_fasta.fa', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_09b.json'], returncode=0)
>>> json_file_new = '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_09b.json'
>>> with open(json_file_new, 'r') as f:
...		json_data_new = json.load(f)
>>> seq_variations = json_data_true['seq_variations']
>>> seq_variations_table_true = pd.DataFrame(seq_variations)
>>> if len(seq_variations_table_true) > 0:
...		seq_variations_new = json_data_new['seq_variations']
...		seq_variations_table_new = pd.DataFrame(seq_variations_new)
...		assert seq_variations_table_true.loc['seq_var'].values.tolist() == seq_variations_table_new.loc['seq_var'].values.tolist(), f'Reference sequence length for ecoli is not equal'
...		assert seq_variations_table_true.loc['ref_end_pos'].values.tolist() == seq_variations_table_new.loc['ref_end_pos'].values.tolist(), f'Ref id for ecoli is not equal'

```

```
>>> os.remove('/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_09b.json')

```

## Software Command for senterica

```
>>> subprocess.run(['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'senterica', '-c', '-ifa', '/home/people/s220672/resfinder/tests/data/test_isolate_09b.fa', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_09b.json'])
CompletedProcess(args=['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'senterica', '-c', '-ifa', '/home/people/s220672/resfinder/tests/data/test_isolate_09b.fa', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_09b.json'], returncode=0)
>>> json_file = '/home/people/s220672/functional_testing/functional_testing_data/senterica/test_isolate_09b/xXxyt/test_isolate_09b_-c.json'
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

### Test seq_regions for senterica

```
>>> seq_regions = json_data_true['seq_regions']
>>> seq_regions_table_true = pd.DataFrame(seq_regions)
>>> seq_region_new = json_data_new['seq_regions']
>>> seq_regions_table_new = pd.DataFrame(seq_region_new)
>>> assert seq_regions_table_true.loc['phenotypes'].values.tolist() == seq_regions_table_new.loc['phenotypes'].values.tolist(), 'Phenotypes are not equal'
>>> assert seq_regions_table_true.loc['alignment_length'].values.tolist() == seq_regions_table_new.loc['alignment_length'].values.tolist(), f'Alignment length for senterica is not equal'
>>> assert seq_regions_table_true.loc['ref_start_pos'].values.tolist() == seq_regions_table_new.loc['ref_start_pos'].values.tolist(), f'Ref start for senterica is not equal'
>>> assert seq_regions_table_true.loc['ref_end_pos'].values.tolist() == seq_regions_table_new.loc['ref_end_pos'].values.tolist(), f'Ref end for senterica is not equal'

```

# Test mutations for senterica and mut_type: chromosomal

```
>>> seq_variations_true = json_data_true['seq_variations']
>>> seq_variations_table_true = pd.DataFrame(seq_variations_true)
>>> seq_variations_new = json_data_new['seq_variations']
>>> seq_variations_table_new = pd.DataFrame(seq_variations_new)
>>> assert seq_variations_table_true.shape == seq_variations_table_new.shape, 'Shape of seq_variations is not equal'
>>> if len(seq_variations_table_true) > 0:
...		assert seq_variations_table_new.loc['key'].values.tolist()[0].split(';')[0] == seq_variations_table_true.loc['key'].values.tolist()[0].split(';')[0], 'First element of key is not equal'
...		assert seq_variations_table_new.loc['ref_start_pos'].values.tolist() == seq_variations_table_true.loc['ref_start_pos'].values.tolist(), 'Mutation start position for chromosomal and senterica is not equal'
...		assert seq_variations_table_new.loc['ref_end_pos'].values.tolist() == seq_variations_table_true.loc['ref_end_pos'].values.tolist(), 'Mutation end position for chromosomal and senterica is not equal'
...		assert seq_variations_table_new.loc['var_codon'].values.tolist() == seq_variations_table_true.loc['var_codon'].values.tolist(), 'Mutation codon for chromosomal and senterica is not equal'
...		assert seq_variations_table_new.loc['ref_codon'].values.tolist() == seq_variations_table_true.loc['ref_codon'].values.tolist(), 'Reference codon for chromosomal and senterica is not equal'

```

### Test if resistant gene is same after splitting in two contigs.

```
>>> with open('/home/people/s220672/resfinder/tests/data/test_isolate_09b.fa', 'r') as f:
...		fasta = f.read()
>>> ModFasta = FastaModifier()
>>> both_seqs, both_headers = ModFasta.modify_seqs(json_data_true, '/home/people/s220672/resfinder/tests/data/test_isolate_09b.fa')
>>> ModFasta.create_fasta('test_temp/new_fasta', both_seqs, both_headers)
>>> subprocess.run(['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'senterica', '-c', '-ifa', 'test_temp/new_fasta.fa', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_09b.json'])
CompletedProcess(args=['python', '/home/people/s220672/functional_testing/.venv/lib/python3.10/site-packages/resfinder/run_resfinder.py', '-s', 'senterica', '-c', '-ifa', 'test_temp/new_fasta.fa', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_09b.json'], returncode=0)
>>> json_file_new = '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_09b.json'
>>> with open(json_file_new, 'r') as f:
...		json_data_new = json.load(f)
>>> seq_variations = json_data_true['seq_variations']
>>> seq_variations_table_true = pd.DataFrame(seq_variations)
>>> if len(seq_variations_table_true) > 0:
...		seq_variations_new = json_data_new['seq_variations']
...		seq_variations_table_new = pd.DataFrame(seq_variations_new)
...		assert seq_variations_table_true.loc['seq_var'].values.tolist() == seq_variations_table_new.loc['seq_var'].values.tolist(), f'Reference sequence length for senterica is not equal'
...		assert seq_variations_table_true.loc['ref_end_pos'].values.tolist() == seq_variations_table_new.loc['ref_end_pos'].values.tolist(), f'Ref id for senterica is not equal'

```

```
>>> os.remove('/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_09b.json')

```

