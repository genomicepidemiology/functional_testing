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
>>> subprocess.run(['/home/people/s220672/virulencefinder/.venv/bin/python', '-m', 'src.virulencefinder.__main__', '-p', '/home/people/s220672/databases/virulencefinder_db', '-ifq', '/home/people/s220672/virulencefinder/tests/data/e_coli/test_isolate_e_coli1_ne.fq.gz', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_e_coli1_ne.json'])
CompletedProcess(args=['/home/people/s220672/virulencefinder/.venv/bin/python', '-m', 'src.virulencefinder.__main__', '-p', '/home/people/s220672/databases/virulencefinder_db', '-ifq', '/home/people/s220672/virulencefinder/tests/data/e_coli/test_isolate_e_coli1_ne.fq.gz', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_e_coli1_ne.json'], returncode=0)
>>> json_file = '/home/people/s220672/functional_testing/functional_testing_data/ecoli/test_isolate_e_coli1_ne/bEZUP/test_isolate_e_coli1_ne.json'
>>> with open(json_file, 'r') as f:
...		json_data_true = json.load(f)
>>> json_file_new = '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_e_coli1_ne.json'
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

### Custom test

```
>>> seq_regions = json_data_true['seq_regions']
>>> seq_regions_table_true = pd.DataFrame(seq_regions)
>>> if seq_regions_table_true.shape[0] > 0:
...		seq_region_new = json_data_new['seq_regions']
...		seq_regions_table_new = pd.DataFrame(seq_region_new)
...		assert seq_regions_table_true.loc['phenotypes'].values.tolist() == seq_regions_table_new.loc['phenotypes'].values.tolist(), 'Phenotypes are not equal'
...		assert seq_regions_table_true.loc['alignment_length'].values.tolist() == seq_regions_table_new.loc['alignment_length'].values.tolist(), f'Alignment length for all is not equal'
...		assert seq_regions_table_true.loc['ref_start_pos'].values.tolist() == seq_regions_table_new.loc['ref_start_pos'].values.tolist(), f'Ref start for all is not equal'
...		assert seq_regions_table_true.loc['ref_end_pos'].values.tolist() == seq_regions_table_new.loc['ref_end_pos'].values.tolist(), f'Ref end for all is not equal'
>>> else:
...		pass
>>> seq_variations_true = json_data_true['seq_variations']
>>> seq_variations_table_true = pd.DataFrame(seq_variations_true)
>>> if seq_variations_table_true.shape[0] > 0:
...		seq_variations_new = json_data_new['seq_variations']
...		seq_variations_table_new = pd.DataFrame(seq_variations_new)
...		assert seq_variations_table_true.loc['ref_end_pos'].values.tolist() == seq_variations_table_new.loc['ref_end_pos'].values.tolist(), 'Ref end position is not equal'
>>> phenotypes_true = json_data_true['phenotypes']
>>> phenotypes_table_true = pd.DataFrame(phenotypes_true)
>>> if phenotypes_table_true.shape[0] > 0:
...		seq_variations_new = json_data_new['phenotypes']
...		seq_variations_table_new = pd.DataFrame(seq_variations_new)
...		assert seq_variations_table_true.loc['phenotypes'].values.tolist() == seq_variations_table_new.loc['phenotypes'].values.tolist(), 'Phenotypes are not equal'

```

```
>>> os.remove('/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_e_coli1_ne.json')

```

## Software Command for all

```
>>> subprocess.run(['/home/people/s220672/virulencefinder/.venv/bin/python', '-m', 'src.virulencefinder.__main__', '-p', '/home/people/s220672/databases/virulencefinder_db', '-ifq', '/home/people/s220672/virulencefinder/tests/data/e_coli/test_isolate_e_coli2_ne.fq.gz', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_e_coli2_ne.json'])
CompletedProcess(args=['/home/people/s220672/virulencefinder/.venv/bin/python', '-m', 'src.virulencefinder.__main__', '-p', '/home/people/s220672/databases/virulencefinder_db', '-ifq', '/home/people/s220672/virulencefinder/tests/data/e_coli/test_isolate_e_coli2_ne.fq.gz', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_e_coli2_ne.json'], returncode=0)
>>> json_file = '/home/people/s220672/functional_testing/functional_testing_data/ecoli/test_isolate_e_coli2_ne/fURcU/test_isolate_e_coli2_ne.json'
>>> with open(json_file, 'r') as f:
...		json_data_true = json.load(f)
>>> json_file_new = '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_e_coli2_ne.json'
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

### Custom test

```
>>> seq_regions = json_data_true['seq_regions']
>>> seq_regions_table_true = pd.DataFrame(seq_regions)
>>> if seq_regions_table_true.shape[0] > 0:
...		seq_region_new = json_data_new['seq_regions']
...		seq_regions_table_new = pd.DataFrame(seq_region_new)
...		assert seq_regions_table_true.loc['phenotypes'].values.tolist() == seq_regions_table_new.loc['phenotypes'].values.tolist(), 'Phenotypes are not equal'
...		assert seq_regions_table_true.loc['alignment_length'].values.tolist() == seq_regions_table_new.loc['alignment_length'].values.tolist(), f'Alignment length for all is not equal'
...		assert seq_regions_table_true.loc['ref_start_pos'].values.tolist() == seq_regions_table_new.loc['ref_start_pos'].values.tolist(), f'Ref start for all is not equal'
...		assert seq_regions_table_true.loc['ref_end_pos'].values.tolist() == seq_regions_table_new.loc['ref_end_pos'].values.tolist(), f'Ref end for all is not equal'
>>> else:
...		pass
>>> seq_variations_true = json_data_true['seq_variations']
>>> seq_variations_table_true = pd.DataFrame(seq_variations_true)
>>> if seq_variations_table_true.shape[0] > 0:
...		seq_variations_new = json_data_new['seq_variations']
...		seq_variations_table_new = pd.DataFrame(seq_variations_new)
...		assert seq_variations_table_true.loc['ref_end_pos'].values.tolist() == seq_variations_table_new.loc['ref_end_pos'].values.tolist(), 'Ref end position is not equal'
>>> phenotypes_true = json_data_true['phenotypes']
>>> phenotypes_table_true = pd.DataFrame(phenotypes_true)
>>> if phenotypes_table_true.shape[0] > 0:
...		seq_variations_new = json_data_new['phenotypes']
...		seq_variations_table_new = pd.DataFrame(seq_variations_new)
...		assert seq_variations_table_true.loc['phenotypes'].values.tolist() == seq_variations_table_new.loc['phenotypes'].values.tolist(), 'Phenotypes are not equal'

```

```
>>> os.remove('/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_e_coli2_ne.json')

```

