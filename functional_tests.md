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
...			with open(output_file, 'w') as fasta_file:
...				SeqIO.write(records, fasta_file, 'fasta')
...		@staticmethod
...		def calc_positions(json_data_true):
...			seq_variations_true = json_data_true['seq_variations']
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
>>> subprocess.run(['/home/people/s220672/virulencefinder/.venv/bin/python', '-m', 'src.virulencefinder.__main__', '-p', '/home/people/s220672/databases/virulencefinder_db', '-ifq', '/home/people/s220672/virulencefinder/tests/data/e_coli/test_isolate_e_coli1_ne.fq.gz', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_e_coli1_ne.json', '-l', '0.6'])
CompletedProcess(args=['/home/people/s220672/virulencefinder/.venv/bin/python', '-m', 'src.virulencefinder.__main__', '-p', '/home/people/s220672/databases/virulencefinder_db', '-ifq', '/home/people/s220672/virulencefinder/tests/data/e_coli/test_isolate_e_coli1_ne.fq.gz', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_e_coli1_ne.json', '-l', '0.6'], returncode=0)
>>> json_file = '/home/people/s220672/functional_testing/functional_testing_data/ecoli/test_isolate_e_coli1_ne.json'
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

## Software Command for all

```
>>> subprocess.run(['/home/people/s220672/virulencefinder/.venv/bin/python', '-m', 'src.virulencefinder.__main__', '-p', '/home/people/s220672/databases/virulencefinder_db', '-ifq', '/home/people/s220672/virulencefinder/tests/data/e_coli/test_isolate_e_coli1_ne.fq.gz', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_e_coli1_ne.json', '-l', '0.7'])
CompletedProcess(args=['/home/people/s220672/virulencefinder/.venv/bin/python', '-m', 'src.virulencefinder.__main__', '-p', '/home/people/s220672/databases/virulencefinder_db', '-ifq', '/home/people/s220672/virulencefinder/tests/data/e_coli/test_isolate_e_coli1_ne.fq.gz', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_e_coli1_ne.json', '-l', '0.7'], returncode=0)
>>> json_file = '/home/people/s220672/functional_testing/functional_testing_data/ecoli/test_isolate_e_coli1_ne.json'
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

## Software Command for all

```
>>> subprocess.run(['/home/people/s220672/virulencefinder/.venv/bin/python', '-m', 'src.virulencefinder.__main__', '-p', '/home/people/s220672/databases/virulencefinder_db', '-ifq', '/home/people/s220672/virulencefinder/tests/data/e_coli/test_isolate_e_coli1_ne.fq.gz', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_e_coli1_ne.json', '-l', '0.8'])
CompletedProcess(args=['/home/people/s220672/virulencefinder/.venv/bin/python', '-m', 'src.virulencefinder.__main__', '-p', '/home/people/s220672/databases/virulencefinder_db', '-ifq', '/home/people/s220672/virulencefinder/tests/data/e_coli/test_isolate_e_coli1_ne.fq.gz', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_e_coli1_ne.json', '-l', '0.8'], returncode=0)
>>> json_file = '/home/people/s220672/functional_testing/functional_testing_data/ecoli/test_isolate_e_coli1_ne.json'
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

## Software Command for all

```
>>> subprocess.run(['/home/people/s220672/virulencefinder/.venv/bin/python', '-m', 'src.virulencefinder.__main__', '-p', '/home/people/s220672/databases/virulencefinder_db', '-ifq', '/home/people/s220672/virulencefinder/tests/data/e_coli/test_isolate_e_coli2_ne.fq.gz', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_e_coli2_ne.json', '-l', '0.6'])
CompletedProcess(args=['/home/people/s220672/virulencefinder/.venv/bin/python', '-m', 'src.virulencefinder.__main__', '-p', '/home/people/s220672/databases/virulencefinder_db', '-ifq', '/home/people/s220672/virulencefinder/tests/data/e_coli/test_isolate_e_coli2_ne.fq.gz', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_e_coli2_ne.json', '-l', '0.6'], returncode=0)
>>> json_file = '/home/people/s220672/functional_testing/functional_testing_data/ecoli/test_isolate_e_coli2_ne.json'
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

## Software Command for all

```
>>> subprocess.run(['/home/people/s220672/virulencefinder/.venv/bin/python', '-m', 'src.virulencefinder.__main__', '-p', '/home/people/s220672/databases/virulencefinder_db', '-ifq', '/home/people/s220672/virulencefinder/tests/data/e_coli/test_isolate_e_coli2_ne.fq.gz', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_e_coli2_ne.json', '-l', '0.7'])
CompletedProcess(args=['/home/people/s220672/virulencefinder/.venv/bin/python', '-m', 'src.virulencefinder.__main__', '-p', '/home/people/s220672/databases/virulencefinder_db', '-ifq', '/home/people/s220672/virulencefinder/tests/data/e_coli/test_isolate_e_coli2_ne.fq.gz', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_e_coli2_ne.json', '-l', '0.7'], returncode=0)
>>> json_file = '/home/people/s220672/functional_testing/functional_testing_data/ecoli/test_isolate_e_coli2_ne.json'
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

## Software Command for all

```
>>> subprocess.run(['/home/people/s220672/virulencefinder/.venv/bin/python', '-m', 'src.virulencefinder.__main__', '-p', '/home/people/s220672/databases/virulencefinder_db', '-ifq', '/home/people/s220672/virulencefinder/tests/data/e_coli/test_isolate_e_coli2_ne.fq.gz', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_e_coli2_ne.json', '-l', '0.8'])
CompletedProcess(args=['/home/people/s220672/virulencefinder/.venv/bin/python', '-m', 'src.virulencefinder.__main__', '-p', '/home/people/s220672/databases/virulencefinder_db', '-ifq', '/home/people/s220672/virulencefinder/tests/data/e_coli/test_isolate_e_coli2_ne.fq.gz', '-o', '/home/people/s220672/functional_testing/test_temp', '-j', '/home/people/s220672/functional_testing/functional_testing_temp/test_isolate_e_coli2_ne.json', '-l', '0.8'], returncode=0)
>>> json_file = '/home/people/s220672/functional_testing/functional_testing_data/ecoli/test_isolate_e_coli2_ne.json'
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

