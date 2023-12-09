## Test read_fasta

```
>>> from src.functional_testing.helper_scripts.split_fasta import FastaModifier
>>> import os
>>> Modifier = FastaModifier()
>>> cwd = os.getcwd()
>>> path_fasta = os.path.join(cwd, 'tests', 'data', 'test_isolate_02.fa')
>>> sequence, header = Modifier.read_fasta(path_fasta)
>>> assert sequence.startswith(">") == False
>>> assert header.startswith(">") == True
>>> assert type(sequence) == str

```

## Test calc_positions

```
>>> import json
>>> with open(os.path.join(cwd,'tests', 'data', 'test_isolate_02_-acq-u.json')) as f:
...     data_json = json.load(f)
>>> start_pos, end_pos, half_pos = Modifier.calc_positions(data_json)
>>> assert type(start_pos) == int
>>> assert type(end_pos) == int
>>> assert type(half_pos) == int
>>> assert start_pos == 1
>>> assert end_pos == 876
>>> assert half_pos == 437

```

## Test modify_seqs

```
>>> both_seqs, both_headers = Modifier.modify_seqs(data_json, path_fasta)
>>> assert len(both_seqs) == 2
>>> assert len(both_headers) == 2
>>> assert type(both_seqs[0]) == str
>>> assert type(both_seqs) == list

```

## Test create_fasta

```
>>> Modifier.create_fasta("new_fasta.fa", both_seqs, both_headers)
>>> sequence,header = Modifier.read_fasta("new_fasta.fa")
>>> assert len(sequence) == (876 - 437)
>>> os.remove("new_fasta.fa")

```