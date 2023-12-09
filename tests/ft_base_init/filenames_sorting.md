## Test sort_filenames

```python
>>> from src.functional_testing.filenames_sorting import CreateOutputs
>>> test_match = ['/path/to/some/dir/aaaaa_1.fastq', '/path/to/some/dir/aaaaa_2.fastq', '/path/to/some/dir/aaaae_1.fastq', '/path/to/some/dir/aaaae_R1','/path/to/some/dir/aaaae_R2']
>>> output_dic = CreateOutputs.sort_filenames(test_match)
>>> output_dic
{'aaaaa_1': ['/path/to/some/dir/aaaaa_1.fastq', '/path/to/some/dir/aaaaa_2.fastq'], 'aaaae_1': ['/path/to/some/dir/aaaae_1.fastq'], 'aaaae_R1': ['/path/to/some/dir/aaaae_R1', '/path/to/some/dir/aaaae_R2']}

```