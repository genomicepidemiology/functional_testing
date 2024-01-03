## Created tests

Here you can download the created tests for different settings of the program. The tests were created with resfinder version 4.3.3

### Test 1 - A test script for all possible parameters


Output: test_all.md

**COMMAND:**
```
pdm run python src/functional_testing/__main__.py -d /home/people/s220672/resfinder/tests/data 

```

### Test 2 - A test script for the mutation parameters for ecoli

Output: test_mutation_ecoli.md

**COMMAND:**
```
pdm run python src/functional_testing/__main__.py -d /home/people/s220672/resfinder/tests/data -t mutation --values_list ecoli

```

### Test 3 - A test script for the species parameter for ecoli

Output: test_species_ecoli.md

**COMMAND:**
```
pdm run python src/functional_testing/__main__.py -d /home/people/s220672/resfinder/tests/data -t species --values_list ecoli

```

### Test 4 - A test script for processing the contigs correctly

Output: test_contigs.md

**COMMAND:**
```
pdm run python src/functional_testing/__main__.py -d /home/people/s220672/resfinder/tests/data -t contigs

```

