# Permutation Test
This is a script to identify gens significantly affected by multiple mutations (MMs) using permutation test.


## Dependency
```
python 2.7.15
python packages
annot_utils
pysam (0.14.1)
pandas (0.23.4)
numpy(1.15.3)
```

## Preparation
```
python lib/permutation_test/preparation2.py
python lib/permutation_test/preparation2.py
python lib/permutation_test/silent.py maf_file
```

## Run 
```
python2 lib/permutation_test/run_simulation.py SNV     maf_file output_file --gene interested_gene
python2 lib/permutation_test/run_simulation.py non_SNV maf_file output_file --gene interested_gene
```

