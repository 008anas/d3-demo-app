# bioroboost_functions
Checking functions for bioroboost project


Example:
```bash

# Scoring:
python3 bioroboost.py -f1 -s CGGCGAGCGGCGAGTAAAATATATAAAATGAGCGGAGAGCG -m ./data/pribnow_mpneumoniae.txt

# Run evaluation of a sequence with default './promoters_TATAAT.fa' training set
python3 TATA_test.py TATAAT
python3 TATA_test.py AAAAAA

# Run evaluation of a sequence with custom training set
python3 TATA_test.py TATAAT ./promoters_TATAAT.fa

```
## Tools needed
```
Bio
Muscle
```
