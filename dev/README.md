# bioroboost
Checking functions for bioroboost project


## Scoring mode

In this mode we evaluate a sequence given a Position-Weight Matrix (PWM). Arguments required are a <sequence> in plain or fasta format and a <matrix> provided as a file with the format MxN where: 
 - M is the number of rows, i.e. nucleotide bases ('DNA', 'RNA') or amino acids ('PROT')
 - N is the number of columns, i.e. positions of the motif to be evaluated:

| residue | i1 | i2 | ... | in-1 | in |
|---------|----|----|-----|------|----|
| A       |    |    |     |      |    |
| C       |    |    |     |      |    |
| G       |    |    |     |      |    |
| T       |    |    |     |      |    |

Example:
```bash
# Scoring:
python3 bioroboost.py -f 1 -s CGGCGAGCGGCGAGTAAAATATATAAAATGAGCGGAGAGCG -m ./data/pribnow_mpneumoniae.txt     # Plain sequence
python3 bioroboost.py -f 1 -s ./data/test_sequence.fa -m ./data/pribnow_mpneumoniae.txt     # Fasta sequence
```

-f = function, use '1' or 'scoring' for this mode


Returns a json format compatible with ProtVista:
```
[{"start": "1", "score": "-46.73098889037074"}, {"start": "2", "score": "-33.03855102401263"}, {"start": "3", "score": "-49.085245529539776"}, {"start": "4", "score": "-31.785128118649908"}, {"start": "5", "score": "-52.44468382286683"}, {"start": "6", "score": "-48.426197633929114"}, {"start": "7", "score": "-51.96061171682282"}, {"start": "8", "score": "-45.11226550507895"}, {"start": "9", "score": "-27.411570857613878"}, {"start": "10", "score": "-52.85855567871745"}, {"start": "11", "score": "-28.07956190046838"}, {"start": "12", "score": "-59.50567644344946"}, {"start": "13", "score": "-43.55430894030511"}, {"start": "14", "score": "-0.9000824024886531"}, {"start": "15", "score": "-44.27981827539081"}, {"start": "16", "score": "-24.860527498733358"}, {"start": "17", "score": "-45.66709740485235"}, {"start": "18", "score": "-44.792096068057575"}, {"start": "19", "score": "-27.08226599105691"}, {"start": "20", "score": "-62.88706507369018"}, {"start": "21", "score": "-20.949590592598348"}, {"start": "22", "score": "-58.07871144625646"}, {"start": "23", "score": "-2.301328597185963"}, {"start": "24", "score": "-44.886764198137826"}, {"start": "25", "score": "-49.47658595818071"}, {"start": "26", "score": "-43.47767155368455"}, {"start": "27", "score": "-50.5174927144411"}, {"start": "28", "score": "-45.98008807248126"}, {"start": "29", "score": "-29.73760362458989"}, {"start": "30", "score": "-67.54897543052283"}, {"start": "31", "score": "-46.31585602360879"}, {"start": "32", "score": "-48.365239441170765"}, {"start": "33", "score": "-46.73098889037074"}, {"start": "34", "score": "-13.691944955409621"}, {"start": "35", "score": "-64.78239667687117"}, {"start": "36", "score": "-17.854577286077223"}, {"start": "37", "score": "-65.18956785472935"}, {"start": "38", "score": "-51.313341283756856"}, {"start": "39", "score": "-36.34036149767725"}, {"start": "40", "score": "-48.244590537425175"}, {"start": "41", "score": "-51.96061171682282"}]
```
