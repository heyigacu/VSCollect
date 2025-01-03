# 1. AI screen

## 1.1 RosettaVS
![flowchart](RosettaVS/flowchart.png)
### 1.1.1 prepare
* files need:  
zinc17.csv # cluster library
zinc22.csv # screen library
please put receptor.pdb, ligand.mol2 (refrence ligand) in task dir

* packages need:
```
git clone https://github.com/RosettaCommons/rosetta
cd rosetta
```

### 1.1.2 run 

* generate druglike-centroid library from zinc library

```
python gen_centroid.py
```
will generate druglike_centroid_library.csv  
if you don't want this step, please change the name of your screened library to druglike_centroid_library.csv

```
python main.py
```