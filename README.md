# 3D-MCTS
3D-MCTS: A general structure-based molecule generation method with MCTS.

## Abstract

We present a novel search-based framework, 3D-MCTS, for structure-based de novo drug design. Distinct from prevailing atom-centric methods, 3D-MCTS employs a **fragment-based** molecular editing strategy. The fragments decomposed from small-molecule drugs are **recombined under predefined retrosynthetic rules**, offering improved drug-likeness and synthesizability, overcoming the inherent limitations of atom-based approaches. The integration of **multi-threaded parallel simulations and real-time energy constraint-based pruning strategy** equips 3D-MCTS with the capability to efficiently generate molecules with superior binding affinity (-2.0 kcal/mol better than state-of-the-art (SOTA) methods, yet at a comparable computational cost) and more reliable binding conformations (with a 43.6 % higher success rate than SOTAs). Such advantage becomes even more pronounced when handling targets that significantly deviate from the training dataset. Moreover, the flexibility of our framework makes it easy to incorporate domain knowledge during the process, thereby enabling the generation of molecules with desirable pharmacophore and enhanced binding affinity. The adaptability of 3D-MCTS is further showcased in metalloprotein applications, highlighting its potential across various drug design scenarios.

![error](https://github.com/Brian-hongyan/3D-MCTS/blob/main/method.jpg)

# Environment
```
# software
gnina         ## https://github.com/gnina/gnina
ADFR          ## https://ccsb.scripps.edu/adfr/

# python
python >= 3.7                ## conda
openbabel >= 3.1.1           ## conda
rdkit >= 2022.03.2           ## conda
func-timeout >= 4.3.5        ## pip


After installing gnina and ADFR, please edit the PATH for these two softwares in 3D-MCTS.py:

''
GNINA = '/home/hongyan/software/gnina'
ADFR = '/home/hongyan/software/ADFR/bin'
''

```


## How to sample molecules for a specific protein

### Step1: Prepare several files needed

file 1. The protein structure file, PDB format. (Without ligand atoms)

file 2. The ligand file, sdf format. (To determine the position of binding site.)

file 3. The pocket file, PDB format. (To speed the calculation. It can be replaced by file 1.)

### Step2: Prepare the fragment library

You can use the library that we provide in frags/fragment.txt or modify it according to your need.

### Step3: Prepare the initial fragment

We provide three initial fragments in the directory ```init/```. You can provide initial fragments according to your need.

### Step4: Run the Code.

```
python 3D-MCTS.py --num_sims 100000 --ligand ./ligand.sdf --protein ./protein.pdb --pocket ./pocket.pdb --score -7 --qed 0.3 --processor 48 --start 1
```

Molecules that meet the criteria are saved in ```record/```.
