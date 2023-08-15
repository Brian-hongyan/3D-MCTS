# 3D-MCTS
3D-MCTS: A general structure-based molecule generation method with MCTS.

![error](https://github.com/Brian-hongyan/3D-MCTS/blob/main/method.jpg)

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
