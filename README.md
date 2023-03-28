# 3D-MCTS
3D-MCTS: A general structure-based molecule generation method with MCTS.

## How to sample molecules for a specific protein

### Step1: Prepare several files needed

file 1. The protein structure file, PDB format. (Without ligand atoms)

file 2. The ligand file, sdf format. (To determine the position of binding site.)

file 3. The pocket file, PDB format. (To speed the calculation, it can be replaced by file 1.)

### Step2: Prepare the fragment library

You can use the library that we provide in frags/fragment.txt or modify it according to your need.

### Step2: Prepare the initial fragment

We provide three initial fragments in the directory '''init/'''. You can provide initial fragments according to your need.
