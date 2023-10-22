# :loudspeaker:3D-MCTS
[3D-MCTS](https://pubs.rsc.org/en/content/articlelanding/2023/SC/D3SC04091G): A Flexible Data-Free Framework for Structure-Based De Novo Drug Design with Reinforcement Learning.

## Abstract

We present a novel search-based framework, 3D-MCTS, for structure-based de novo drug design. Distinct from prevailing atom-centric methods, 3D-MCTS employs a **fragment-based** molecular editing strategy. The fragments decomposed from small-molecule drugs are **recombined under predefined retrosynthetic rules**, offering improved drug-likeness and synthesizability, overcoming the inherent limitations of atom-based approaches. The integration of **multi-threaded parallel simulations and real-time energy constraint-based pruning strategy** equips 3D-MCTS with the capability to efficiently generate molecules with superior binding affinity (-2.0 kcal/mol better than state-of-the-art (SOTA) methods, yet at a comparable computational cost) and more reliable binding conformations (with a 43.6 % higher success rate than SOTAs). 3D-MCTS is capable of achieving thirty times more hits with high binding affinity than traditional virtual screening methods, which demonstrates the superior ability of 3D-MCTS to explore chemical space.

![error](https://github.com/Brian-hongyan/3D-MCTS/blob/main/method.jpg)

# DataSet
The main data for benchmark is CrossDock2020, which is utilized by most of the methods. You can download the processed data from [this link](https://drive.google.com/drive/folders/1CzwxmTpjbrt83z_wBzcQncq84OVDPurM). This is the processed version of original files, which is processed by [Luoshi Tong](https://github.com/luost26/3D-Generative-SBDD/tree/main/data). 
# Environment
```
# software
gnina         ## https://github.com/gnina/gnina
ADFR          ## https://ccsb.scripps.edu/adfr/

# python
python >= 3.7                ## 
openbabel >= 3.1.1           ## conda install openbabel -c openbabel
rdkit >= 2022.03.2           ## conda install rdkit -c rdkit
func-timeout >= 4.3.5        ## pip install func-timeout


After installing gnina and ADFR, please edit the PATH for these two softwares in 3D-MCTS.py:

'''
GNINA = '/home/hongyan/software/gnina'
ADFR = '/home/hongyan/software/ADFR/bin'
'''

```


## How to sample molecules for a specific protein

### Step1: Prepare several files needed

file 1. The protein structure file, ```pdb``` format. (Without ligand atoms)

file 2. The ligand file, ```sdf``` format. (To determine the position of binding site.)

file 3. The pocket file, ```pdb``` format. (To speed the calculation. It can be replaced by file 1.)

### Step2: Prepare the fragment library

We provide a fragment library comes from small molecule drugs: ```frags/fragment.txt```. We recommend users to modify it according specific needs.
We also provide a script to help users transform their customized fragments to building blocks needed by 3D-MCTS using ```prepare_building_blocks.py```:

```
python prepare_building_blocks.py --frag customized_frags.smi --o customized_building_blocks.smi
```

### Step3: Prepare the initial fragment

We provide three initial fragments in the directory ```init/```. Users can prepare other starting fragments (```pdbqt format```) according to their needs.

### Step4: Run the Code.

```
python 3D-MCTS.py --num_sims 100000 --ligand ./ligand.sdf --protein ./protein.pdb --pocket ./pocket.pdb --score -7 --qed 0.3 --processor 48 --start 1
```

Molecules that meet the criteria are saved in ```record/```.

## Cite us
```
@Article{D3SC04091G,
author ="Du, Hongyan and Jiang, Dejun and Zhang, Odin and Wu, Zhenxing and Gao, Junbo and Zhang, Xujun and Wang, Xiaorui and Deng, Yafeng and Kang, Yu and Li, Dan and Pan, Peichen and Hsieh, Chang-Yu and Hou, Tingjun",
title  ="A flexible data-free framework for structure-based de novo drug design with reinforcement learning",
journal  ="Chem. Sci.",
year  ="2023",
pages  ="-",
publisher  ="The Royal Society of Chemistry",
doi  ="10.1039/D3SC04091G",
url  ="http://dx.doi.org/10.1039/D3SC04091G",
}
```
