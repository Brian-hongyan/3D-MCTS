# Convert fragments into labeled building blocks according to BRICS rules

from rdkit import Chem
import argparse

def assign_iso(smi):

    '''
    Using isotopes to label types of heavy atoms. C N O S
    '''

    mol_h = Chem.AddHs(Chem.MolFromSmiles(smi))
    for at in mol_h.GetAtoms():
        if True:
            label = 0
            if at.GetAtomicNum() == 6:
                if at.IsInRing():
                    if at.GetIsAromatic():
                        label = 16
                        neis = at.GetNeighbors()
                        for nei in neis:
                            if nei.GetAtomicNum() == 7 or 8 or 16:
                                label = 14
                    else:
                        label = 15
                        neis = at.GetNeighbors()
                        for nei in neis:
                            if nei.GetAtomicNum() == 7 or 8 or 16:
                                label = 13
                else:
                    neis = at.GetNeighbors()
                    label = 8
                    for nei in neis:
                        bond = mol_h.GetBondBetweenAtoms(at.GetIdx(), nei.GetIdx())
                        if bond.GetBondType() != Chem.rdchem.BondType.SINGLE:
                            if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE and nei.GetAtomicNum() == 8:
                                label = 6
                                break
                            else:
                                label = 0
                                break
            if at.GetAtomicNum() == 7:
                label = 0
                if at.GetIsAromatic():
                    label = 9
                else:
                    if at.IsInRing():
                        neis = at.GetNeighbors()
                        for nei in neis:
                            if nei.GetAtomicNum() == 6:
                                neis2 = nei.GetNeighbors()
                                for nei2 in neis2:
                                    if nei2.GetAtomicNum() == 8 and mol_h.GetBondBetweenAtoms(nei.GetIdx(),
                                                                                              nei2.GetIdx()).GetBondType() == Chem.rdchem.BondType.DOUBLE:
                                        label = 10
                if label != 10:
                    if at.GetHybridization() == Chem.rdchem.HybridizationType.SP3:
                        label = 5
            if at.GetAtomicNum() == 8:
                label = 3
            if at.GetAtomicNum() == 16:
                if at.GetHybridization() == Chem.rdchem.HybridizationType.SP3:
                    label = 11
                else:
                    neis = at.GetNeighbors()
                    count = 0
                    for nei in neis:
                        bond = mol_h.GetBondBetweenAtoms(at.GetIdx(), nei.GetIdx())
                        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE and nei.GetAtomicNum() == 8:
                            count += 1
                    if count == 2:
                        label = 12
            if label != 0:
                neis = at.GetNeighbors()
                for at in neis:
                    if at.GetAtomicNum() == 1:
                        at.SetIsotope(label)
    
    mol_rm_h = Chem.RemoveHs(mol_h)

    for at in mol_rm_h.GetAtoms():
        if at.GetAtomicNum() == 1 and at.GetIsotope() != 0:
            at.SetAtomicNum(0)
            pass
    return Chem.MolToSmiles(mol_rm_h)

parser = argparse.ArgumentParser(description='To prepare user customed building blocks for 3D-MCTS.')
parser.add_argument('--frag', action="store", type=str, default='./test.smi',help='smi file that contains smiles of fragments.')
parser.add_argument('--o', action="store", type=str, default='./customized_BBs.smi', help='output file contains prepared building blocks')

args = parser.parse_args()

smis = open(args.frag,'r').readlines()
out = open(args.o,'w')
for smi in smis:
    smi = smi.split()[0]
    try:
        prepared = assign_iso(smi)
        out.write(f'{prepared}\n')
    except:
        continue

out.close()