import random
import math
import hashlib
import argparse
from rdkit import Chem as ch
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolAlign
from rdkit.Chem import rdMolTransforms
import numpy as np
import os
from openbabel.pybel import *
import time
import multiprocessing, traceback
from scipy.spatial.distance import cdist
from func_timeout import func_set_timeout
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import QED
from rdkit.Chem import Descriptors

# Please replace the PATH for GNINA and ADFR first.
GNINA = '/home/hongyan/software/gnina'
ADFR = '/home/hongyan/software/ADFR/bin'

parser = argparse.ArgumentParser(description='3D-MCTS code, for molecular generation')
parser.add_argument('--num_sims', action="store", type=int, default=1000000, help='Number of simulation steps')
parser.add_argument('--ligand', action="store", type=str, help='sdf file to determine the position of pocket',
                    default='ligand.sdf')
parser.add_argument('--protein', action="store", type=str, help='protein, PDB format',
                    default='protein.pdb')
parser.add_argument('--pocket', action="store", type=str, help='pocket, PDB format',
                    default='pocket.pdb')
parser.add_argument('--score', action="store", type=float, help='threshold for vina score', default=-7)
parser.add_argument('--qed', action="store", type=float, help='threshold for qed', default=0.3)
parser.add_argument('--processor', action="store", type=int, help='number of processor for multiprocessing', default=48)
parser.add_argument('--start', action="store", type=str, help='start fragment', default='1')
parser.add_argument('--frag_lib', action="store", type=str, help='fragment library', default='frags/fragment.txt')
parser.add_argument('--gnina', action="store", type=str, help='the path for GNINA',
                    default='/home/hongyan/software/gnina')
parser.add_argument('--adfr', action="store", type=str, help='the path for adfr',
                    default='/home/hongyan/software/ADFR/bin')

args = parser.parse_args()

os.system('mkdir node vinaa record tmp init_pose unique')
os.system(rf'{ADFR}/prepare_receptor -r {args.protein} -o pro_test.pdbqt')

# Used to balance exploration and exploitation
SCALAR = 1 / (2 * math.sqrt(2.0))
FRAGMENT_LIB = []
with open(f'{args.frag_lib}','r') as f:
    lines = f.readlines()
    for line in lines:
        FRAGMENT_LIB.append(line.strip())

GNINA = args.gnina
ADFR = args.adfr

# Rules for reassembly of fragments
Iso_dic = {
    1: [3, 5, 10],  ###
    2: [],
    3: [1, 6, 4, 8, 13, 14, 15, 16],
    4: [3, 5, 11, ],  ####
    5: [1, 6, 4, 8, 12, 14, 16, 13, 15],
    6: [13, 14, 15, 16, 3, 5, 10],
    7: [],
    8: [3, 5, 11, 9, 10, 13, 14, 15, 16],
    9: [8, 13, 14, 15, 16],
    10: [1, 6, 8, 13, 14, 15, 16],
    11: [4, 8, 13, 14, 15, 16],
    12: [5],
    13: [3, 5, 6, 8, 9, 10, 11, 14, 15, 16],
    14: [3, 5, 6, 8, 9, 10, 11, 13, 14, 15, 16],
    15: [3, 5, 6, 8, 9, 10, 11, 13, 14, 16],
    16: [3, 5, 6, 8, 9, 10, 11, 13, 14, 15, 16]
}
pro = 'pro_test.pdbqt'
# init = ''
NODE_ID = 1
RECORD = 1
SCORES = []
GOOD_SCORES = []
NO_QED = []
NO_QED_GOOD = []
mol_dic = []
mol = Chem.MolFromPDBFile(rf"{args.pocket}")
writer = Chem.SDWriter(rf'./pocket.sdf')
writer.write(mol)
writer.close()
POCKET = Chem.SDMolSupplier('pocket.sdf')[0]
P = POCKET.GetConformer()
POCKET_CORS = [P.GetPositions()[i] for i in range(len(POCKET.GetAtoms()))]



def dock(init):

    '''
    Generate the conformation of the starting fragment.
    '''
    mymol = list(readfile('sdf', rf'init/{init}.sdf'))[0]
    mymol.write("pdbqt", rf'init/{init}.pdbqt', overwrite=True)
    os.system(
        rf'{GNINA}/gnina --receptor {pro} --ligand init/{init}.pdbqt --autobox_ligand {args.ligand} --cnn_scoring=none --out ./init_pose/{init}.pdbqt')

    modes = open(rf'./init_pose/{init}.pdbqt', 'r').read().split('ENDMDL\n')
    valids = []
    for i in range(9):

        mode = open(rf'./init_pose/{init}_{i}.pdbqt', 'w')
        score = float(modes[i].split('\n')[1].split()[2].replace('REMARK', ''))
        mode.write(modes[i] + 'ENDMDL\n')
        mode.close()

        os.system(rf'obabel ./init_pose/{init}_{i}.pdbqt -O ./init_pose/{init}_{i}.sdf -D -B -R -C')

        f = open(rf'./init_pose/{init}_{i}.sdf', 'r').readlines()
        f_ = open(rf'./init_pose/{init}_{i}_repair.sdf', 'w')
        for line in f[:4]:
            f_.write(line)
        for j in range(4, len(f)):
            if f[j].startswith(' ') and len(f[j]) >= 50:
                new_line = f[j][:50] + '0' + f[j][51:]
                f_.write(new_line)
            else:
                br = j
                break
        for j in range(br, len(f)):
            f_.write(f[j])
        f_.close()

        try:
            mol_h = Chem.AddHs(Chem.SDMolSupplier(rf'./init_pose/{init}_{i}_repair.sdf')[0], addCoords=True)
            writer = Chem.SDWriter(rf'./init_pose/{init}_{i}_repair_H.sdf')
            writer.write(mol_h)
            writer.close()
            valids.append([init, i, score])
        except:
            traceback.print_exc()
            continue
    return valids


def assign_iso(frag_id, Del_plus, start=False):
    '''
    Isotope labeling of fragment attachment sites
    '''
    if start == False:
        lig_name = rf'{frag_id}_{Del_plus}'
        mol_h = Chem.AddHs(Chem.SDMolSupplier(rf'./vinaa/{lig_name}_repair.sdf')[0], addCoords=True)
    else:
        lig_name = rf'{frag_id}_{Del_plus}'
        mol_h = Chem.SDMolSupplier(rf'./init_pose/{lig_name}_repair_H.sdf', removeHs=False)[0]
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
    if start == False:
        writer = Chem.SDWriter(rf'./tmp/{lig_name}.sdf')
        writer.write(mol_h)
        writer.close()
    else:
        writer = Chem.SDWriter(rf'./start.sdf')
        writer.write(mol_h)
        writer.close()


def get_id_bysymbol(combo, symbol):
    for at in combo.GetAtoms():
        if at.GetSymbol() == symbol:
            return at.GetIdx()


def get_neiid_bysymbol(combo, symbol):
    for at in combo.GetAtoms():
        if at.GetSymbol() == symbol:
            at_nei = at.GetNeighbors()[0]
            return at_nei.GetIdx()


def get_neiid_byisotope(combo, symbol):
    for at in combo.GetAtoms():
        if at.GetIsotope() != 0:
            at_nei = at.GetNeighbors()[0]
            return at_nei.GetIdx()


def frag_avail(h1):

    lib_avail = []
    iso_avail = Iso_dic[h1]
    for i in iso_avail:
        for smi in FRAGMENT_LIB:
            if rf'[{i}*]' in smi:
                lib_avail.append(smi)

    return lib_avail


def mol_connect(par, h1, frag, frag_id, smile):
    '''
    Connect new molecular fragments
    '''
    @func_set_timeout(5)
    def mol_connect2(par, h1, frag, frag_id, smile):
        for atom in par.GetAtoms():
            atom.SetIntProp('FromPar', 1)
            if atom.GetIsotope() != 0 and atom.GetIdx() == h1:
                Iso_type = atom.GetIsotope()
                atom.GetNeighbors()[0].SetIntProp('Nei', 1)
                atom.SetAtomicNum(37)
            if atom.GetIsotope() != 0 and atom.GetIdx() != h1:
                atom.SetAtomicNum(1)

        Iso_types = Iso_dic[Iso_type]

        for atom in frag.GetAtoms():
            if atom.GetIsotope() in Iso_types:
                atom.SetAtomicNum(87)
                atom.GetNeighbors()[0].SetIntProp('Nei', 1)
                break

        for atom in frag.GetAtoms():
            atom.SetIntProp('FromPar', 0)
            if atom.GetIsotope() != 0 and atom.GetAtomicNum() == 0:
                atom.SetAtomicNum(1)

        combo = ch.CombineMols(par, frag)
        Rb_neiid = get_neiid_bysymbol(combo, 'Rb')
        Fr_neiid = get_neiid_bysymbol(combo, 'Fr')
        edcombo = ch.EditableMol(combo)
        edcombo.AddBond(Rb_neiid, Fr_neiid, order=Chem.rdchem.BondType.SINGLE)

        Rb_index = get_id_bysymbol(combo, 'Rb')
        edcombo.RemoveAtom(Rb_index)
        back = edcombo.GetMol()

        Fr_index = get_id_bysymbol(back, 'Fr')
        edcombo = Chem.EditableMol(back)
        edcombo.RemoveAtom(Fr_index)
        back = edcombo.GetMol()

        Dihedral = []
        for at in back.GetAtoms():
            if at.GetProp('FromPar') == '1' and at.HasProp('Nei') == 1:
                for nei in at.GetNeighbors():
                    if nei.GetProp('FromPar') == '1':
                        Dihedral.append(str(nei.GetIdx()))
                        Dihedral.append(str(at.GetIdx()))
                        break
                break

        for at in back.GetAtoms():
            if at.GetProp('FromPar') == '0' and at.HasProp('Nei') == 1:
                for nei in at.GetNeighbors():
                    if nei.GetProp('FromPar') == '0':
                        Dihedral.append(str(at.GetIdx()))
                        Dihedral.append(str(nei.GetIdx()))
                        break
                break

        par_back = par
        par_back.GetAtomWithIdx(get_id_bysymbol(par_back, 'Rb')).SetAtomicNum(1)
        par_back_rm_H = Chem.RemoveHs(par_back)
        for at in par_back_rm_H.GetAtoms():
            at.SetIsotope(0)
        par_back_rm_H = Chem.RemoveHs(par_back_rm_H)
        par_index = par_back.GetSubstructMatch(par_back_rm_H)
        combo_index = back.GetSubstructMatch(par_back_rm_H)
        try:
            cmap = {combo_index[j]: par_back.GetConformer().GetAtomPosition(par_index[j]) for j in
                    range(len(par_index))}
        except:
            traceback.print_exc()
            print(
                [par_index, combo_index, ch.MolToSmiles(par_back), ch.MolToSmiles(back), ch.MolToSmiles(par_back_rm_H)]
            )
        cids = AllChem.EmbedMultipleConfs(back, numConfs=30, coordMap=cmap, maxAttempts=1000, numThreads=4,
                                          randomSeed=1)
        tag = 1
        if len(cids) > 0:

            rms_min = 10
            for i in range(len(cids)):

                rms = rdMolAlign.AlignMol(back, par_back, prbCid=i, atomMap=list(zip(combo_index, par_index)))
                if rms < 2:
                    if rms < rms_min:
                        rms_min = rms
                        writer = Chem.SDWriter(rf'tmp/{frag_id}.sdf')
                        writer.SetProps(['DihAtoms', 'DihDeg'])
                        back.SetProp('DihAtoms', ','.join(Dihedral))
                        Deg = rdMolTransforms.GetDihedralDeg(back.GetConformer(id=i), int(Dihedral[0]),
                                                             int(Dihedral[1]),
                                                             int(Dihedral[2]), int(Dihedral[3]))
                        back.SetProp('DihDeg', str(Deg))
                        writer.write(back, confId=i)
                        writer.close()
                        tag = 0
        return [frag_id, tag, smile]

    try:
        a = mol_connect2(par, h1, frag, frag_id, smile)
        return a
    except:
        return [frag_id, 1, smile]


def mol_rotate(frag_id, Del_plus, smile):
    '''
    Rotate newly introduced dihedral angles to get multiple molecular conformations
    '''
    mol = ch.SDMolSupplier(rf'tmp/{frag_id}.sdf', removeHs=False)[0]
    Dihedral_atoms = mol.GetProp("DihAtoms").split(',')
    Deg = float(mol.GetProp("DihDeg")) + Del_plus
    Comf = mol.GetConformer()
    rdMolTransforms.SetDihedralDeg(Comf, int(Dihedral_atoms[0]), int(Dihedral_atoms[1]), int(Dihedral_atoms[2]),
                                   int(Dihedral_atoms[3]), Deg)
    writer = Chem.SDWriter(rf'tmp/{frag_id}_{Del_plus}.sdf')
    writer.write(mol)
    writer.close()


def score0(frag_id, Del_plus, smile):
    '''
    Determine whether newly connected fragments can cause collisions between atoms
    '''
    time1 = time.time()
    lig = ch.SDMolSupplier(rf'tmp/{frag_id}_{Del_plus}.sdf', removeHs=False)[0]
    for atom in lig.GetAtoms():
        if atom.GetAtomicNum() == 37:
            atom.SetAtomicNum(1)
    lig_name = rf'{frag_id}_{Del_plus}'
    L = lig.GetConformer()
    lig_cors = [L.GetPositions()[i] for i in range(len(lig.GetAtoms()))]

    dis_matrix = cdist(POCKET_CORS, lig_cors, metric='euclidean')
    pt = Chem.GetPeriodicTable()
    radis_matrix = np.zeros((len(POCKET.GetAtoms()), len(lig.GetAtoms())))
    for i in range(len(POCKET.GetAtoms())):
        at1 = pt.GetRcovalent(POCKET.GetAtomWithIdx(i).GetAtomicNum())
        for j in range(len(lig.GetAtoms())):
            at2 = pt.GetRcovalent(lig.GetAtomWithIdx(j).GetAtomicNum())
            radis_matrix[i][j] = at1 + at2
    judge = (radis_matrix > dis_matrix).sum()

    lig_intra_matrix = cdist(lig_cors, lig_cors, metric='euclidean')
    covalent = np.ones((len(lig.GetAtoms()), len(lig.GetAtoms())))
    bonds = lig.GetBonds()
    for bd in bonds:
        idx1 = bd.GetBeginAtomIdx()
        idx2 = bd.GetEndAtomIdx()
        covalent[idx1][idx2] = 0
        covalent[idx2][idx1] = 0
        covalent[idx1][idx1] = 0
        covalent[idx2][idx2] = 0
    lig_intra_radis = np.zeros((len(lig.GetAtoms()), len(lig.GetAtoms())))
    for i in range(len(lig.GetAtoms())):
        at1 = pt.GetRcovalent(lig.GetAtomWithIdx(i).GetAtomicNum())
        for j in range(len(lig.GetAtoms())):
            at2 = pt.GetRcovalent(lig.GetAtomWithIdx(j).GetAtomicNum())
            lig_intra_radis[i][j] = at1 + at2
    judge2 = ((lig_intra_matrix < lig_intra_radis) * covalent).sum()
    # print('score0:',time.time()-time1)

    return [frag_id, Del_plus, judge + judge2, smile]


def score1(frag_id, Del_plus, smile):
    '''
    Prepare files for scoring.
    '''
    time1 = time.time()
    pro = 'pro_test.pdbqt'
    mol = ch.SDMolSupplier(rf'tmp/{frag_id}_{Del_plus}.sdf', removeHs=False)[0]
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 37:
            atom.SetAtomicNum(1)
    lig_name = rf'{frag_id}_{Del_plus}'
    writer = Chem.SDWriter(rf'vinaa/{lig_name}_score.sdf')
    writer.write(mol)
    writer.close()
    mymol = list(readfile('sdf', rf'vinaa/{lig_name}_score.sdf'))[0]
    mymol.write("pdbqt", 'vinaa/' + lig_name + '.pdbqt', overwrite=True)


def score3(frag_id, Del_plus, smile):
    '''
    Obtain the binding affinity of the molecule
    '''
    time1 = time.time()
    pro = 'pro_test.pdbqt'
    lig_name = rf'{frag_id}_{Del_plus}'

    try:
        score = float(os.popen(
            rf'{GNINA}/gnina --minimize --receptor {pro} --ligand vinaa/{lig_name}.pdbqt --cnn_scoring=none --out ./vinaa/{lig_name}_minimize.sdf').read().split(
            'Affinity: ')[1].split(' (kcal/mol)')[0].split()[0])
    except:
        score = 0

    return [frag_id, Del_plus, score, smile]


def repair(frag_id, Del_plus, score, smile):
    '''
    Repair the chemical valence of atoms after minimization
    '''
    lig_name = rf'{frag_id}_{Del_plus}'
    f = open(rf'./vinaa/{lig_name}_minimize.sdf', 'r').readlines()
    f_ = open(rf'./vinaa/{lig_name}_repair.sdf', 'w')
    for line in f[:4]:
        f_.write(line)
    for i in range(4, len(f)):
        if f[i].startswith(' ') and len(f[i]) >= 50:
            new_line = f[i][:50] + '0' + f[i][51:]
            f_.write(new_line)
        else:
            br = i
            break
    for i in range(br, len(f)):
        f_.write(f[i])
    f_.close()

    try:
        assign_iso(frag_id, Del_plus)
    except:
        mol_parent = Chem.SDMolSupplier(rf'./tmp/{lig_name}.sdf', removeHs=False)[0]
        for at in mol_parent.GetAtoms():
            if at.GetIsotope() != 0:
                at.SetIsotope(0)
        writer = Chem.SDWriter(rf'./tmp/{lig_name}.sdf')
        writer.write(mol_parent)
        writer.close()
        traceback.print_exc()


def qed_score(frag_id, Del_plus, score, smile):
    '''
    Obtain the drug-like properties of a molecule
    '''
    lig_name = rf'{frag_id}_{Del_plus}'
    mol = Chem.SDMolSupplier(rf'./tmp/{lig_name}.sdf', removeHs=False)[0]
    qed = QED.qed(mol)
    return [frag_id, Del_plus, score, smile, qed]


def roulette(select_list):
    '''
    roulette algorithm
    '''
    sum_val = sum(select_list)
    random_val = random.random()
    probability = 0
    if sum_val != 0:
        for i in range(len(select_list)):
            probability += select_list[i] / sum_val
            if probability >= random_val:
                return i
            else:
                continue
    else:
        return random.choice(range(len(select_list)))


class State():
    '''
    The status of the node. There are two types of nodes in 3D-MCTS, one is used to determine the connection position of fragments (type 1),
    and the other is used to determine the type and conformation of fragments (type 0).
    '''
    def __init__(self, state_type=0, sdf='start.sdf', h1=None, Frag_Deg=None, frag=None, sco=-3.55479, ter_tag=None,
                 sta=[], choices=[]):
        self.type = state_type
        self.sdf = sdf

        if self.type == 0:

            self.Frag_Deg = Frag_Deg
            self.h1s_avail = []
            mol = ch.SDMolSupplier(self.sdf, removeHs=False)[0]
            if self.sdf == 'start.sdf':
                os.system(rf'cp start.sdf state0.sdf')
            self.score = sco
            for atom in mol.GetAtoms():
                if atom.GetIsotope() != 0:
                    self.h1s_avail.append(atom.GetIdx())
            self.states = sta + [self.score]
            self.choices = choices + [self.Frag_Deg]

        elif self.type == 1:

            self.h1 = h1

            self.score = 0
            self.states = sta + [self.score]
            self.choices = choices + [self.h1]
            self.ter_tag = 0

    def next_states(self):
        '''
        Get all possibilities for the next state of the current state.

        '''
        if self.type == 0:
            pass

        elif self.type == 1:
            mol = ch.SDMolSupplier(self.sdf, removeHs=False)[0]
            h1_isotop = mol.GetAtomWithIdx(self.h1).GetIsotope()
            try:
                self.frags_avail = frag_avail(h1_isotop)
            except:
                print(self.sdf, self.h1, h1_isotop)

            frags = []
            mols = []
            h1s = []
            ids = []
            j = 0
            time1 = time.time()
            for i in self.frags_avail:
                frags.append(ch.AddHs(ch.MolFromSmiles(i)))
                mols.append(mol)
                h1s.append(self.h1)
                j += 1
                ids.append(j)
            # fragment connection
            pool = multiprocessing.Pool(args.processor)
            ids_tags_smiles = pool.starmap(mol_connect, zip(mols, h1s, frags, ids, self.frags_avail))
            pool.close()
            pool.join()
            # print(time.time()-time1)
            time1 = time.time()
            ids_degs_smiles = []
            for id_tag_smile in ids_tags_smiles:
                if id_tag_smile[1] == 0:
                    for deg in np.arange(0, 360, 15):
                        ids_degs_smiles.append([id_tag_smile[0], deg, id_tag_smile[2]])

            # Rotate dihedral angles to obtain multiple conformations
            pool = multiprocessing.Pool(args.processor)
            pool.starmap(mol_rotate, ids_degs_smiles)
            pool.close()
            pool.join()
            # print(time.time() - time1)
            # nextmove = random.choice([x for x in self.frags_avail])

            time1 = time.time()
            pool = multiprocessing.Pool(args.processor)
            ids_degs_judge_smiles = pool.starmap(score0, ids_degs_smiles)
            pool.close()
            pool.join()
            # print(time.time() - time1)

            ids_degs_smiles = [[i[0], i[1], i[3]] for i in ids_degs_judge_smiles if i[2] == 0]
            if len(ids_degs_smiles) > 0:
                pool = multiprocessing.Pool(args.processor)
                pool.starmap(score1, ids_degs_smiles)
                pool.close()
                pool.join()

                pool = multiprocessing.Pool(args.processor)
                ids_degs_scores_smiles = pool.starmap(score3, ids_degs_smiles)
                pool.close()
                pool.join()
                ids_degs_scores = [i for i in ids_degs_scores_smiles if i[2] <= (self.states[-2] + 0.5)]

                # Sort multiple conformations according to binding affinity
                ids_degs_scores_smiles = ids_degs_scores
                ids_degs_scores_smiles_selected = []
                ids_degs_scores_sorted = sorted(ids_degs_scores, key=lambda x: x[2])
                ids_degs_scores_smiles_dic = {}
                for i in ids_degs_scores_sorted:
                    try:
                        if len(ids_degs_scores_smiles_dic[i[0]]) < 1:
                            ids_degs_scores_smiles_dic[i[0]].append(i)
                        else:
                            continue
                    except:
                        ids_degs_scores_smiles_dic[i[0]] = [i]
                for value in ids_degs_scores_smiles_dic.values():
                    ids_degs_scores_smiles_selected.extend(value)


                # Repair the structure after minimization
                pool = multiprocessing.Pool(args.processor)
                pool.starmap(repair, ids_degs_scores_smiles_selected)
                pool.close()
                pool.join()

                # Save the molecules that match the criteria
                for i in ids_degs_scores_smiles_selected:
                    if i[2] < (args.score):
                        mol = ch.SDMolSupplier(rf'tmp/{i[0]}_{i[1]}.sdf', removeHs=False)[0]
                        global RECORD
                        for at in mol.GetAtoms():
                            at.SetIsotope(0)
                        writer = Chem.SDWriter(rf'tmp/{i[0]}_{i[1]}_.sdf')
                        writer.write(mol)
                        writer.close()
                        mol = ch.SDMolSupplier(rf'tmp/{i[0]}_{i[1]}_.sdf', removeHs=False)[0]
                        if Chem.MolToSmiles(mol) not in mol_dic:
                            mol_dic.append(Chem.MolToSmiles(mol))
                            qed = QED.qed(mol)
                            NO_QED.append(i[2])
                            if qed > args.qed:
                                writer = Chem.SDWriter(rf'record/record_{RECORD}.sdf')
                                writer.SetProps(['Score'])
                                mol.SetProp('Score', str(i[2]))
                                writer.write(mol)
                                writer.close()
                                RECORD += 1
                                SCORES.append(i[2])
                                if i[2] < (args.score + 0):
                                    GOOD_SCORES.append(i[2])
                            if i[2] < (args.score + 0):
                                NO_QED_GOOD.append(i[2])

                if False:
                    self.ids_degs_scores = ids_degs_scores_smiles[:50]
                else:
                    pool = multiprocessing.Pool(args.processor)
                    ids_degs_scores = pool.starmap(qed_score, ids_degs_scores_smiles_selected)
                    pool.close()
                    pool.join()


                    self.ids_degs_scores = [i for i in ids_degs_scores if i[-1] > args.qed]

                if len(self.ids_degs_scores) == 0:
                    self.ter_tag = 1
                    print('node_end')
                    return []
                else:
                    self.ter_tag = 0
                    print('expand')
                    return self.ids_degs_scores
            else:
                self.ter_tag = 1
                return []

    def terminal(self, tag):

        # Determine whether the current state meets the termination conditions
        if self.type == 1:
            if self.ter_tag == 1:
                return True
            else:
                return False


        else:
            mol = ch.SDMolSupplier(self.sdf, removeHs=False)[0]

            if len(self.h1s_avail) == 0:
                return True

            if self.score > 0:
                self.stop_normal = 0
                return True
            else:
                try:
                    if (self.score - 1.5) > self.states[-3]:
                        return True
                    else:
                        if rdMolDescriptors.CalcNumLipinskiHBA(
                                mol) > 9 or rdMolDescriptors.CalcNumLipinskiHBD(mol) > 4 or Descriptors.MolWt(
                            mol) > 500:
                            self.stop_normal = 1
                            return True
                        else:
                            return False
                except:
                    if rdMolDescriptors.CalcNumLipinskiHBA(
                            mol) > 9 or rdMolDescriptors.CalcNumLipinskiHBD(mol) > 4 or Descriptors.MolWt(mol) > 500:
                        self.stop_normal = 1
                        return True
                    else:
                        return False

    def reward(self):
        r1 = 0
        r2 = max((np.mean(SCORES) - self.score2) / abs(np.mean(SCORES + [self.score2])), 0)
        if self.score2 < np.mean(SCORES):
            SCORES.append(self.score2)
        print(self.score, self.score2)
        print(r1, r2)
        return r1, r2

    def __hash__(self):
        if self.type == 0:
            return int(hashlib.md5(str(self.Frag_Deg).encode('utf-8')).hexdigest(), 16)
        elif self.type == 1:
            return int(hashlib.md5(str(self.h1).encode('utf-8')).hexdigest(), 16)
    def __eq__(self, other):
        if hash(self) == hash(other):
            return True
        return False

class Node():
    def __init__(self, state, parent=None, node_id=1, best_score=-3.55479, reward=0, qed=1):
        self.visits = 0
        self.reward = reward
        self.state = state
        self.children = []
        self.parent = parent
        self.best_score = best_score
        self.longest_path = 0
        self.id = node_id
        self.qed = qed

    def add_child(self, child_state, node_id, bestscore, qed=1):
        sdf_name = child_state.sdf
        os.system(rf'cp {sdf_name} node/{node_id}.sdf')
        child_state.sdf = rf'node/{node_id}.sdf'
        child = Node(child_state, node_id=node_id, parent=self, best_score=bestscore, qed=qed)
        self.children.append(child)

    def update(self, reward):
        self.reward += reward
        self.visits += 1

    def fully_expanded(self, num_moves_lambda):
        if self.state.type == 0:
            num_moves = len(self.state.h1s_avail)
        elif self.state.type == 1:
            num_moves = len(self.state.ids_degs_scores)
        if num_moves_lambda != None:
            num_moves = num_moves_lambda(self)
        if len(self.children) == num_moves:
            return True
        return False

    def __repr__(self):
        s = "Node; children: %d; visits: %d; reward: %f" % (len(self.children), self.visits, self.reward)
        return s


def UCTSEARCH(budget, root, start_score=0, num_moves_lambda=None):
    # Begin the MCTS
    start = time.time()
    global GOOD_SCORES
    global NO_QED_GOOD

    for iter in range(int(budget)):
        front = TREEPOLICY(root, start_score, num_moves_lambda)
        BACKUP2(front)


def TREEPOLICY(node, start_score, num_moves_lambda):
    # Choose whether to expand the node based on the status of the current node
    while node.state.terminal('state1.sdf') == False:
        if len(node.children) == 0:
            node = EXPAND(node, start_score)
        else:
            node = BESTCHILD(node, SCALAR, start_score)
    return node


def EXPAND(node, start_score):

    # Get the children of a node and add them to the tree
    global NODE_ID
    if node.state.type == 0:
        for nextmove in node.state.h1s_avail:
            next = State(state_type=1, sdf=rf'{node.state.sdf}', h1=nextmove, sta=node.state.states)
            NODE_ID += 1
            node.add_child(next, node_id=NODE_ID, bestscore=start_score)
        return node.children[-1]
    elif node.state.type == 1:
        new_states = node.state.next_states()
        if len(new_states) == 0:
            return node
        else:
            scores = []
            for nextmove in new_states:
                os.system(rf'cp tmp/{nextmove[0]}_{nextmove[1]}.sdf state0.sdf')
                next = State(state_type=0, sdf=rf'state0.sdf', Frag_Deg=nextmove[:2], sco=nextmove[2],
                             sta=node.state.states)
                NODE_ID += 1
                best_score = min(start_score, nextmove[2])
                scores.append(abs(nextmove[4]))
                node.add_child(next, node_id=NODE_ID, bestscore=best_score, qed=abs(nextmove[4]))
            return node.children[roulette(scores)]

def BESTCHILD(node, scalar, start_score):

    # Select child nodes based on the node's UCB
    scores = []
    scores2 = []
    for c in node.children:

        exploit = start_score - c.best_score
        explore = math.sqrt(2.0 * math.log(node.visits + 0.000001) / float(c.visits + 0.000001))
        score = exploit + scalar * explore
        scores.append(score)
        scores2.append(c.qed)
    if True:
        idx1 = roulette(scores)
        idx2 = roulette(scores2)
        idx = random.choice([idx1, idx2])
    else:
        idx = random.choice(range(len(scores)))
    return node.children[idx]


def DEFAULTPOLICY(node):
    state = node.state
    num_states = 0

    while state.terminal('state1.sdf') == False:
        state = state.next_state()
        num_states += 1
    if state.type == 1:
        if num_states != 0:
            num_states -= 1
    num_nodes = len(state.states) - num_states
    print(state.type)
    return state.states, num_nodes, num_states


def BACKUP2(node):
    '''
    After a search is completed, backpropagation is performed to update the number of visits to the node.
    '''
    parent_node = node
    while parent_node != None:
        parent_node.visits += 1
        if len(parent_node.children) == 0:
            x = parent_node
            parent_node = node.parent
            son_node = x
        else:
            if parent_node.best_score > son_node.best_score:
                parent_node.best_score = son_node.best_score
            x = parent_node
            parent_node = parent_node.parent
            son_node = x


def BACKUP(node, states, num_nodes, num_states):
    print(states)
    i = 1
    if node.longest_path == 0:
        node.longest_path = len(states)
    while node != None:
        node.visits += 1
        best_score = min(states[num_nodes - i:])
        i += 1
        if best_score not in SCORES:
            if best_score < node.best_score:
                node.best_score = best_score
                reward = max(-3.55479 - best_score, 0)
            else:
                reward = 0
            if best_score < np.mean(SCORES):
                SCORES.append(best_score)
        else:
            if best_score < node.best_score:
                node.best_score = best_score
                reward = max(-3.55479 - best_score, 0)
            else:
                reward = 0
        node.reward += reward
        node = node.parent
    return


if __name__ == "__main__":

    results = []
    i = args.start
    # Obtain binding modes of starting fragments using molecular docking
    results.extend(dock(i))
    # Choose the best docking conformation
    results.sort(key=lambda results: results[2])
    for i in range(len(results)):
        # Label starting fragments with isotopes for new fragment ligation
        assign_iso(results[i][0], results[i][1], start=True)
        current_node = Node(State(sta=[], sco=results[i][2]), best_score=results[i][2])
        result = UCTSEARCH(args.num_sims, current_node, start_score=results[i][2])
