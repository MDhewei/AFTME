from __future__ import division
import os
aftDir = os.getcwd()


# the file folder containing PDB files
PDBDir = os.path.join(aftDir, "PDB_all")
# the file folder containing txt files written ligand nonredundant pdbID for each of 233 ligand
nonredID = os.path.join(aftDir, "NonredundantLigID")
# the excel file recording all 233 ligands information
ligInfo = os.path.join(aftDir, "lig_233.xlsx")
# the folder where distance matrix file of the ligand lies
distFolder = os.path.join(aftDir, "DistanceMatrix_use")
# the folder where clustered heatmap of the ligand lies
heatFolder = os.path.join(aftDir, "heatmaps")
# the folder where clustering results txt recording atoms of function group and motif binding lies
statisFolder = os.path.join(aftDir, "statistic")
# the folder where excel files recording all 233 ligands clustering quality lie
scoreFolder = os.path.join(aftDir, "scores")
# the excel file recording inherent atom composition for each predefined function group of each ligand
atomsLigand = os.path.join(aftDir, "atoms_233.xlsx")
# the folder placing excel files recording the function group names' assignment to the clustering number of ligand atoms
assignFolder = os.path.join(aftDir, "assign")
# the folder placing plots that display motif features for all ligands
motifFolder = os.path.join(aftDir, "motif_analyse")
# the excel file recording the frequencies of unique function group appearing in different ligands and the corresponding
# ligands, the binding motif cluster number
FG_inLigand = os.path.join(aftDir, "fg_inLig.xlsx")
# the folder placing plots that display motif Reusability for 3 FGs of ATP
atpFolder = os.path.join(aftDir, "picForATP")
# the folder placing excel files recording the reusability function groups which count >= 3
fg3mFolder = os.path.join(aftDir, "fg_reuse")


def ChangePDBCodeToList(LigandFile):
    '''this function is to obtain list of pdbcode_chainID entries after reading the ligand.txt
       o input: (1)LigandFile: the location of ligand file recording non-redundant pdbID entries
       o output: return a list of pdbcode_chainID which is shaped: ['pdbID_chainID', ...]
    '''

    PDBList = []
    with open(LigandFile, 'rt') as f:
        data = f.readlines()
        for ln in data:
            PDBList.extend(ln.strip().split(';')[:-1])
    return PDBList


'''the dict describe the chemical property of the atom in all kinds of amino acid in proteins,
different digit represent different chemical property of the atom:
1:hydrophibic
2:acceptor
3:donor
4:hydrophobic
5:aromatic
6:neutral
'''
atom_type = {'ALA': {'O': 2, 'C': 6, 'CA': 6, 'N': 3, 'CB': 4},
             'GLY': {'O': 2, 'C': 6, 'CA': 6, 'N': 3},
             'PRO': {'O': 2, 'C': 6, 'CA': 6, 'CB': 4, 'N': 3, 'CG': 4, 'CD': 4},
             'ASN': {'O': 2, 'C': 6, 'CA': 6, 'N': 3, 'CB': 4, 'CG': 6, 'OD1': 2, 'ND2': 3},
             'ASP': {'O': 2, 'C': 6, 'CA': 6, 'N': 3, 'CB': 4, 'CG': 6, 'OD1': 2, 'OD2': 1},
             'PHE': {'O': 2, 'C': 6, 'CA': 6, 'N': 3, 'CB': 4, 'CG': 5, 'CD1': 5, 'CD2': 5, 'CE1': 5, 'CE2': 5,
                     'CZ': 5},
             'LYS': {'O': 2, 'C': 6, 'CA': 6, 'N': 3, 'CG': 4, 'CD': 4, 'CE': 6, 'NZ': 3},
             'ILE': {'O': 2, 'C': 6, 'CA': 6, 'N': 3, 'CB': 4, 'CG1': 4, 'CG2': 4, 'CD1': 4},
             'LEU': {'O': 2, 'C': 6, 'CA': 6, 'N': 3, 'CB': 4, 'CG': 4, 'CD1': 4, 'CD2': 4},
             'ARG': {'O': 2, 'C': 6, 'CA': 6, 'N': 3, 'CB': 4, 'CG': 4, 'CD': 6, 'NE': 3, 'CZ': 4, 'NH1': 3, 'NH2': 1},
             'CYS': {'O': 2, 'C': 6, 'CA': 6, 'N': 3, 'CB': 6, 'SG': 6},
             'MET': {'O': 2, 'C': 6, 'CA': 6, 'N': 3, 'CB': 4, 'CG': 6, 'SD': 6, 'CE': 6},
             'THR': {'O': 2, 'C': 6, 'CA': 6, 'N': 3, 'CB': 6, 'CG2': 4, 'OG1': 1},
             'TYR': {'O': 2, 'C': 6, 'CA': 6, 'N': 3, 'CB': 4, 'CG': 5, 'CD1': 5, 'CD2': 5, 'CE1': 5, 'CE2': 5, 'CZ': 5,
                     'OH': 1},
             'HIS': {'O': 2, 'C': 6, 'CA': 6, 'N': 3, 'CB': 4, 'CG': 5, 'ND1': 3, 'CD2': 5, 'CE1': 5, 'NE2': 2},
             'VAL': {'O': 2, 'C': 6, 'CA': 6, 'N': 3, 'CB': 4, 'CG1': 4, 'CG2': 4},
             'SER': {'O': 2, 'C': 6, 'CA': 6, 'N': 3, 'CB': 4, 'OG': 1},
             'GLU': {'O': 2, 'C': 6, 'CA': 6, 'N': 3, 'CB': 4, 'CG': 3, 'CD': 6, 'OE1': 2, 'OE2': 1},
             'GLN': {'O': 2, 'C': 6, 'CA': 6, 'N': 3, 'CB': 4, 'CG': 3, 'CD': 6, 'OE1': 2, 'NE2': 3},
             'TRP': {'O': 2, 'C': 6, 'CA': 6, 'N': 3, 'CB': 4, 'CG': 6, 'CD1': 4, 'CD2': 5, 'NE1': 3, 'CE2': 5,
                     'CE3': 5, 'CZ2': 5, 'CZ3': 5, 'CH2': 5}}

prop_dic = {1: 'hydrophilic', 2: 'acceptor', 3: 'donor', 4: 'hydrophobic', 5: 'aromatic', 6: 'neutral'}

