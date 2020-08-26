import os
from Bio.PDB.PDBParser import PDBParser
import pandas as pd
from filepath import *


def GetLigand(PDBCode, Ligand, chain, l):
    '''Given a PDB code this function is to obtain the target ligand
       o input: (1) PDBCode: a valid pdb code which is named in scheme pdbcode.pdb
                (2) Ligand: the ligand id e.g.,ATP
                (3) chain: the chain ID of the target ligand
                (4) l: the number of target ligand heavy atoms
       o output: return the target ligand structure
    '''
    p = PDBParser(PERMISSIVE=1, QUIET=True)
    PDBFileName = os.path.join(PDBDir, PDBCode)
    s = p.get_structure(PDBCode, PDBFileName)
    chain = s[0][chain]
    for residue in chain.get_residues():
        residue_id = residue.get_id()
        if residue_id[0] == "H_" + Ligand:
            # filter the ligand with modified atoms
            if len(residue) == l:
                return residue


def GetFunctionAtomsSet(PDBCode, Ligand, chain, fix_lenth):
    '''this function is to extract protein atoms excluded main chain atoms within 5 Å of the ligand atoms
       atoms within R to any atoms of the ligand are determined as functional atoms
       o input: (1) PDBCode: a valid pdb code which is named in scheme pdbcode.pdb
                (2) Ligand: the ligand id e.g.,ATP
                (3) chain: the chain ID of the target ligand
                (4) fix_lenth: the number of target ligand heavy atoms
       o output: return a list of atoms within 5 Å of the target ligand
    '''
    FAS = []
    ligand = GetLigand(PDBCode, Ligand, chain, fix_lenth)
    try:
        atoms1 = ligand.get_parent().get_atoms()
        # to exclude main chain atoms of the protein
        atom_list = ["C", "N", "O", "CA", "CB"]
        # traverse the protein atoms
        for atom1 in atoms1:
            # exclude the non-protein atoms and the main chain atoms
            if atom1.get_full_id()[3][0] == ' ' and atom1.get_name() not in atom_list:
                # traverse the ligand atoms
                for atom2 in ligand:
                    # add the protein atom closer than cutoff distance 5
                    if 0 < atom1-atom2 < 5.0:
                        FAS.append(atom1)
    except AttributeError as e:
        print(e)

    return FAS


def GetAtomToLigandDistanceDict(PDBCode, Ligand, fix_lenth, Chain):
    '''this function is to obtain the distance of functional atoms (fas) to all the atoms of the ligand
       o input: (1) PDBCode: a valid pdb code which is named in scheme pdbcode.pdb
                (2) Ligand: the ligand id e.g.,ATP
                (3) fix_lenth: the number of target ligand heavy atoms
                (4) Chain: the chain ID of the target ligand
       o output: return dict recording the distance to all atoms of the ligand for each fas,
                 and the functional atom is named in scheme: "pdbID_AminoAcid+AminoAcidID_atom"
                 e.g.,{'1A0I_LYS10_CE': {'PG': 6.995315, ..., 'C5': 10.834808}
                        ...
                       '1A0I_GLU32_CG': {'PG': 18.274075, ..., 'C5': 6.933704}}
    '''
    Dic_all = {}
    # get ligand object
    ligand = GetLigand(PDBCode, Ligand, Chain, fix_lenth)
    # get fas given pdb code
    FAS = GetFunctionAtomsSet(PDBCode, Ligand, Chain, fix_lenth)
    if FAS:
        for i in sorted(FAS):
            # to label atom with "pdbID_AminoAcid+AminoAcidID_atom"
            atom_label = PDBCode[:4] + "_" + i.get_parent().get_resname() + str(i.get_parent().get_id()[1]) \
                         + "_" + i.get_name()
            DisToLigandDic = {}
            for j in range(fix_lenth):
                distance = i - ligand.get_list()[j]
                DisToLigandDic[ligand.get_list()[j].get_name()] = distance
            Dic_all[atom_label] = DisToLigandDic
    return Dic_all


def writeLigFAS(LigInfo=ligInfo, NonRedunID=nonredID, DistFolder=distFolder):
    '''this function is to obtain functional atoms distance to their binding ligand for all ligands.
       o input: (1) ligInfo: the location of file recording all 233 ligands information
                (2) nonRedunID: the location of file folder containing txt files written
                                ligand nonredundant pdbID for each of 233 ligand
       o output: for each of 233 ligands, write fas distance to the ligand binding into  _dismat.xlsx file,
                 and there will be 233 .xlsx files
                 each .xlsx file includes:
                 fas name: all fas binding the ligand, which are named in scheme: "pdbID_AminoAcid+AminoAcidID_atom"
                 ligand atoms: all the atom names consisting of the ligand
                 distance: the distance from fas to each of the ligand atoms.
    '''
    df = pd.read_excel(LigInfo)
    for lig in df["Ligand ID"]:
        # get the number of ligand heavy atoms
        length = df[df["Ligand ID"] == lig]["Heavy Atom Count"].values[0]
        Dic = {}
        pdbList = ChangePDBCodeToList(os.path.join(NonRedunID, lig))
        for pdbcode_chain in pdbList:
            pdbcode = pdbcode_chain[:4] + ".pdb"
            if pdbcode not in os.listdir(PDBDir):
                print("Reference PRO structure %s is not available,Download and save it to PDBDir" % pdbcode)
                cmd = "curl https://files.rcsb.org/download/" + pdbcode + " -o " + PDBDir + "/" + pdbcode
                os.system(cmd)
            Dic_temp = GetAtomToLigandDistanceDict(pdbcode, lig, length, pdbcode_chain[-1])
            Dic.update(Dic_temp)
        data = pd.DataFrame(Dic)
        data.T.to_excel(os.path.join(DistFolder, lig + '_dismat.xlsx'))



