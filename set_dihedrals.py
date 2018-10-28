import itertools
import sys
import os

from rdkit import Chem
from rdkit.Chem import rdMolTransforms, rdMolAlign
import openbabel

from qmconftool import QMMol


def find_dihedral_idx(mol,smarts_patt):
    
    patt_mol = Chem.MolFromSmarts(smarts_patt)

    matches = mol.GetSubstructMatches(patt_mol)

    unique_match = list()
    match_list = list()

    for m in matches:
        if m[:3] not in match_list:
            unique_match.append(m)
            match_list.append(m[:3])

    if len(unique_match) != 2:
        print("more than two dihedrals in " + filename)
        quit()

    return unique_match


def changeAndOpt(rdkit, theta): 
    
    Chem.SanitizeMol(rdkit)
    initconf = rdkit.GetConformer()
    
    # set outer most dihedral to 180 degrees.
    smarts_patt = "C-S-C-[C,Si,Ge;H0]"
    outer_dihedral_idx = find_dihedral_idx(rdkit, smarts_patt)

    for k, i, j, l in outer_dihedral_idx:
        rdMolTransforms.SetDihedralDeg(initconf, k,i,j,l, 180.0)

    # change second outmost dihedral with +-120 degrees.
    patt = "S-C-[C,Si,Ge;H0]-[C,Si,Ge]"
    dihedral_idx = find_dihedral_idx(rdkit, patt)
    
    new_angles = list()
    for k, i, j, l in dihedral_idx:
        init_dihedral_angle = rdMolTransforms.GetDihedralDeg(initconf, k,i,j,l)
        new_angles.append([init_dihedral_angle + x*theta for x in range(int(360./theta))])
    
    angle_combinations = list(itertools.product(*new_angles)) # all combinations.
    

    for dihedrals in angle_combinations:
        
        for (k,i,j,l), angle in zip(dihedral_idx, dihedrals):
            rdMolTransforms.SetDihedralDeg(initconf, k,i,j,l, angle )


        rdkit.AddConformer(initconf, assignId=True)
    
    rdMolAlign.AlignMolConformers(rdkit)
    
    mol_list = list()
    for idx, conf in enumerate(rdkit.GetConformers()):

        if idx == 0:
            continue
        
        sdf_txt = Chem.SDWriter.GetText(rdkit, conf.GetId())
        m = Chem.MolFromMolBlock(sdf_txt, removeHs=False)
        
        conf_name =  m.GetProp("_Name") + "-" + str(idx-1)
        m.SetProp("_Name", conf_name)

        mol_list.append(m)

    # Optimize structures with new dihedrals.
    confqmmol = QMMol(mol_list, fmt="mol_list", charge=0, multi=1, charged_fragments=True)
    confqmmol.optimize(program="xtb", method="opt", cpus=24, babelAC=True)


    # Write xyz files of conformers
    for newConf in confqmmol.GetConformers():
       
        obConversion = openbabel.OBConversion()
        obConversion.SetInAndOutFormats("sdf", "xyz")

        newConfm = openbabel.OBMol()
        obConversion.ReadString(newConfm, Chem.MolToMolBlock(newConf))

        new_xyz = obConversion.WriteString(newConfm)
        
        with open(newConf.GetProp("_Name") + ".xyz", 'w') as f:
           f.write(new_xyz)
        

if __name__ == "__main__":

    mols = list()
    for fname in os.listdir('.'):
    
        if fname.endswith("sdf"):
            m = Chem.MolFromMolFile(fname, removeHs=False)
            m.SetProp("_Name", fname.split('.')[0])
            
            mols.append(m)
    
    # optimize mol with xTB.
    qmmol = QMMol(mols, fmt="mol_list", charge=0, multi=1, charged_fragments=True)
    qmmol.optimize(program="xtb", method="opt", cpus=47, babelAC=True)
    
    
    theta_change = 120.
    
    # Change dihedrals 
    for c in qmmol.GetConformers():
        changeAndOpt(c, theta_change)
