import copy

from itertools import islice


from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolTransforms

def change_molecule_recursive(smiles_list,substitutions,rxn_smarts_list):
    raw_smiles = copy.copy(smiles_list)
    for i in range(substitutions):
        for smiles in islice(raw_smiles,i,len(raw_smiles)):
            mol = Chem.MolFromSmiles(smiles)

            for rxn_smarts in rxn_smarts_list:
                try:
                    rxn = AllChem.ReactionFromSmarts(rxn_smarts)
                    new_mols = rxn.RunReactants((mol,))
                except:
                    continue

                for new_mol in new_mols:
                    new_smiles = Chem.MolToSmiles(new_mol[0])
                    if new_smiles not in raw_smiles:
                        raw_smiles.append(new_smiles)

    for smiles in raw_smiles:
        try:
            mol = Chem.MolFromSmiles(smiles)
        except:
            continue
        if mol != None:
            smiles = Chem.MolToSmiles(Chem.MolFromSmiles(smiles))
            if smiles not in smiles_list:
                smiles_list.append(smiles)
        
    return smiles_list


def get_molecules(smiles_list):
    molecules = []
    for smiles in smiles_list:
        mol = Chem.MolFromSmiles(smiles)
        molecules.append(mol)
    
    return molecules


def change_molecule_saturate(smiles_list,smarts):
    new_smiles_list = []
    for smiles in smiles_list:
        mol = Chem.MolFromSmiles(smiles)
        patt = Chem.MolFromSmarts(smarts.split(">>")[0])
        if mol.HasSubstructMatch(patt):
            while mol.HasSubstructMatch(patt):
                rxn = AllChem.ReactionFromSmarts(smarts)
                ps = rxn.RunReactants((mol,))
                mol = ps[0][0]
                mol = Chem.MolFromSmiles(Chem.MolToSmiles(mol))
            new_smiles_list.append(Chem.MolToSmiles(mol))
        else:
            new_smiles_list.append(smiles)
    
    return new_smiles_list


def write_xyz_file(m, name):
    
    number_of_atoms = m.GetNumAtoms()
    atom_symbols = [a.GetSymbol() for a in m.GetAtoms()]

    for i, conf in enumerate(m.GetConformers()):
        
        file_name = name+"+"+str(i)+".xyz"

        xyz_string = str(number_of_atoms) + '\n'
        xyz_string += file_name + '\n'
        

        for symbol, pos in zip(atom_symbols, conf.GetPositions()):
           xyz_string +=  "{}  {:10.5f} {:10.5f} {:10.5f}\n".format(symbol, *tuple(pos))
    
        with open(file_name, 'w') as f:
            f.write(xyz_string)
        

def write_files(molecules):
    patt = Chem.MolFromSmarts('[C]([H])([H])([H])S[CH2][*;!#1]')
    
    for i,mol in enumerate(molecules):
        name = 'comp'+str(i)
        mol = Chem.AddHs(mol)
        title = mol.GetSubstructMatches(patt)
        
        number_of_atoms = mol.GetNumAtoms()
        first7atoms = list(title[0])
        last7atoms = list(title[1])
        missing = set(range(number_of_atoms)).difference(first7atoms+last7atoms)
        missing = list(missing)
        order = first7atoms + missing + last7atoms
        mol = Chem.RenumberAtoms(mol, order)

        AllChem.EmbedMultipleConfs(mol,numConfs=50)
        AllChem.UFFOptimizeMoleculeConfs(mol,maxIters=1000)
     
        write_xyz_file(mol,name)


if __name__ == '__main__':
    
    smiles_list = ['NC1(CC2)CCC2(N)CC1']
        
    change_elements = ['[C:1]>>[Si:1]','[C:1]>>[Ge:1]'] # change C to Si or Ge
    
    SiH_to_SiMe = '[Si;H2,H1:1]>>[Si:1]C' # methylate Si
    GeH_to_GeMe = '[Ge;H2,H1:1]>>[Ge:1]C' # methylate Ge

    N2C = ['[N;!$(N[CH2][CH2]):1][C:2]>>[N:1]C[C:2]']
    N2term = '[N:1]>>[C:1]SC'
        
    substitutions = 8
    
    smiles_list = change_molecule_recursive(smiles_list,substitutions,change_elements)    
    smiles_list = change_molecule_saturate(smiles_list,SiH_to_SiMe)
    smiles_list = change_molecule_saturate(smiles_list,GeH_to_GeMe)
    smiles_list = change_molecule_saturate(smiles_list,N2term)
    
    
    molecules = get_molecules(smiles_list)
    write_files(molecules)
    
    
