# Conversion of SMILES to XYZ format using RDKit
import os
import time
from rdkit import Chem
from rdkit.Chem import AllChem

def get_xyz(mol):
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    AllChem.EmbedMultipleConfs(mol, numConfs=2)
    AllChem.MMFFOptimizeMolecule(mol)
    
    n_atoms = mol.GetNumAtoms()
    conf = mol.GetConformer()
    xyz = ''
    for i in range(mol.GetNumAtoms()):
        pos = conf.GetAtomPosition(i)
        xyz += '{} {} {} {}\n'.format(mol.GetAtomWithIdx(i).GetSymbol(), pos.x, pos.y, pos.z)
    return xyz, n_atoms

def smiles_to_xyz(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    xyz, n_atoms = get_xyz(mol)
    return xyz, n_atoms

# Save a xyz file in the current directory with the name of the molecule
def create_xyz_from_smiles(smiles, xyzname):
    xyz, n_atoms = smiles_to_xyz(smiles)
    if xyz is None:
        return None
    # If the file already exists, overwrite it and save the previous one with a different name
    if os.path.exists(smiles + '.xyz'):
        os.rename(smiles + '.xyz', smiles + '_{}.xyz'.format(str(time.time())))

    if xyzname is None:
        smiles = smiles.replace('/', '_')
        smiles = smiles.replace('\\', '__')
        xyzname = smiles
    with open(xyzname + '.xyz', 'w') as f:
        f.write(f'{n_atoms}\n\n')
        f.write(xyz)
