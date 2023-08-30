import os
from rdkit import Chem
from rdkit.Chem import AllChem

from rdkit.Chem.rdmolops import AssignStereochemistry

#mol = Chem.MolFromXYZFile("~/test_directory/jcc_Jul/Schlegel/SN2_test/1/R.xyz")
mol = Chem.MolFromSmiles("[H]C([H])(P(C([H])([H])[H])(F)(C([H])([H])[H])Cl)[H]")
Chem.SanitizeMol(mol)
mol = Chem.AddHs(mol)
"""
for x in Chem.FindPotentialStereo(mol):
    print(x.type, x.centeredOn)
    atom = mol.GetAtomWithIdx(x.centeredOn)
    print(atom.GetSymbol())
"""

desired_chirality = Chem.ChiralType.CHI_TRIGONALBIPYRAMIDAL
Chem.AssignStereochemistry(mol, 1, desired_chirality)

AllChem.EmbedMultipleConfs(mol, numConfs=10, params=Chem.rdDistGeom.srETKDGv3())
for conf in mol.GetConformers():
    AllChem.UFFOptimizeMolecule(mol,confId=conf.GetId())
    block = Chem.MolToXYZFile(mol, f"{conf.GetId()}.xyz", conf.GetId())
