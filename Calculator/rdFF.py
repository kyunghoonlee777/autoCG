import os
from copy import deepcopy

import numpy as np
import rdkit.ForceField.rdForceField as FF
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdForceFieldHelpers import UFFGetMoleculeForceField

### ace-reaction libraries ###
from autoCG import chem
from autoCG.utils import process

def get_rd_mol3D(ace_mol):
    pos = ace_mol.get_coordinate_list()

    rd_mol = ace_mol.get_rd_mol()
    Chem.SanitizeMol(rd_mol)
    AllChem.EmbedMolecule(rd_mol)
    for i in range(rd_mol.GetNumAtoms()):
        rd_mol.GetConformer().SetAtomPosition(i, pos[i])

    return rd_mol

class UFFOptimizer:
    """UFF Optimizer Class (utilizing RDKit)"""

    MAXITS = 200
    ETOL = 0.1

    def __init__(self):
        self.rd_mol = None
        self.opt = None
        self.save_directory = ""  # rdkit does not need this
        self.working_directory = ""  # rdkit does not need this
        self.command = "rdkit.ForceField"  # rdkit does not need this

    def setUFFoptimizer(self, rd_mol):
        opt = UFFGetMoleculeForceField(rd_mol)
        opt.Initialize()
        self.opt = opt

    def rdoptimize(self, rd_mol, max_cycles=MAXITS, e_tol=ETOL):
        assert self.opt is not None

        for i in range(max_cycles):
            energy = self.opt.CalcEnergy()
            self.opt.Minimize(maxIts=1)
            if abs(energy - self.opt.CalcEnergy()) < e_tol:
                break

        return rd_mol.GetConformer().GetPositions()

    ### Common functions

    def get_energy(self, molecule):
        #assert self.opt is not None
        opt = UFFGetMoleculeForceField(get_rd_mol3D(molecule))
        opt.Initialize()

        return opt.CalcEnergy()

    def optimize_geometry(
        self,
        molecule,
        max_cycles: int = MAXITS,
        e_tol: float = ETOL,
    ):
        rd_mol = get_rd_mol3D(molecule)

        self.setUFFoptimizer(rd_mol)
        new_coords = self.rdoptimize(rd_mol, max_cycles, e_tol)

        process.locate_molecule(molecule, np.array(new_coords))

        return

    def relax_geometry(
        self,
        molecule,
        contraints=None,
        chg=None,
        multiplicity=None,
        file_name=None,
        num_relaxation=MAXITS,
        maximal_displacement=1000,
        save_directory=None,
    ):
        return self.optimize_geometry(
            molecule, bo_matrix=None, max_cycles=max_cycles, e_tol=e_tol
        )

    def change_working_directory(self, path):
        pass

    def clean_scratch(self):
        pass


if __name__ == "__main__":
    """
    mol = Chem.MolFromSmiles("CCCC")
    mol = Chem.AddHs(mol)
    Chem.SanitizeMol(mol)
    """

    mol = Chem.MolFromXYZFile("initial_ts.xyz")
    tsPos = mol.GetConformer().GetPositions()

    # mol = Chem.MolFromXYZFile("P.xyz")
    ace_mol = chem.Molecule("P.xyz")
    bo = process.get_bo_matrix_from_adj_matrix(ace_mol, chg=-1)
    ace_mol.bo_matrix = bo

    mol = ace_mol.get_rd_mol()
    Chem.SanitizeMol(mol)
    AllChem.EmbedMolecule(mol)

    for i in range(mol.GetNumAtoms()):
        mol.GetConformer(0).SetAtomPosition(i, tsPos[i])
    # Chem.SanitizeMol(mol)
    print(len(mol.GetBonds()))
    for bond in mol.GetBonds():
        print(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())

    # AllChem.EmbedMultipleConfs(mol, 1)

    for c in mol.GetConformers():
        before = c.GetPositions()

    Chem.MolToXYZFile(mol, "before.xyz")

    UFF = UFFGetMoleculeForceField(mol)
    UFF.Initialize()

    MAXITS = 1000
    ETOL = 10.0
    # ETOL = 0.1
    # ETOL=1e-6
    XYZ = ""

    energy = UFF.CalcEnergy()
    print(f"Energy before opts(kcal/mol):", energy)
    XYZ += Chem.MolToXYZBlock(mol)

    for i in range(1, MAXITS+1):
        UFF.Minimize(maxIts=1)
        print(f"Energy after {i} opts(kcal/mol):", UFF.CalcEnergy())
        XYZ += Chem.MolToXYZBlock(mol)
        if abs(energy - UFF.CalcEnergy()) < ETOL:
            print(f"***Convergence condition satisfied. (dE < {ETOL} kcal/mol)")
            break
        energy = UFF.CalcEnergy()
        # print(type((mol.GetConformer(-1).GetPositions())))

    with open("UFF.xyz", "w") as f:
        f.write(XYZ)

    print(dir(UFF))
    print(UFF.mol)
