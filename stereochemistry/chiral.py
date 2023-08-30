### python provided modules ###
import os
import sys
import time
import typing

### extra common libraries ###
import numpy as np

### ace-reaction libraries ###
from autoCG import chem
from autoCG.utils import process


def mirrorPoint(position, nVec):
    # reflect coordinate by a plane with normal vector of nVec,
    # which containing origin.
    proj = np.dot(np.array(position), np.array(nVec))
    return position - 2 * proj


def mirrorReflection(coordinate_list, nVec, center_of_mirror):
    coords = np.array(coordinate_list)
    tr_coords = coords - np.array(center_of_mirror)
    normV = np.array(nVec) / np.linalg.norm(nVec)
    projs = np.reshape(tr_coords @ normV, (-1, 1)) * normV

    return coords - 2 * projs


def changeRS(
    mol: chem.Molecule, chiral_center_index: int, adj_index1: int, adj_index2: int
):
    coordinate_list = mol.get_coordinate_list()

    com = coordinate_list[chiral_center_index]
    nVec = np.cross(
        coordinate_list[chiral_center_index] - coordinate_list[adj_index1],
        coordinate_list[chiral_center_index] - coordinate_list[adj_index2],
    )
    
    new_coords = mirrorReflection(coordinate_list, nVec, com)
    process.locate_molecule(mol, new_coords)

    return 


def changeEZ(coordinate_list, double_bond_info, adj_matrix):
    return


if __name__ == "__main__":
    mol = chem.Molecule("CC(Cl)Br")
    cl = mol.make_3d_coordinate()
    process.locate_molecule(mol, cl)
    mol.write_geometry("mirror_before.xyz")

    nVect = np.cross(np.array(cl[1] - cl[2]), np.array(cl[1] - cl[3]))
    cc = cl[1]

    new = mirrorReflection(cl, nVect, cc)
    process.locate_molecule(mol, new)
    mol.write_geometry("mirror_after.xyz")
