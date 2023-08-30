### python provided modules ###
import itertools
from collections.abc import Iterable
from copy import deepcopy

### extra common libraries ###
import numpy as np
from rdkit.Chem import FindMolChiralCenters, FindPotentialStereoBonds
from scipy.spatial.transform import Rotation

### ace-reaction libraries ###
from autoCG import chem
from autoCG.utils import conformation, process


def powerset(l):
    return itertools.chain.from_iterable(
        itertools.combinations(l, s) for s in range(len(l) + 1)
    )


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
    mol: chem.Molecule,
    chiral_center_index: int,
    adj_index1: int,
    adj_index2: int,
):
    coords = np.array(mol.get_coordinate_list())

    com = coords[chiral_center_index]
    nVec = np.cross(
        coords[chiral_center_index] - coords[adj_index1],
        coords[chiral_center_index] - coords[adj_index2],
    )
    new_coords = mirrorReflection(coords, nVec, com)

    # select atoms to be reflected
    adj = mol.get_adj_matrix().copy()
    adj[chiral_center_index][adj_index1] -= 1
    adj[adj_index1][chiral_center_index] -= 1
    adj[chiral_center_index][adj_index2] -= 1
    adj[adj_index2][chiral_center_index] -= 1
    target = list(process.get_molecule_group(adj, chiral_center_index))
    coords[target] = new_coords[target]

    process.locate_molecule(mol, coords)

    return


def sampleRS(
    mol: chem.Molecule, scope: Iterable[int] = None, check_chirality=False
):
    if scope is None:
        scope = range(len(mol.get_z_list()))

    atom_ids = np.arange(len(mol.get_z_list()))
    # print("atom_ids", atom_ids)
    lookups = np.array(scope)
    # print("lookups", lookups)
    # print("get_valency_list", np.flatnonzero(mol.get_valency_list() == 4))
    ccs = lookups[np.isin(lookups, np.flatnonzero(mol.get_valency_list() == 4))]
    if check_chirality:
        rd_mol = mol.get_rd_mol()
        new_ccs = [
            x
            for x, chi_type in FindMolChiralCenters(
                rd_mol, includeCIP=False, useLegacyImplementation=False
            )
            if x in ccs
        ]
        ccs = np.array(new_ccs)
    # print("ccs", ccs)
    parameters = []
    for center in ccs:
        adj = np.copy(mol.get_adj_matrix())
        neighbors = np.flatnonzero(adj[center])
        for neighbor in neighbors:
            adj[center][neighbor] -= 1
            adj[neighbor][center] -= 1
        frags = sorted(
            process.group_molecules(adj), key=lambda x: len(x), reverse=True
        )
        # print("center: ", center)
        # print("neighbors: ", neighbors)
        # print("frags: ", frags)
        for neighbor in neighbors:
            adj[center][neighbor] += 1
            adj[neighbor][center] += 1

        adj1 = None
        adj2 = None
        for frag in frags:
            for neighbor in neighbors:
                if adj1 is None and neighbor in frag:
                    adj1 = neighbor
                elif adj2 is None and neighbor in frag:
                    adj2 = neighbor

                if adj1 is not None and adj2 is not None:
                    break

        parameters.append((center, adj1, adj2))

    conformers = []
    # count = 0
    # print("parameters", parameters)
    for combi in powerset(parameters):
        # count += 1
        # print("count", count)
        conformer = deepcopy(mol)
        for param in combi:
            changeRS(conformer, *param)
        conformers.append(conformer)

    return conformers


def changeEndoExo(
    mol: chem.Molecule, four_point_info: tuple[int, int, int, int]
):
    """four point info should be given in the right order, representing sequential two bonds"""
    coords = np.array(mol.get_coordinate_list())
    com = np.average(coords[np.array(four_point_info)], axis=0)

    # define mirror plane by LEAST SQUARE fitting
    mat = np.hstack((coords, np.ones((coords.shape[0], 1))))
    nVect, _, _, _ = np.linalg.lstsq(
        mat[np.ix_(list(four_point_info), [0, 1, 3])],
        mat[list(four_point_info), 2],
        rcond=None,
    )
    new_coords = mirrorReflection(coords, nVect, com)

    # update only either side
    adj = np.copy(mol.get_adj_matrix())
    target = list(process.get_molecule_group(adj, four_point_info[0]))
    coords[target] = new_coords[target]

    process.locate_molecule(mol, coords)

    return


def changeEZ(mol: chem.Molecule, double_bond_info: tuple[int, int]):
    """The bond should not be the part of a ring"""
    conformation.rotate_bond(
        mol,
        double_bond_info[0],
        double_bond_info[1],
        180,
        allow_double_bond=True,
    )

    return


def sampleEZ(
    mol: chem.Molecule, scope: Iterable[tuple[int, int]] = None, check_EZ=False
):
    if scope is None:
        dbs = np.asarray(
            [
                (a, b)
                for a, b in np.stack(
                    np.asarray(mol.get_bo_matrix() == 2).nonzero(), axis=-1
                )
                if a < b
            ]
        )
    else:
        dbs = scope

    if check_EZ:
        new_dbs = []
        rd_mol = mol.get_rd_mol()
        FindPotentialStereoBonds(rd_mol)
        for bond in rd_mol.GetBonds():
            start = bond.GetBeginAtomIdx()
            end = bond.GetEndAtomIdx()
            if end < start:
                start, end = end, start
            if (
                start,
                end,
            ) in dbs and bond.GetStereo() == Chem.BondStereo.STEREOANY:
                new_dbs.append((start, end))
        dbs = new_dbs

    conformers = []
    for combi in powerset(dbs):
        conformer = deepcopy(mol)
        for doublebond in combi:
            changeEZ(conformer, doublebond)
        conformers.append(conformer)

    return conformers


if __name__ == "__main__":
    import os

    from rdkit import Chem

    mol = chem.Molecule("/home/leejinwon/autoCG/stereochemistry/initial_ts.xyz")
    # mol = chem.Molecule("/home/leejinwon/autoCG/stereochemistry/DACP2_ts.xyz")
    # mol = chem.Molecule("C/C=C\C=CC=CC")
    # process.locate_molecule(mol, mol.make_3d_coordinate())
    # mol = Chem.MolFromXYZFile("/home/leejinwon/autoCG/stereochemistry/initial_ts.xyz")
    # mol.write_geometry("/home/leejinwon/autoCG/stereochemistry/tetra_before.xyz")
    # mol.write_geometry("/home/leejinwon/autoCG/stereochemistry/_before.xyz")
    # mol.write_geometry("/home/leejinwon/autoCG/stereochemistry/db_before.xyz")
    changeRS(mol, 1, 2, 7)
    # changeEZ(mol,(1, 2))
    # conformers = sampleEZ(mol, scope=None, check_EZ=False)
    # for conformer in conformers:
    #    conformer.write_geometry("/home/leejinwon/autoCG/stereochemistry/db_sample.xyz")

    mol.write_geometry("/home/leejinwon/autoCG/stereochemistry/tetra_after.xyz")
    # mol.write_geometry("/home/leejinwon/autoCG/stereochemistry/endo_after.xyz")
    # mol.write_geometry("/home/leejinwon/autoCG/stereochemistry/db_after.xyz")
