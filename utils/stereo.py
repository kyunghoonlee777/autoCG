import os
import copy
import numpy as np
import random
import sys
import itertools

from autoCG import chem
from autoCG.utils import frag
from autoCG.utils.ic import *
from autoCG.utils import process

import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Geometry import Point3D
from scipy.spatial.transform import Rotation as R
from rdkit.Chem import rdmolops

def find_potential_chirals_and_stereo_bonds(ace_mol):
    rd_mol = ace_mol.get_rd_mol()
    AllChem.EmbedMolecule(rd_mol)
    chiral_center_infos = Chem.FindMolChiralCenters(rd_mol,includeUnassigned=True)
    chiral_indices = [info[0] for info in chiral_center_infos]
    chiral_bonds = []
    Chem.FindPotentialStereoBonds(rd_mol)
    for bond in rd_mol.GetBonds():
        if bond.GetStereo() == Chem.BondStereo.STEREOANY:
            start = bond.GetBeginAtomIdx()
            end = bond.GetEndAtomIdx()
            if start>end:
                chiral_bonds.append((end,start))
            else:
                chiral_bonds.append((start,end))
    return chiral_indices,chiral_bonds

def find_complex_center_indices(ace_mol):
    adj_matrix = ace_mol.get_adj_matrix()
    valency_list = np.sum(adj_matrix,axis=1)
    indices = np.where(valency_list>4)[0].tolist()
    return indices


# Only works for molecules that have the same atom ordering ...
def get_atom_stereo(mol,atom_index):
    adj_matrix = mol.get_adj_matrix()
    neighbor_indices = np.where(adj_matrix[atom_index]>0)[0].tolist()
    neighbor_indices.sort()
    if len(neighbor_indices) != 4:
        return None
    coordinate_list = mol.get_coordinate_list()
    # Make mirror plane
    a,b,c,d = neighbor_indices
    v1 = get_vector(coordinate_list,atom_index,a)
    v2 = get_vector(coordinate_list,atom_index,b)
    n = get_cross_vector(v1,v2,True)
    # Compare sign with the third vector ...
    v3 = get_vector(coordinate_list,atom_index,c)
    v3 /= np.linalg.norm(v3)
    c = np.dot(n,v3)
    if c > 0:
        c = 1
    else:
        c = -1
    return c


def get_bond_stereo(mol,bond): # Must be double ...
    bo_matrix = mol.get_bo_matrix()
    if bo_matrix is None:
        print ('Cannot assign E/Z !!!')
        return None
    start,end = bond
    if bo_matrix[start][end] == 1:
        print ('The given bond is not multiple !!!')
        return None
    start_neighbor_indices = np.where(bo_matrix[start]>0)[0].tolist()
    end_neighbor_indices = np.where(bo_matrix[end]>0)[0].tolist()
    start_neighbor_indices.remove(end)
    end_neighbor_indices.remove(start)
    if len(start_neighbor_indices) != 2 and len(end_neighbor_indices) != 2:
        print ('Cannot have stereo isomer ...')
        return None
    start_neighbor_indices.sort()
    end_neighbor_indices.sort()
    coordinate_list = mol.get_coordinate_list()
    # Take first element for each neighbor
    a = start_neighbor_indices[0]
    d = end_neighbor_indices[0]
    angle = get_dihedral_angle(coordinate_list,a,start,end,d) 
    return np.cos(angle)


def invert_molecule(ace_mol, center_idx, freezing_indices, moving_indices, degree, step, save = False):

    adj_matrix = ace_mol.get_adj_matrix()
    coord_list = ace_mol.get_coordinate_list()    
    try:
        rd_mol = ace_mol.get_rd_mol()
    except:
        chg_list, bo_matrix = process.get_chg_and_bo(ace_mol, ace_mol.chg)
        ace_mol.set_atom_feature(chg_list, 'chg')
        ace_mol.chg_list = chg_list
        ace_mol.bo_matrix = bo_matrix
        rd_mol = ace_mol.get_rd_mol()

    zero_point = coord_list[center_idx]

    idx1, idx2 = moving_indices

    AllChem.EmbedMolecule(rd_mol)
    for i in range(len(coord_list)):
        x, y, z = coord_list[i]-zero_point
        rd_mol.GetConformer().SetAtomPosition(i, Point3D(x, y, z))
    v1 = coord_list[idx1] - coord_list[center_idx]
    v2 = coord_list[idx2] - coord_list[center_idx]
    v1 /= np.linalg.norm(v1)
    v2 /= np.linalg.norm(v2)
    rotation_axis = v1 + v2

    #print (target_atoms)        
    rotation_axis = rotation_axis / np.linalg.norm(rotation_axis)
    rotation = R.from_rotvec(rotation_axis * degree / 180 * np.pi)

    if save:
        f = open('rotate.xyz','w')
    n = len(ace_mol.atom_list)
    freezing_indices = freezing_indices + [center_idx]
    for i in range(step):
        conformer = rd_mol.GetConformer()
        #zero_point = coord_list[center[0]]
        coord_list = conformer.GetPositions()
        #coord_list = coord_list - zero_point
        target_atom_coords = coord_list[moving_indices]
        target_atom_coords = rotation.apply(target_atom_coords)

        for j in range(len(moving_indices)):
            x, y, z = target_atom_coords[j]
            conformer.SetAtomPosition(moving_indices[j], Point3D(x, y, z))

        ff = AllChem.UFFGetMoleculeForceField(rd_mol)
        ff.Initialize()

        for j in freezing_indices + moving_indices:
            ff.AddFixedPoint(j)

        if len(freezing_indices) + len(moving_indices) < len(coord_list):
            ff.Minimize()
        if save:
            f.write(f'{n}\n\n')
            content = ''
            for i in range(n):
                x, y, z = conformer.GetAtomPosition(i)
                element = rd_mol.GetAtomWithIdx(i).GetSymbol()
                content = content + f'{element:2} {x:12.8f} {y:12.8f} {z:12.8f}\n'
            f.write(content)
            #print (n)
            #print ()
            #print (content)
                
    if save:
        f.write('\n')
        f.close()

    ff = AllChem.UFFGetMoleculeForceField(rd_mol)
    ff.AddFixedPoint(center_idx)
    ff.Initialize()
    ff.Minimize()
    energy = ff.CalcEnergy()

    coord_list = rd_mol.GetConformer().GetPositions()
    for i in range(len(coord_list)):
        x, y, z = coord_list[i] + zero_point
        ace_mol.atom_list[i].x = x
        ace_mol.atom_list[i].y = y
        ace_mol.atom_list[i].z = z
    ace_mol.energy = energy/627.5095

def change_RS(ace_mol, center_idx,degree = 10,save=False):

    if degree == 0:
        degree = 10
    step = int(180/degree)

    # Make combinations
    reference_coordinate_list = ace_mol.get_coordinate_list()
    adj_matrix = ace_mol.get_adj_matrix()
    neighbor_indices = np.where(adj_matrix[center_idx]>0)[0].tolist()
    neighbor_indices.sort()

    previous_atom_stereo = get_atom_stereo(ace_mol,center_idx)
    converted = False

    n = len(neighbor_indices)
    nums = list(range(n))
    pairs = itertools.combinations(nums,2)

    nums = set(range(n))
    combinations = []

    for pair in pairs:
        combination = [list(nums-set(pair)),pair]
        combinations.append(combination)
    
    for combination in combinations:
        freezing_indices = combination[0]
        moving_indices = combination[1]
        freezing_indices = [neighbor_indices[i] for i in freezing_indices]
        moving_indices = [neighbor_indices[i] for i in moving_indices]
        invert_molecule(ace_mol, center_idx,freezing_indices,moving_indices, degree, step,True)
        new_atom_stereo = get_atom_stereo(ace_mol,center_idx)
        if new_atom_stereo != previous_atom_stereo:
            converted = True
            break
        else:
            process.locate_molecule(ace_mol,reference_coordinate_list)
            degree = -degree
            invert_molecule(ace_mol, center_idx,freezing_indices,moving_indices, degree, step,save)
            new_atom_stereo = get_atom_stereo(ace_mol,center_idx)
            if new_atom_stereo != previous_atom_stereo:
                converted = True
                break
    if converted:
        print ('Well converted !!!')
    else:
        print ('Not well converted ...', center_idx)    

def change_EZ(ace_mol, idx_1, idx_2,degree = 10,save=False):

    if degree == 0:
        degree = 10
    step = int(180/degree)

    # Make combination
    adj_matrix = ace_mol.get_adj_matrix()
    neighbor_indices_1 = np.where(adj_matrix[idx_1]>0)[0].tolist()
    neighbor_indices_1.sort()

    neighbor_indices_2 = np.where(adj_matrix[idx_2]>0)[0].tolist()
    neighbor_indices_2.sort()

    freezing_indices = neighbor_indices_1 + [idx_1,idx_2]
    moving_indices = neighbor_indices_2
    moving_indices.remove(idx_1)
    center_idx = idx_1

    previous_bond_stereo = get_bond_stereo(ace_mol,[idx_1,idx_2])
    invert_molecule(ace_mol, center_idx,freezing_indices,moving_indices, degree, step,save)
    new_bond_stereo = get_bond_stereo(ace_mol,[idx_1,idx_2])
    #print (new_bond_stereo, previous_bond_stereo)
    
    if new_bond_stereo == previous_bond_stereo:
        process.locate_molecule(ace_mol,reference_coordinate_list)
        degree = -degree
        invert_molecule(ace_mol, center_idx,freezing_indices,moving_indices, degree, step,save)
        new_bond_stereo = get_bond_stereo(ace_mol,[idx_1,idx_2])
        #print (new_bond_stereo, previous_bond_stereo)


def get_valid_atom_stereo_infos(ace_mol,atom_indices):
    atom_stereo_infos = dict()
    for atom_index in atom_indices:
        current_atom_stereo_info = get_atom_stereo(ace_mol,atom_index)
        if current_atom_stereo_info is None:
            continue
        atom_stereo_infos[atom_index] = current_atom_stereo_info
    return atom_stereo_infos


def get_valid_bond_stereo_infos(ace_mol,bonds):
    bond_stereo_infos = dict()
    for bond in bonds:
        current_bond_stereo_info = get_bond_stereo(ace_mol,bond)
        if current_bond_stereo_info is None:
            continue
        else:
            bond_stereo_infos[bond] = current_bond_stereo_info
    return bond_stereo_infos


def get_desired_stereoisomer(ace_mol,atom_stereo_infos, bond_stereo_infos):
    considering_atom_indices = []
    considering_bonds = []
    # Collect indices that have different stereo info
    current_atom_stereo_infos = get_valid_atom_stereo_infos(ace_mol,list(atom_stereo_infos.keys()))
    current_bond_stereo_infos = get_valid_bond_stereo_infos(ace_mol,list(bond_stereo_infos.keys()))
    for atom_index in current_atom_stereo_infos:
        current_atom_stereo_info = current_atom_stereo_infos[atom_index]
        desired_atom_stereo_info = atom_stereo_infos[atom_index]
        if current_atom_stereo_info * desired_atom_stereo_info < 0:
            considering_atom_indices.append(atom_index)
    for bond in current_bond_stereo_infos:
        current_bond_stereo_info = current_bond_stereo_infos[bond]
        desired_bond_stereo_info = bond_stereo_infos[bond]
        if current_bond_stereo_info * desired_bond_stereo_info < 0:
            considering_bonds.append(bond)
    
    invert_stereos(ace_mol,considering_atom_indices,considering_bonds)
    return ace_mol


def invert_stereos(ace_mol,atom_indices,bonds):

    #print ('before')
    #ace_mol.print_coordinate_list()
    for atom_index in atom_indices:
        change_RS(ace_mol,atom_index)
    #print ('after')
    #ace_mol.print_coordinate_list()

    for bond in bonds:
        change_EZ(ace_mol,bond[0],bond[1])
   
 
def enumerate_stereoisomers_for_complex_center(ace_mol,center_idx,degree = 10, save = False):

    if degree == 0:
        degree = 10
    step = int(180/degree)

    # Make combinations
    reference_coordinate_list = ace_mol.get_coordinate_list()
    adj_matrix = ace_mol.get_adj_matrix()
    neighbor_indices = np.where(adj_matrix[center_idx]>0)[0].tolist()
    neighbor_indices.sort()
    
    n = len(neighbor_indices)
    nums = list(range(n))
    pairs = itertools.combinations(nums,2)

    nums = set(range(n))
    combinations = []

    for pair in pairs:
        combination = [list(nums-set(pair)),pair]
        combinations.append(combination)

    molecules = [ace_mol.copy()]  
    for combination in combinations:
        freezing_indices = combination[0]
        moving_indices = combination[1]
        freezing_indices = [neighbor_indices[i] for i in freezing_indices]
        moving_indices = [neighbor_indices[i] for i in moving_indices]
        new_mol = ace_mol.copy()
        invert_molecule(new_mol, center_idx,freezing_indices,moving_indices, degree, step,True)
        molecules.append(new_mol)

    return molecules

def enumerate_stereoisomers_for_complex_centers(ace_mol,complex_indices,degree = 10):
    stereoisomers = [ace_mol.copy()]
    for center_idx in complex_indices:
        new_stereoisomers = []
        for stereoisomer in stereoisomers:
            new_stereoisomers += enumerate_stereoisomers_for_complex_center(stereoisomer,center_idx)
        stereoisomers = new_stereoisomers
    return stereoisomers  

def enumerate_stereoisomers(ace_mol,atom_indices=[],bonds=[],complex_indices=[]):
    stereo_atom_infos = get_valid_atom_stereo_infos(ace_mol, atom_indices)
    stereo_bond_infos = get_valid_bond_stereo_infos(ace_mol, bonds)
    considering_atom_indices = [idx for idx in atom_indices if idx in stereo_atom_infos]
    considering_bond_indices = [bond for bond in bonds if bond in stereo_bond_infos]
     
    #considering_atom_indices = [idx for idx in atom_indices]
    #considering_bond_indices = [bond for bond in bonds]
   
    n1= len(considering_atom_indices)
    n2= len(considering_bond_indices)
    m = 2**(n1+n2)
    bins = [bin(i)[2:].rjust(n1+n2,'0') for i in range(m)]
    bins = [(bin_num[:n1], bin_num[n1:]) for bin_num in bins]


    # Enumerate R/S and E/Z
    stereoisomers = []
    for atom_bins, bond_bins in bins:
        new_mol = ace_mol.copy()
        # atom
        changed_atom_indices = list()
        for atom_idx, atom_bin in zip(considering_atom_indices, list(atom_bins)):
            if atom_bin == '1':
                changed_atom_indices.append(atom_idx)
        # bond
        changed_bond_indices = list()
        for bond_idx, bond_bin in zip(considering_bond_indices, list(bond_bins)):
            if bond_bin == '1':
                changed_bond_indices.append(bond_idx)
        invert_stereos(new_mol, changed_atom_indices, changed_bond_indices)
        stereoisomers.append(new_mol)

    # Enumerate complex center indices
    final_stereoisomers = []
    for stereoisomer in stereoisomers:
        final_stereoisomers += enumerate_stereoisomers_for_complex_centers(stereoisomer,complex_indices)

    return final_stereoisomers


if __name__ == '__main__':
    molecule = chem.Molecule(sys.argv[1])
    reference_molecule = molecule.sample_conformers(1)[0] 
    reference_molecule.sanitize()
    atom_indices, bonds = find_potential_chirals_and_stereo_bonds(molecule) 
    complex_indices = find_complex_center_indices(molecule) 
    #print (atom_indices, bonds)
    #'''
    atom_indices, bonds = find_potential_chirals_and_stereo_bonds(molecule) 
    print ('reference')
    reference_molecule.print_coordinate_list()
    atom_stereo_infos = dict()
    bond_stereo_infos = dict()
    for atom_index in atom_indices:
        atom_stereo_infos[atom_index] = -get_atom_stereo(reference_molecule,atom_index)
    for bond in bonds:
        bond_stereo_infos[bond] = -get_bond_stereo(reference_molecule,bond)
    
    #print (atom_stereo_infos, bond_stereo_infos)
    #new_molecule = get_desired_stereoisomer(reference_molecule,atom_stereo_infos,bond_stereo_infos)
    stereoisomers = enumerate_stereoisomers(reference_molecule,atom_indices,bonds,complex_indices)

    for i,stereoisomer in enumerate(stereoisomers):
        print (f'{i+1}th stereoisomer ...')
        stereoisomer.print_coordinate_list()


