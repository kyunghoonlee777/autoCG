
"""
---generate_molecule.py---
Generate 2D graph molecule in ACE-Reaction format from xyz files or SMILES
"""

import os 
import subprocess
import copy
import numpy as np 
import random

from rdkit import Chem

from autoCG import chem

from autoCG.utils import frag

def locate_molecule(ace_molecule,coordinate_list,update = False):
    """ Locates atoms according to coordinate_list, be cautious on ordering of atoms
    Args:
        |  ace_molecule (<class 'Molecule'>): class instance of Molecule
        |  coordinate_list (list of list (x,y,z)): list of 3d coordinate, where each x,y,z are float type
    Returns:
        |  No return, it direcly modifies the geometry of a given molecule
    """
    atom_list = ace_molecule.atom_list
    for i in range(len(atom_list)):
        atom = atom_list[i]
        atom.set_coordinate(coordinate_list[i])
    if update:
        ace_molecule.set_adj_matrix(None)

def translate_molecule(ace_molecule,vector):
    atom_list = ace_molecule.atom_list
    for atom in atom_list:
        translate_atom(atom,vector)


def locate_atom(atom,coordinate):
    """ Locates a single atom to input 'coordinate'
    Args:
        |  atom (<class 'Atom'>): class instance of Atom
        |  coordinate (list of float): coordinate in form of [x,y,z]
    Returns:
        |  No return, it directly locates the atom to the given coordinate
    """
    atom.x = coordinate[0]
    atom.y = coordinate[1]
    atom.z = coordinate[2]

def translate_atom(atom,vector):
    atom.x += vector[0]
    atom.y += vector[1]
    atom.z += vector[2]


def get_ace_mol_from_rd_mol(rd_molecule,add_hydrogen = True,include_stereo = False):
    """ It converts rd_molecule type info ace_molecule type
    Args:
        |  rd_molecule (<class 'rdkit.Molecule>')
    Returns:
        |  ace_molecule(<class Molecule>)
    """
    # Kekulize molecule
    rd_molecule_copy = copy.deepcopy(rd_molecule)
    try:
        if add_hydrogen:
            rd_molecule = Chem.AddHs(rd_molecule) 
        Chem.rdmolops.Kekulize(rd_molecule)
    except:
        rd_molecule = rd_molecule_copy
    bond_types = {Chem.BondType.SINGLE:1, Chem.BondType.DOUBLE:2, Chem.BondType.TRIPLE:3} 
    n = rd_molecule.GetNumAtoms()
    atom_list = []
    chg_list = []
    atom_feature = dict()
    # Make atom_list
    for i in range(n):
        rd_atom = rd_molecule.GetAtomWithIdx(i)
        ace_atom = chem.Atom()
        chg_list.append(rd_atom.GetFormalCharge())
        '''
        position = rd_molecule.GetAtomPosition(i)
        if position!=None:
            ace_atom.x = position[0]
            ace_atom.y = position[1]
            ace_atom.z = position[2]
        '''
        ace_atom.atomic_number = rd_atom.GetAtomicNum()
        atom_list.append(ace_atom)
    atom_feature['chg'] = np.array(chg_list)
    # Make bond order matrix
    bonds = rd_molecule.GetBonds()
    bo_matrix = np.zeros((n,n))
    for bond in bonds:
        begin = bond.GetBeginAtomIdx()
        end = bond.GetEndAtomIdx()
        bond_order = bond_types[bond.GetBondType()]
        bo_matrix[begin][end] = bo_matrix[end][begin] = bond_order
    ace_molecule = chem.Molecule()
    ace_molecule.atom_list = atom_list
    ace_molecule.bo_matrix = bo_matrix
    ace_molecule.atom_feature = atom_feature
    return ace_molecule

def get_permuted_molecule(molecule,permutation):
    atom_list = molecule.atom_list
    bo_matrix = molecule.bo_matrix
    adj_matrix = molecule.adj_matrix
    atom_feature = molecule.atom_feature
    new_molecule = chem.Molecule()
    if bo_matrix is not None:
        bo_matrix = get_permuted_matrix(bo_matrix,permutation)
        adj_matrix = np.where(bo_matrix>0,1,0)
    elif adj_matrix is not None:
        adj_matrix = get_permuted_matrix(adj_matrix,permutation)
    new_molecule.adj_matrix = adj_matrix
    new_molecule.bo_matrix = bo_matrix
    new_molecule.atom_list = get_permuted_atom_list(atom_list,permutation)
    new_molecule.atom_feature = get_permuted_atom_feature(atom_feature,permutation)
    return new_molecule

def get_permuted_atom_list(atom_list,permutation):
    n = len(permutation)
    new_atom_list = [None] * n
    for i in range(n):
        new_atom_list[permutation[i]] = atom_list[i]
    return new_atom_list

def get_permuted_atom_feature(atom_feature,permutation):
    new_atom_feature = dict()
    n = len(permutation)
    if atom_feature == None:
        return atom_feature
    for feature in atom_feature:
        feature_value = atom_feature[feature]
        new_feature_value = copy.deepcopy(feature_value)
        for i in range(n):
            value = permutation[i]
            new_feature_value[value] = feature_value[i]
        new_atom_feature[feature] = new_feature_value
    return new_atom_feature
     
def get_permuted_matrix(matrix,permutation):
    n = len(matrix)
    permuted_matrix = np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            permuted_matrix[permutation[i]][permutation[j]] = matrix[i][j]
    return permuted_matrix

def get_adj_matrix_from_distance(molecule,coeff = 1.35):
    """
    Returns adj_matrix from 3d coordinate of given molecule
    It recognizes bond between two atoms, if the sum of radius * coeff is less than distance between two atoms

    :param coeff(float):
        criteria for recognizing bond. If criteria gets higher, more and more bonds are generated between atoms, 
        since criteria distance for bond distance gets higher.
        Appropriate criteria value is between 0.8 ~ 1.3, here we set default value as 1.10
    
    :return adj(pyclass 'numpy.ndarray'):
        connectivity matrix between atoms
         
    """
    ratio_matrix = molecule.get_ratio_matrix()
    adj = np.where(ratio_matrix<coeff,1,0)
    np.fill_diagonal(adj,0)
    # Here, consider additional condition
    max_valency_list = molecule.get_max_valency_list()
    valency_list = np.sum(adj,axis=1)
    over_octet_indices = np.where(valency_list>max_valency_list)[0].tolist()
    # Remove bond with maximum ratio for overoctet indices that is bonded ...
    while len(over_octet_indices) > 0:
        index = over_octet_indices[0]
        removing_index = np.argmax(ratio_matrix[index])
        adj[index][removing_index] -= 1 # Remove bond
        adj[removing_index][index] -= 1
        ratio_matrix[index][removing_index] = 0 # To label that bond is removed 
        ratio_matrix[removing_index][index] = 0
        valency_list = np.sum(adj,axis = 1)
        over_octet_indices = np.where(valency_list>max_valency_list)[0].tolist() 
    return adj 

def get_bo_matrix_from_adj_matrix(molecule,chg=None,method='SumofFragments',obtain_all_resonance = False):
    """
    Returns bo_matrix from adj_matrix stored in pyclass 'Molecule'
    
    :param chg(int):
        total charge of the molecule

    :param method(str):
        Mainly 'SumofFragments' and 'Ionic' is possible
        'SumofFragments' use user defined fragment charge to evaluate bo_matrix 
        (also uses some chemical heuristics)
        'Ionic' method uses chemical heuristics to evaluate bo_matrix

    :param obtain_all_resonance(boolean):
        If True, it returns multiple possible resonance structures, therefore return as list of numpy.ndarray
        (Normally this function is not used)

    :return bo_matrix(pyclass 'numpy.ndarray' or list of pyclass 'numpy.ndarray'):
        possible bo_matrices for obtain_all_resonance=True, otherwise just single bo_matrix
    """
    atom_list = molecule.atom_list
    adj_matrix = molecule.get_adj_matrix()
    if chg is None:
        chg = molecule.get_chg()
        if chg is None:
            print ('Total charge is not specified! Provide charge information!')
            return None
    if obtain_all_resonance:
        bo_candidates = frag.AdjtoBO(atom_list,adj_matrix,chg,method,True)
        return bo_candidates
    else:
        bo_matrix = frag.AdjtoBO(atom_list,adj_matrix,chg,method,False)
        return bo_matrix

def get_chg_list_from_bo_matrix(molecule,chg,bo_matrix,method = 'SumofFragments'):
    """
    Returns chg_list from a given bo_matrix stored in pyclass 'Molecule'
    
    :param chg(int):
        total charge of the molecule

    :param bo_matrix(pyclass 'numpy.ndarray'):
        Possible bond order matrix of given molecule

    :param method(str):
        Mainly 'SumofFragments' and 'Ionic' is possible
        'SumofFragments' use user defined fragment charge to evaluate bo_matrix 
        (also uses some chemical heuristics)
        'Ionic' method uses chemical heuristics to evaluate bo_matrix

    :return chg_list(pyclass 'numpy.ndarray' (1D-array)):
        formal charge of each atom
    """
    atom_list = molecule.atom_list
    #bo_matrix_before = np.copy(bo_matrix)
    chg_list = frag.getFC(atom_list,bo_matrix,chg,method)
    return np.array(chg_list)

def get_chg_list_and_bo_matrix_from_adj_matrix(molecule,chg=None,method='SumofFragments'):
    atom_list = molecule.atom_list
    adj_matrix = molecule.get_adj_matrix()
    if chg is None:
        chg = molecule.get_chg()
        if chg is None:
            print ('Total charge is not specified! Provide charge information!')
            return None,None
    bo_matrix = get_bo_matrix_from_adj_matrix(molecule,chg,method)
    chg_list = get_chg_list_from_bo_matrix(molecule,chg,bo_matrix,method)
    return chg_list,bo_matrix

def get_reduced_intermediate(intermediate,reduce_function):
    reduced_intermediate = chem.Intermediate()
    adj_matrix = intermediate.get_adj_matrix()
    bo_matrix = intermediate.get_bo_matrix()
    chg_list = intermediate.get_chg_list()
    n = len(reduce_function)
    if type(reduce_function) is list:
        reduced_atom_list = []
        reduced_bo_matrix = None
        index_function = np.ix_(reduce_function,reduce_function)
        reduced_chg_list = []
        reduced_atom_list = [intermediate.atom_list[index] for index in reduce_function]
        reduced_adj_matrix = adj_matrix[index_function]
        if bo_matrix is not None:
            reduced_bo_matrix = bo_matrix[index_function]
        if chg_list is not None:
            reduce_chg_list = np.array([chg_list[index] for index in reduce_function]) 
    else:
        reduced_atom_list = [None] * n
        reduced_adj_matrix = np.zeros((n,n))
        reduced_chg_list = []
        reduced_bo_matrix = None
        if bo_matrix is not None:
            reduced_bo_matrix = np.zeros((n,n))
        if chg_list is not None:
            reduced_chg_list = np.zeros((n))
        for original_i in reduce_function:
            i = reduce_function[original_i]
            reduced_atom_list[i] = intermediate.atom_list[original_i]
            if chg_list is not None:
                reduced_chg_list[i] = chg_list[original_i]
            for original_j in reduce_function:
                j = reduce_function[original_j]
                reduced_adj_matrix[i][j] = adj_matrix[original_i][original_j]
                if reduced_bo_matrix is not None:
                    reduced_bo_matrix[i][j] = bo_matrix[original_i][original_j]
    reduced_intermediate.atom_list = reduced_atom_list
    reduced_intermediate.adj_matrix = reduced_adj_matrix
    reduced_intermediate.bo_matrix = reduced_bo_matrix
    if len(reduced_chg_list) > 0:
        reduced_intermediate.atom_feature['chg'] = reduced_chg_list
    return reduced_intermediate
 
def get_molecule_group(adj_matrix,index=0):
    current_list = set([index])
    total_list = set([index])
    while len(current_list) > 0:
        new_current_list = set([])
        for i in current_list:
            neighbor_list = np.where(adj_matrix[i]>0)[0].tolist()
            new_current_list = new_current_list | set(neighbor_list)
        current_list = new_current_list - total_list
        total_list = total_list | new_current_list
    return total_list

def group_molecules(adj_matrix):
    n = len(adj_matrix)
    all_indices = set(range(n))
    groups = []
    index = 0
    while len(all_indices)>0:
        indices = get_molecule_group(adj_matrix,index)
        all_indices = all_indices - indices
        groups.append(list(indices))
        if len(all_indices)>0:
            index = min(all_indices)
        else:
            break
    return groups


def is_same_connectivity(original_molecule,new_molecule,max_coeff=1.3,min_coeff=0.95,space=0.05):
    coeff = min_coeff
    is_same = False
    #new_molecule.print_coordinate_list()
    while coeff < max_coeff:
        adj_matrix = get_adj_matrix_from_distance(new_molecule,coeff)
        new_molecule.set_adj_matrix(adj_matrix)
        is_same = new_molecule.is_same_molecule(original_molecule,False) 
        #print (np.sum(adj_matrix - original_molecule.get_adj_matrix()))
        if is_same:
            break
        coeff += space
    return is_same,coeff


