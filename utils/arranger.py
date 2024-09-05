import numpy as np
import time
import pickle

from autoCG.utils.mapping import ReactionMap


def get_net_reacted_molecules(intermediate_i,intermediate_j):
    """ 
    Returns lists of indices where same molecules within intermediate_i and intermediate_j 
    are removed

    :param intermediate_i,intermediate_j(pyclass 'Intermediate', pyclass 'Intermediate'):

    :return molecule_idx_i,molecule_idx_j(list of integer):
        molecule_indices within the molecule_list intermediates

    """
    molecule_list_i = intermediate_i.molecule_list
    molecule_list_j = intermediate_j.molecule_list
    if molecule_list_i is None:
        molecule_list_i = intermediate_i.get_molecule_list()
    if molecule_list_j is None:
        molecule_list_j = intermediate_j.get_molecule_list()
    len_i = len(molecule_list_i)
    is_molecule_list_i = [True] * len_i
    len_j = len(molecule_list_j)
    is_molecule_list_j = [True] * len_j
    index = 0
    while index < len_i:
        index_prime = 0
        molecule_i = molecule_list_i[index]
        while index_prime < len_j:
            molecule_j = molecule_list_j[index_prime]
            if molecule_i.is_same_molecule(molecule_j,False):
                if is_molecule_list_i[index] and is_molecule_list_j[index_prime]:
                    is_molecule_list_i[index] = False
                    is_molecule_list_j[index_prime] = False
                break
            index_prime += 1
        index += 1
    molecule_idx_i = []
    molecule_idx_j = []
    for i in range(len_i):
        if is_molecule_list_i[i]:
            molecule_idx_i.append(i)
    for i in range(len_j):
        if is_molecule_list_j[i]:
            molecule_idx_j.append(i)
    return molecule_idx_i,molecule_idx_j

def find_permutation_with_respect_to_z(intermediate_i,intermediate_j):
    """ 
    Returns permutation function that well orders the given two atom lists with respect to atomic number
    It matches intermediate_j ordering to intermediate_i ordering
    Ex. If if the atom order of two intermediates are given as: [6,1,1,6,1,1] and [1,1,6,6,1,1]
    Then, the return value, which is the permutation function can have
    {0:1, 1:2, 2:0, 3:3, 4:4, 5:5}, which makes [1,1,6,6,1,1] -> [6,1,1,6,1,1]

    :param intermediate_i, intermediate_j (pyclass 'Intermediate', pyclass 'Intermediate')

    :return permutation(dict):
        permutation function written in dict which matches atom ordering with respect to atomic number (or element)

    """
    permutation = dict()
    formula_by_dict_list_i = intermediate_i.get_formula_as_list()
    formula_by_dict_list_j = intermediate_j.get_formula_as_list()
    for atomic_number in formula_by_dict_list_i:
        list_i = formula_by_dict_list_i[atomic_number]
        list_j = formula_by_dict_list_j[atomic_number]
        if len(list_i) != len(list_j):
            print ('do not have same stoichiometry!!!')
            exit()
        m = len(list_i)
        for k in range(m):
            permutation[list_j[k]] = list_i[k]
    return permutation

def find_z_permutation(z_list1,z_list2):
    """
    This function does same function as 'find_permutation_with_respect_to_z', however inputs are different

    :param z_list1,z_list2(list of integer):
        list of atomic numbers for atom_list

    :return permuation(dict):
        permutation function written in dict which matches atom ordering with respect to atomic number (or element)

    """
    if len(z_list1) != len(z_list2):
        print ('wrong!')
        return None
    type_list1 = dict()
    type_list2 = dict()
    # Find types
    for i in range(len(z_list1)):
        z1 = z_list1[i]
        z2 = z_list2[i]
        if z1 in type_list1:
            type_list1[z1].append(i)
        else:
            type_list1[z1] = [i]
        if z2 in type_list2:
            type_list2[z2].append(i)
        else:
            type_list2[z2] = [i]
    permutation = dict()
    for z in type_list1:
        indices1 = type_list1[z]
        indices2 = type_list2[z]
        for i in range(len(indices1)):
            permutation[indices1[i]] = indices2[i]
    return permutation

def get_permuted_bond_list(bond_list,permutation):
    """
    Returns permuted bond_list permuted by given permutation

    :param bond_list(list of tuple of length 2):
        bond_list that represents the bond connectivity of intermediates

    :param permutation(dict):
        Given one-to-one permutation function

    :return permuted_bond_list(list of tuple of length 2):
        permuted bond list 

    """
    permuted_bond_list = []
    for bond in bond_list:
        permuted_start = permutation[bond[0]]
        permuted_end = permutation[bond[1]]
        #### Must consider order!
        if permuted_start < permuted_end:
            permuted_bond_list.append((permuted_start,permuted_end)) 
        else:
            permuted_bond_list.append((permuted_end,permuted_start))
    return permuted_bond_list


def get_permutation_composition(permutation_1,permutation_2):
    n = len(permutation_1)
    final_permutation = dict()
    for i in permutation_1:
        final_permutation[i] = permutation_2[permutation_1[i]]
    return final_permutation


def get_bond_change(intermediate_i,intermediate_j,num_solution = 1):
    """ 
    Returns the chemical distance defined in Chem.Sci.2018 for given two intermediates intermediate_i and intermediate_j

    :param intermediate_i,intermediate_j(pyclass 'Intermediate', pyclass 'Intermediate'):

    :return final_broken(list of tuple of length 2):
        Bonds that are required to be broken to make intermediate_i to intermediate_j
        Note that bond is not the bond order, which means changing bond order 2 to 1 does not imply that the bond is broken

    :return final_formed(list of tuple of length 2):
        Bonds that are required to be formed to make intermediate_i to intermediate_j
        Note that bond is not the bond order, which means changing bond order 2 to 1 does not imply that the bond is broken

    """
    # First check whether comparison is possible
    if len(intermediate_i.atom_list) != len(intermediate_j.atom_list):
        print ('impossible to calculate CD!!!')
        exit()
    z_list = None
    bond_list_i = None
    bond_list_j = None
    index_function_i = dict()
    index_inverse_function_i = dict()
    index_function_j = dict()
    index_inverse_function_j = dict()
    z_permutation = find_permutation_with_respect_to_z(intermediate_i,intermediate_j)
    z_list_i = intermediate_i.get_z_list()
    z_list_j = intermediate_j.get_z_list()
    bond_list_i = intermediate_i.get_bond_list(False)
    bond_list_j = intermediate_j.get_bond_list(False)
    bond_list_j = get_permuted_bond_list(bond_list_j,z_permutation)
    n = len(intermediate_i.atom_list)
    for i in range(n):
        index_function_i[i] = i
        index_function_j[i] = i
        index_inverse_function_i[i] = i
        index_inverse_function_j[i] = i
    final_permutation = ReactionMap(z_list_i,z_list_i,bond_list_j,bond_list_i,num_solution)
    real_final_permutation = dict()
    #print ('z',z_permutation)
    #print ('final',final_permutation)
    #print ('final_permutation: ', final_permutation)
    inverse_z_permutation = dict()
    inverse_mapping = dict()
    for i in z_permutation:
        inverse_z_permutation[z_permutation[i]] = i
    for i in range(len(final_permutation)):
        real_final_permutation[i] = final_permutation[i]
        inverse_mapping[final_permutation[i]] = i
    bond_list_j = get_permuted_bond_list(bond_list_j,real_final_permutation)
    bond_list_i = set(bond_list_i)
    bond_list_j = set(bond_list_j)
    bond_broken = bond_list_i - bond_list_j # i to j
    bond_formed = bond_list_j - bond_list_i # j to i
    final_broken = get_permuted_bond_list(bond_broken,index_inverse_function_i)
    final_formed = get_permuted_bond_list(bond_formed,index_inverse_function_i)
    real_final_permutation = get_permutation_composition(z_permutation,real_final_permutation)
    return final_broken,final_formed


if __name__ == '__main__':
    from autoCG import chem
    im1 = chem.Intermediate('N.CCl')
    im2 = chem.Intermediate('C[NH3+].[Cl-]')
    #print (chem_dist_calculator.get_chem_distance(im1,im2))
    final_broken,final_formed = get_bond_change(im1,im2)
    print (final_broken,final_formed)
