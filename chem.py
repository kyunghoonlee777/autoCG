'''
--chem.py--
Declaring classes.
(Intermediates, molecules, atoms, ...)
Some parts can be replaced by ASE formats, but now we are using these settings for convenient customization.
'''

from scipy import spatial
import numpy as np
import itertools

from autoCG.utils import ic
from autoCG.utils import process
from autoCG.utils import make_smiles


class Atom:
    """
    :class Atom:
        class Atom mainly contains atomic_number, element, and x,
        Other attributes are not widely used
        molecule_index shows on which molecule that this atom is contained in. For example, if we consider Intermediate C.C, 
        molecule_index can have either 0 or 1, and every atom within atom_list of Intermediate can be assigned by checking
        which molecule does the atom belong to.

    :param data(str or integer):
        Data is either symbol that represents element or atomic number
    
    """
    global periodic_table
    periodic_table = ['H','He','Li','Be','B','C','N','O','F','Ne',\
    'Na','Mg','Al','Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn',\
    'Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y','Zr',\
    'Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I','Xe','Cs','Ba',\
    'La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu',\
    'Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn']

    def __init__(self,data = None):
        self.atomic_number = None
        self.element = None
        self.x = None
        self.y = None
        self.z = None
        if data is not None:
            if type(data) == str:
                self.element = data
            else:
                self.atomic_number = data


    def set_atomic_number(self,atomic_number):
        self.atomic_number = atomic_number
        self.element = periodic_table[atomic_number-1]

    def set_element(self, element):
        """ Type of an atom. e.g. 'C', 'H', 'O', and so on.""" 
        self.element = element
        self.atomic_number = periodic_table.index(element) + 1

    def set_x(self, x):
        """ X-coordinate """ 
        self.x = x

    def set_y(self, y):
        """ Y-coordinate """ 
        self.y = y

    def set_z(self, z):
        """ Z-coordinate """ 
        self.z = z

    def set_coordinate(self, position):
        """ Set Cartesian coordinates of an atom """ 
        dim = len(position)
        if dim == 2:
            x = position[0]
            y = position[1]
            z = 0
        elif dim == 3:
            x = position[0]
            y = position[1]
            z = position[2]
        self.x = x
        self.y = y
        self.z = z



    def get_atomic_number(self):
        """
        Returns the atomic number (number of protons) of a given atom.

        :param:

        :return integer:
            Directly the atomic number of a given atom is returned
        """
        if self.atomic_number == None:
            element = self.element
            if element == None:
                print ('atom is not specified!')
            if len(element)>1:
                end_part = element[1:]
                end_part = str.lower(end_part)
                element = element[0]+end_part
                self.element = element
            if element in periodic_table:
                index = periodic_table.index(element)
                self.atomic_number = index + 1
            else:
                print ('element',element)
                print ('modify periodic table!!!')
        return self.atomic_number


    def get_element(self):
        """
        Returns symbol of a given atom.

        :param:

        :return str:
            Directly the symbol of a given atom is returned
        """
        if self.element is None:
            atomic_number = self.atomic_number
            if atomic_number is None:
                print ('atom is not specified!')
            z = int(self.atomic_number)-1
            self.element = periodic_table[z]
        return self.element

    def get_x(self):
        return self.x

    def get_y(self):
        return self.y

    def get_z(self):
        return self.z

    
    def get_coordinate(self):
        return np.array([self.x,self.y,self.z])    

    
    def get_radius(self):
        """
        Returns a radius information of a given atom. Reference is given here: Dalton Trans., 2008, 2832-2838
        
        :param:

        :return float:
            It directly returns the reference values
        """
        element = self.get_element()       
        a = str.lower(element)
        #reference : Dalton Trans., 2008, 2832-2838
        if a=='h': 
            return 0.31
            #return 0.38
        elif a=='li': 
            return 1.28
        elif a=='be': 
            return 0.96
        elif a=='b': 
            return 0.84
        elif a=='c': 
            return 0.76
        elif a=='n': 
            return 0.71
        elif a=='o': 
            return 0.66
        elif a=='f': 
            return 0.57
        elif a=='na': 
            return 1.66
        elif a=='mg': 
            return 1.41
        elif a=='al': 
            return 1.21
        elif a=='si': 
            return 1.11
        elif a=='p': 
            return 1.07
        elif a=='s': 
            return 1.05
        elif a=='cl': 
            return 1.02
        elif a=='ar': 
            return 0.76
        elif a=='k': 
            return 2.03
        elif a=='ca': 
            return 1.76
        elif a=='co': 
            return 1.38 #1/2*(lowspin+highspin)
        #elif a=='co': return 1.26 #lowspin
        #elif a=='co': return 1.50 #highspin
        elif a=='fe': 
            return 1.42 #1/2*(lowspin+highspin)
        elif a=='ni': 
            return 1.24
        #elif a=='cr': return 1.39
        elif a=='ti': 
            return 1.60
        elif a=='br': 
            return 1.20
        elif a=='rh': 
            return 1.42
        elif a=='pd': 
            return 1.39
        elif a=='i': 
            return 1.39
        elif a=='hf': 
            return 1.75
        else: 
            return 0

        #reference : J. Chem. Phys. 41, 3199 (1964)
        '''
        if a=='h': return 0.25
        elif a=='li': return 1.45
        elif a=='be': return 1.05
        elif a=='b': return 0.85
        elif a=='c': return 0.70
        elif a=='n': return 0.65
        elif a=='o': return 0.60
        elif a=='f': return 0.50
        elif a=='na': return 1.80
        elif a=='mg': return 1.50
        elif a=='al': return 1.25
        elif a=='si': return 1.10
        elif a=='p': return 1.00
        elif a=='s': return 1.00
        elif a=='cl': return 1.00
        elif a=='ar': return 0.71
        elif a=='k': return 2.20
        elif a=='ca': return 1.80
        elif a=='co': return 1.35
        else: return 0
        '''

    def get_period_group(self):
        """
        Returns a period,group information from a given atom. It finds period, group by identifying 
        electron configuration (orbital configuration)
        
        :param:

        :return period,group(int,int):
            Note that these values are actual period and group. If C atomm is given, it returns 2,4

        """
        atomic_number = self.get_atomic_number()
        num_of_electrons = atomic_number
        # Orbital: [n,l,num_of_electrons]
        sum_of_n_and_l=1
        orbital_configuration=[]
        while num_of_electrons>0:
            # Generate orbitals within sum_of_n_and_l=k
            # New orbitals are introduced for (k+1)/2
            maximal_l=int((sum_of_n_and_l-1)/2)
            for l in range(maximal_l,-1,-1):
                # Start with lowest l
                if num_of_electrons>4*l+2:
                    num_of_electrons-=(4*l+2)
                    orbital_configuration.append([sum_of_n_and_l-l,l,4*l+2])
                else:
                    orbital_configuration.append([sum_of_n_and_l-l,l,num_of_electrons])
                    num_of_electrons=0
                    break
            sum_of_n_and_l+=1
        # Get maximal n and l
        period=0
        for orbital in orbital_configuration:
            if orbital[0]>period:
                period = orbital[0]
        # If transition metal, we add 9, Sc has group 9 for else, we do not consider ...
        last_orbital = orbital_configuration[-1]
        if last_orbital[1]<2:
            group = 2 * last_orbital[1] ** 2 + last_orbital[2]
        else:
            group = 8 + last_orbital[2]
        return period,group

    def get_max_valency(self): 
        """
        Returns maximal valency of a given atom. Examples of those values are shown in the main code.  

        :param:

        :return integer(nonnegative integer):
            possible maximal valency
        """
        element = self.get_element()
        a = str.lower(element)
        if a=='c':   
            return 4
        elif a=='si':
            return 4
        elif a=='h': 
            return 1
        elif a=='be':
            return 2
        elif a=='b': 
            return 3
        #elif a=='b': return 4
        elif a=='o': 
            #return 2 # 3 is also possible...
            return 3
        elif a=='n': 
            return 4
        elif a=='li': 
            return 1
        #elif a=='p': return 4
        #elif a=='s': return 4
        elif a=='p': 
            return 5 #valence shell expansion
        elif a=='s': 
            return 6 #valence shell expansion
        elif a=='f': 
            return 1
        elif a=='na': 
            return 1
        elif a == 'mg':
            return 2
        elif a=='co':
            return 6
        elif a=='rh':
            return 6
        elif a=='ni':
            return 6 
        elif a=='ti':
            return 6
        elif a=='fe':
            return 6
        elif a=='cl': 
            return 1
        elif a=='br': 
            return 1
        elif a=='bb': 
            return 3
        elif a=='lg': 
            return 2
        elif a=='pd': 
            return 6
        elif a=='i': 
            return 3

    def copy(self):
        new_atom = Atom()
        # Copy all attributes
        new_atom.atomic_number = self.atomic_number
        new_atom.element = self.element
        new_atom.x = self.x
        new_atom.y = self.y
        new_atom.z = self.z
        return new_atom 

    def is_same_atom(self,atom):
        """
        Returns whether the two atoms have same type by comparing atomic number

        :param atom(pyclass 'Atom'):
            Our defined class 'Atom'

        :return:
            True: Two atoms are same type
            False: Two atoms are different type
        """
        atomic_number = self.get_atomic_number()
        atomic_number_prime = atom.get_atomic_number()
        return atomic_number == atomic_number_prime

    def get_content(self,option='element',criteria = 1e-4):
        x = self.x
        y = self.y
        z = self.z
        if abs(x) < criteria:
            x = 0.00
        if abs(y) < criteria:
            y = 0.00
        if abs(z) < criteria:
            z = 0.00
        content = ' ' + str(x) + ' ' + str(y) + ' ' + str(z) + '\n'
        if option=='element':
            content = self.get_element() + content
        else:
            content = str(self.get_atomic_number()) + content
        return content

    def __eq__(self,atom):
        return self.is_same_atom(atom)



#bookmark - Jinwon Lee
         
class Molecule:
    """
    :class Molecule:
        class Molecule mainly contains atom_list, atom_feature, chg, bo_matrix, adj_matrix, energy, smiles, c_eig_list
        atom_list is a list of atom (pyclass 'Atom')
        atom_feature is a dict, where features of atom such as formal charge, number of pi bonds, etc, are stored.
        Those information can be freely added by giving new dict
        c_eig_list can be used to identify whether the given two molecules are the same. This c_eig_list is invariant
        to the permutation of atom indexing, identification between any two generated molecules can be easily checked.

    :param data(str or xyz file or None):
        data should be feeded as either smiles(str), xyz file(file) or None
        * If smiles is used, rd_mol is generated using the given smiles and converted into ace_mol (pyclass 'Molecule)
        * If xyz file is used, directly the 3d geometry of the given molecule is generated. If you want to generate adj_matrix,
        bo_matrix, chg_list, etc, refer following method contained in pyclass 'Molecule' (charge of molecule should be given!!!)
        i. Generate adj_matrix by using 'get_adj_matrix_from_distance' stored in class 'Molecule'
        ii. Generate bo_matrix by using 'get_adj_matrix_from_adj_matrix' stored in class 'Molecule'
        iii. Then, using those bo_matrix, get chg_list by using 'get_chg_list_from_bo_matrix' stored n class 'Molecule'
        * If None is used, only blank virtual molecule is generated

    """
    def __init__(self,data = None):
        self.atom_list = []
        self.atom_feature = dict()
        self.adj_matrix = None
        self.bo_matrix = None
        self.chg = None
        self.multiplicity = None
        self.energy = None        
        self.smiles = None
        self.c_eig_list = None
        self.formula_id = None
        self.molecule_id = None
        self.atom_id_list = None
        
        if data == None:
            pass

        elif type(data) == str:
            if data[-4] == '.':
                conformer = Conformer(data)
                # At least make adjacency
                self.atom_list = conformer.atom_list
                self.chg = conformer.chg
                self.multiplicity = conformer.multiplicity
                self.energy = conformer.energy
                self.adj_matrix = conformer.get_adj_matrix()
            else:
                # Generate data with rdkit
                from rdkit import Chem            
                try:
                    rd_mol = Chem.MolFromSmiles(data,sanitize = False)
                    Chem.SanitizeMol(rd_mol,sanitizeOps = Chem.SanitizeFlags.SANITIZE_ALL^Chem.SanitizeFlags.SANITIZE_ADJUSTHS)
                    ace_mol = process.get_ace_mol_from_rd_mol(rd_mol)
                except:
                    from openbabel import pybel
                    ob_mol = pybel.readstring('smi',data)
                    ace_mol = process.get_ace_mol_from_ob_mol(ob_mol)
                self.atom_list = ace_mol.atom_list
                self.atom_feature = ace_mol.atom_feature
                self.adj_matrix = ace_mol.adj_matrix
                self.bo_matrix = ace_mol.bo_matrix
                self.chg = np.sum(self.atom_feature['chg'])
                self.smiles = data
          
        else:
            try:
                len(data[0]) # Then, it's some kind of list,tuple,numpy
                make = True
            except:
                make = False
                #print ('Wrong input type!!!')
            if make:
                # Data: (z_list,adj_matrix,bo_matrix,chg_list)
                atom_list = []
                atom_info_list = data[0]
                for atom_info in atom_info_list:
                    atom = Atom(atom_info)
                    atom_list.append(atom)
                self.atom_list = atom_list
                self.adj_matrix = data[1]
                self.bo_matrix = data[2]
                self.atom_feature['chg'] = data[3]
                if self.adj_matrix is None and self.bo_matrix is not None:
                    self.adj_matrix = np.where(self.bo_matrix>0,1,0)
            else:
                #print ('Something wrong with the input data!!! Check your data again!')
                a = 1
             
        #if data != None:
        #    self.formula_id = int(np.sum(self.get_z_list()**3))

    def set_adj_matrix(self,adj_matrix):
        self.atom_feature['chg'] = None
        if adj_matrix is not None:
            self.adj_matrix = np.copy(adj_matrix)
        else:
            self.adj_matrix = None
        self.bo_matrix = None

    def set_bo_matrix(self,bo_matrix):
        self.atom_feature['chg'] = None
        self.bo_matrix = np.copy(bo_matrix)
        self.adj_matrix = np.where(self.bo_matrix>0,1,0)

    
    def get_adj_matrix(self):
        if self.adj_matrix is not None:
            return self.adj_matrix
        if self.bo_matrix is not None:
            adj_matrix = np.where(self.bo_matrix>0,1,0)
            return adj_matrix
        return None
    
    def get_bo_matrix(self):            
        return self.bo_matrix
    
    def get_chg(self):
        if self.chg is None:
            chg_list = self.get_chg_list()
            if chg_list is None:
                return None
            else:
                return int(np.sum(chg_list))
        else:
            return int(self.chg)

    def get_multiplicity(self):
        if self.multiplicity is None:
            try:
                chg = self.get_chg()
            except:
                print ('Total charge is not provided! At least charge information!')
                return None
            try:
                e_list = self.get_num_of_lone_pair_list()
                num_electron = len(np.where((2*e_list) % 2 == 1)[0])    
            except:
                z_sum = np.sum(self.get_z_list())
                num_electron = z_sum - chg
            multiplicity = num_electron % 2 + 1
            return int(multiplicity)
        else:
            return int(self.multiplicity)


    def get_energy(self):
        return self.energy


    def get_c_eig_list(self,c_sum = False):
        """
        Returns the eigenvalues of coulomb matrix (see the instruction 'get_matrix' within class 'Molecule'.

        :param c_sum(boolean):
            If True, it also returns the sum of eigenvalues

        :return c_eig_list(pyclass 'numpy.ndarray' or pyclass 'numpy.ndarray',int(c_matrix)):
            eigenvalue of coulomb matrix of a given molecule and also it's sum if c_sum = True

        """
        if self.c_eig_list is None:        
            c_matrix = self.get_matrix('coulomb')
            #### Check symmetry of matrix
            check_sum = np.sum(np.abs(c_matrix-np.transpose(c_matrix)))
            if check_sum > 0.01:
                print ('something wrong')
            c_eig_list,vec = np.linalg.eig(c_matrix)
            if np.isnan(c_eig_list[0]):
                print ('nanvalue!!!')
                import sys;sys.exit()
            if c_sum:
                return np.sort(c_eig_list.real), int(np.sum(c_matrix))
            else:
                return np.sort(c_eig_list.real)
        else:
            return self.c_eig_list

    def get_radius_list(self):
        """
        Returns radius list of a given molecule
        Ex. For CH4 with atom order [C,H,H,H,H], then radius_list is given as
        [0.8,0.5,0.5,0.5,0.5], if radius of C and H are given as 0.8 and 0.5 (unit is Angstrom)

        :param:

        :return radius_list(list of float):
            list of radius of each atom
        """
        atom_feature = self.atom_feature
        if 'radius' in atom_feature:
            return atom_feature['radius']
        radius_list = list(map(lambda x:x.get_radius(),self.atom_list))
        return radius_list

    def get_smiles(self, method="none", find_stereocenter="N"):
        """
        Returns smiles of a given molecule. It could return different smiles according to the atom order.

        :param method(str):
            method can either be 'ace' or 'rdkit'
            If 'ace' is used, smiles is generated throughout our defined method
            If 'rdkit' is used, it generates smiles using rdkit module

        :param find_stereocenter(str):
            It can have either 'Y' or 'N'. If 'Y'(yes), stereo information such as E/Z or R/S are desginated in the smiles

        :return smiles(str):
            Returns a smiles of given molecule generated by input option
        """
        from rdkit import Chem

        if method == "ace":
            # Implement ACE methods
            atom_list = self.atom_list
            bo_matrix = self.get_matrix("bo")
            fc_list = self.get_chg_list()
            bond_list = self.get_bond_list(False)
            if bo_matrix is None:
                adj_matrix = self.get_matrix("adj")
                if adj_matrix is None or self.chg is None:
                    print("We need to know both adjacency and charge!!!")
                    return None
                else:
                    bo_matrix = process.get_bo_matrix_from_adj_matrix(self, self.chg)
                    fc_list = process.get_chg_list_from_bo_matrix(
                        self, self.chg, bo_matrix
                    )
            elif fc_list is None:
                fc_list = process.get_chg_list_from_bo_matrix(self, self.chg, bo_matrix)
            fc_list = fc_list.tolist()
            # Check element is ready!
            for atom in atom_list:
                atom.set_element(atom.get_element())
            smiles_list = make_smiles.GetSMILES(
                atom_list, bo_matrix, bond_list, fc_list, find_stereocenter
            )
            return smiles_list[0]
        elif method == "rdkit":
            if self.smiles != None:
                return self.smiles
            else:
                rd_mol = self.get_rd_mol()
                SMILES = Chem.MolToSmiles(rd_mol)
                return SMILES
        elif self.smiles != None:
            return self.smiles
        else:
            return None


            
    def get_matrix(self,type_of_matrix = 'bo'):
        """
        Returns a matrix that contains some information of a given molecule that is widely used. 
        Distance matrix can be easily evaluated by using rdkit module. 
        :param type_of_matrix(str):
            Mainly 4 inputs are possible: 'bo', 'adj', 'coulomb', 'distance'
            'bo' corresponds to bond order matrix (B_{ij})
            'adj' corresponds to adjacency matrix whose elements (A_{ij}) are 1 for bonded i,jth atom and 0 for nonbonded i,jth atom 
            'coulomb' corresponds to the coulomb matrix whose elements are represented as C_{ij} = A_{ij}*Z_{i}*Z_{j}
            'distance' corresponds to the graphical distance matrix (not actual distance matrix). 
            D_{ij} is the length of shortest path between ith atom and jth atom

        :return matrix(pyclass numpy.ndarray):
         
        """
        atom_list = self.atom_list
        if type_of_matrix == 'bo':
            return self.get_bo_matrix()
        elif type_of_matrix == 'adj':
            return self.get_adj_matrix()
        elif type_of_matrix == 'coulomb':
            adj_matrix = self.get_adj_matrix()
            z_list = self.get_z_list()
            if z_list is None:
                return None
            new_adj_matrix = adj_matrix + np.diag([1]*len(z_list))
            diagonal_matrix = np.diag(z_list)
            coulomb_matrix = np.matmul(np.matmul(diagonal_matrix,new_adj_matrix),diagonal_matrix)
            return coulomb_matrix
        elif type_of_matrix == 'distance':
            adj_matrix = self.get_adj_matrix()
            n = len(self.atom_list)
            distance_matrix = 1000 * np.ones((n,n))
            np.fill_diagonal(distance_matrix,0)
            update_matrix = np.identity(n)
            for d in range(n):
                cur_d = d + 1
                update_matrix = np.matmul(update_matrix,adj_matrix)
                d_update_matrix = cur_d * np.where(update_matrix>0,1,1000)
                indices = np.where(d_update_matrix<distance_matrix)
                distance_matrix = np.where(d_update_matrix<distance_matrix,d_update_matrix,distance_matrix)
                if len(indices[0]) == 0:
                    break
            return distance_matrix

    def check_matrix(self,type_of_matrix='bo'):
        if type_of_matrix == 'adj':
            if self.get_adj_matrix() is None:
                return False
            return True
        if type_of_matrix == 'bo':
            if self.get_bo_matrix() is None:
                print ('bo matrix is not prepared!!! It is necessary to define bo matrix')
                return False
            return True

    def get_distance_matrix(self):
        return self.get_matrix('distance')


    def __eq__(self,molecule):
        return self.is_same_molecule(molecule,True)

       
    def is_same_molecule(self,molecule,option=False):
        """
        Checks whether the given two molecules are the same by comparing the c_eig (see instructions in 'get_c_eig_list' within 
        class 'Molecule'
            
        :param molecule(pyclass 'Molecule):
            
        :return True/False(boolean):
            True: Two molecules are the same
            False: Two molecules are different
        """
        if len(self.atom_list) != len(molecule.atom_list):
            return False
        if option:
            if self.get_chg() != molecule.get_chg():
                return False
            #elif self.get_multiplicity() != molecule.get_multiplicity():
            #    return False
        '''
        id1 = self.get_molecule_id()
        id2 = molecule.get_molecule_id()
        if id1 != id2:
            return False
        else:
            return True
        '''
        c_eig_list1 = self.get_c_eig_list()
        c_eig_list2 = molecule.get_c_eig_list()
        delta_c_eig_list = np.abs(c_eig_list1-c_eig_list2)
        total_delta = np.sum(delta_c_eig_list)
        return total_delta < 1e-8 

    def copy(self,copy_all = False):
        new_molecule = Molecule()
        atom_list = self.atom_list
        # First copy atoms
        new_atom_list = []
        for atom in atom_list:
            new_atom_list.append(atom.copy())
        new_molecule.atom_list = new_atom_list
        # Copy connectivity information
        bo_matrix = self.get_matrix('bo')
        if bo_matrix is not None:
            new_molecule.bo_matrix = np.copy(bo_matrix)
        else:
            adj_matrix = self.get_matrix('adj')
            if adj_matrix is not None:
                new_molecule.adj_matrix = np.copy(adj_matrix)
            else:
                print ('Warning: Connectivity information is not included in the molecule!!!')
        # Finally, copy charge
        new_molecule.chg = self.get_chg()
        #new_molecule.multiplicity = self.get_multiplicity()
        if 'chg' in self.atom_feature:
            new_molecule.atom_feature['chg'] = np.copy(self.atom_feature['chg'])
        # Above things are essential for copy
        if copy_all:
            # Copy other attributes
            new_molecule.energy = self.energy
            new_molecule.smiles = self.smiles
            new_molecule.c_eig_list = self.c_eig_list
        return new_molecule

    def get_z_list(self):
        """
        Returns atomic number list of a given molecule
         
        :param:

        :return z_list(pyclass 'numpy.ndarray'):
            list of atomic number
        """
        atom_feature = self.atom_feature
        if atom_feature is not None and 'atomic number' in atom_feature:
            return atom_feature['atomic number']
        else:
            z_list = list(map(lambda x:x.get_atomic_number(),self.atom_list))
            return np.array(z_list)

    def get_element_list(self):
        """
        Returns element list of a given molecule
         
        :param:

        :return element_list(list of string):
            list of element written as capital letters
        """
        atom_feature = self.atom_feature
        if atom_feature is not None and 'element' in atom_feature:
            return atom_feature['element']
        else:
            element_list = list(map(lambda x:x.get_element(),self.atom_list))
            return element_list

    def get_group_list(self):
        atom_feature = self.atom_feature
        if 'group' in atom_feature:
            return atom_feature['group']
        else:
            group_list = list(map(lambda x:x.get_period_group()[1],self.atom_list))
            return np.array(group_list)

    def get_period_list(self):
        atom_list = self.atom_list
        atom_feature = self.atom_feature
        if 'period' in atom_feature:
            return atom_feature['period']
        else:
            period_list = list(map(lambda x:x.get_period_group()[0],self.atom_list))
            return np.array(period_list)

    def get_period_group_list(self):
        """ 
        Returns period_group_list for atoms within molecule. See details for obtaining those period/group
        in 'get_period_group' defined within class 'Atom'
        
        :param:

        :return period_list,group_list(pyclass 'numpy.ndarray',pyclass 'numpy.ndarray'):

        """
        atom_list = self.atom_list
        atom_feature = self.atom_feature
        if 'period' in atom_feature and 'group' in atom_feature:
            return atom_feature['period'],atom_feature['group']
        n = len(atom_list)
        period_list = []
        group_list = []
        for i in range(n):
            period,group = atom_list[i].get_period_group()
            period_list.append(period)
            group_list.append(group)
        return np.array(period_list),np.array(group_list)

    def get_mass_list(self):
        atom_feature = self.atom_feature
        if 'mass' in atom_feature:
            return atom_feature['mass']
        else:
            mass_list = list(map(lambda x:x.get_mass(),self.atom_list))
            return np.array(mass_list)

    def get_chg_list(self):
        """
        Returns chg list of a given molecule
         
        :param:

        :return chg_list(pyclass 'numpy.ndarray'):
            list of formal charge
        """
        atom_feature = self.atom_feature
        try:
            return atom_feature['chg']
        except:
            #print ('charge lists are not prepared!!!')
            return None
    
    def get_valency_list(self):
        """
        Returns the valency list of a given molecule. Valency is defined as the number of neighbors 
        (in otherwords, the number of atoms that are bonded with each other, therefore uses adj_matrix to compute valency_list)

        :param:

        :return valency_list(list of integer):
            list of valency of all atoms within a given molecule
        """
        adj_matrix = self.get_adj_matrix()
        if adj_matrix is None:
            print ('You need to define adj matrix!!!')
            return None
        valency_list = np.sum(adj_matrix,axis = 1)
        return valency_list

    def get_total_bond_order_list(self):
        """ 
        Returns the list of number of bonds for each atoms. For example of C2H4 ethene, 
        the valency of C is 3, since it is bonded to H, H, C. However, for total bond order, its value is 4, 
        since C-C is double bond. If atom order is given as [C,C,H,H,H,H], each valency, total_bond_order is given as
        valency_list: [3,3,1,1,1,1]
        total_bond_orer_list: [4,4,1,1,1,1]

        :param:

        :return total_bond_order_list(list of integer):
            list of bond order of all atoms within a given molecule

        """
        atom_list = self.atom_list
        bo_matrix = self.get_bo_matrix()
        if bo_matrix is None:
            print ('We cannot obtain bond order, since the bond order matrix is not prepared!!!')
            return None
        n = len(atom_list)        
        total_bond_order_list = np.sum(bo_matrix,axis = 1)
        return total_bond_order_list

    def get_max_valency_list(self):
        """
        Returns max_valency_list for atoms within molecule. Max valency is defined in 'get_max_valency' within class 'Atom'

        :param:

        :return max_valency_list(list of integer):
            list of maximal valency of atoms within molecule 
        """
        atom_feature = self.atom_feature
        if atom_feature is not None and 'max valency' in atom_feature:
            return atom_feature['max valency']
        else:
            ### Get default max valency
            atom_list = self.atom_list
            n = len(atom_list)
            max_valency_list = []
            for i in range(n):
                atom = atom_list[i]
                max_valency_list.append(atom.get_max_valency())
            return np.array(max_valency_list)


    def get_bond_list(self,contain_bond_order = True):
        """
        Returns the total bond list as list of tuples
        For example, if CH4 is given with atom order [C,H,H,H,H], if contain_bond_order = False, the output is given as
        [(0,1),(0,2),(0,3),(0,4)]
        if contain_bond_order = True, the output is given as
        [(0,1,1),(0,2,1),(0,3,1),(0,4,1)]

        :param contain_bond_order(boolean):
            If contain_bond_order is False, it only returns bonds represented as (i,j) between atoms within the given intermediate.
            If contain_bond_order is True, it returns bonds with bond order included (i,j,bond_order), where bond_order can only have 1,2,3.
            Therefore, a given molecule(self) should be kekulized.

        :return bond_list(either list of tuple with size 2 or 3):
            bond_list
        """
        atom_list = self.atom_list
        n = len(atom_list)
        check_matrix = self.bo_matrix
        total_bond_list = []
        if contain_bond_order:
            check_matrix = self.get_matrix('bo')
            if check_matrix is None:
                print ('we cannot give bond order!!!')
                print ('We will automatically give only bond list!')
                contain_bond_order = False
                check_matrix = self.get_matrix('adj')
                if check_matrix is None:
                    print ('matrix',self.atom_list)
                    print ('hahahahahaha',check_matrix)
                    print ('Give connectivity! We cannot find the bond!')
                    return None
        if contain_bond_order:            
            bond_type = [1,2,3]
            check_matrix = self.get_matrix('bo')
        else:
            bond_type = [1]
            check_matrix = self.get_matrix('adj')
        #check_matrix = self.adj_matrix
        for bond_order in bond_type:
            bond_list = np.where(check_matrix == bond_order)
            bond_list = np.stack(bond_list,axis = 1)
            for array in bond_list:
                if array[0] < array[1]:
                    if contain_bond_order:
                        bond_tuple = (int(array[0]),int(array[1]),int(bond_order))
                    else:
                        bond_tuple = (int(array[0]),int(array[1]))
                    total_bond_list.append(bond_tuple) 
        return total_bond_list 

    def get_formula_as_list(self):
        """
        Returns stoichiometry (chemical formula) of a given molecule including atom indices
        For example, for given CH4 molecule with atom ordering [C,H,H,H,H], it returns dict form
        {C:[0],H:[1,2,3,4]}
        For C2H4 molecule with atom ordering [C,H,H,C,H,H], it returns dict form
        {C:[0,3],H:[1,2,4,5]}
        This function is used in arranger for evaluating chemical distance

        :param:
       
        :return element_idx_list(dict):
            chemical formula with atom indices

        """
        atom_list = self.atom_list
        element_idx_list = dict()
        if atom_list is None:
            print ('No atoms!!! We cannot get formula!')
        for i in range(len(atom_list)):
            atom = atom_list[i]
            element_type = atom.get_element()
            if element_type in element_idx_list:
                element_idx_list[element_type].append(i)
            else:
                element_idx_list[element_type] = [i]
        return element_idx_list


    def get_rd_mol(self,atom_stereos = dict(), bond_stereos = dict()):
        """
        Returns molecule with type pyclass 'rdkit.Chem.rdchem.Mol' from our type pyclass 'Molecule'
        Note that atom ordering and bond order is well preserved

        :param include_stereo(boolean):
            Do not touch this option, we have not develop options for molecule that considers stereocenter 
        
        :return rd_mol(pyclass 'rdkit.Chem.rdchem.Mol'):

        """
        from rdkit import Chem
        rd_mol = None
        bond_types = {1:Chem.BondType.SINGLE, 2:Chem.BondType.DOUBLE, 3:Chem.BondType.TRIPLE}
        #bond_stereo = {0:Chem.rdchem.BondStereo.STEREONONE,1:Chem.rdchem.BondStereo.STEROE,2:Chem.rdchem.BondStereo.StEREOZ}
        rd_mol = Chem.Mol()
        rde_mol = Chem.EditableMol(rd_mol)
        atom_list = self.atom_list
        #print("atom_list:", list(map(lambda x: x.get_element(), self.atom_list)))
        atom_feature = self.atom_feature
        chg_list = self.get_chg_list()
        #if chg_list is None:
        #    print ('Charge information is not included in the RDKit molecule!')
        n = len(atom_list)
        # Add atoms
        for i in range(n):
            atom = atom_list[i]
            rd_atom = Chem.Atom(int(atom.get_atomic_number()))
            if chg_list is not None:
                rd_atom.SetFormalCharge(int(chg_list[i]))
            rde_mol.AddAtom(rd_atom)
        # Add bonds
        bond_list = self.get_bond_list(True)
        for bond in bond_list:
            if bond[0] < bond[1]:
                rde_mol.AddBond(bond[0],bond[1],bond_types[bond[2]])
        rd_mol = rde_mol.GetMol()
        Chem.SanitizeMol(rd_mol)
        return rd_mol

    def get_ob_mol(self,include_stereo=False):
        """
        Returns molecule with type pyclass 'openbabel.OBMol' from our type pyclass 'Molecule'
        Note that atom ordering and bond order is well preserved
        :param include_stereo(boolean):
            Do not touch this option, we have not develop options for molecule that considers stereocenter 
        
        :return rd_mol(pyclass 'openbabel.OBMol'):
        """
        from openbabel import openbabel

        ob_mol = openbabel.OBMol()
        atom_list = self.atom_list
        bond_list = self.get_bond_list(True)
        n = len(atom_list)
        z_list = self.get_z_list()
        chg_list = None
        atom_feature = self.atom_feature
        if atom_feature is not None and 'chg' in atom_feature:
            chg_list = atom_feature['chg']
        # Generate atoms
        for i in range(n):
            atom = atom_list[i]
            ob_atom = ob_mol.NewAtom()
            z = int(atom.get_atomic_number())
            ob_atom.SetAtomicNum(z)
            ob_atom.SetFormalCharge(int(chg_list[i]))
        # Generate bonds
        for bond in bond_list:
            ob_mol.AddBond(bond[0]+1,bond[1]+1,bond[2])                
        '''
        else:
            import pybel
            ob_mol = pybel.readstring('smi',self.smiles)
        '''
        return ob_mol

    def sanitize(self):
        adj_matrix = self.get_adj_matrix()
        if adj_matrix is None:
            print ('Cannot sanitize molecule because there is no adj matrix information !!!')
        else:
            chg = self.get_chg()
            if chg is None:
                print ('Cannot sanitize molecule because there is no charge information !!!')
            bo_matrix = process.get_bo_matrix_from_adj_matrix(self,chg)
            fc_list = process.get_chg_list_from_bo_matrix(self,chg,bo_matrix)
            self.bo_matrix = bo_matrix
            self.atom_feature['chg'] = fc_list


    def get_valid_molecule(self):
        #bo_matrix = self.get_matrix('bo')
        #chg_list = self.get_chg_list()
        z_list = self.get_z_list()
        original_z_list = np.copy(z_list)
        n = len(z_list)
        chg = self.get_chg()
        if chg is None:
            virtual_chg = 0
        else:
            virtual_chg = chg
        """
        if bo_matrix is None:
            chg = self.chg
            if chg is None:
                chg = 0
            bo_matrix = process.get_bo_matrix_from_adj_matrix(self,chg)
            chg_list = process.get_chg_list_from_bo_matrix(self,chg,bo_matrix)
        """
        #n = len(chg_list)
        period_list, group_list = self.get_period_group_list()
        adj_matrix = self.get_adj_matrix()
        adjacency_list = np.sum(adj_matrix, axis=0)
        
        #Compute SN
        problem_indices = np.flatnonzero(adjacency_list > self.get_max_valency_list())
        #print("z_list", z_list)
        #print("problem_indices", problem_indices)
        new_z_list = np.copy(z_list)
        for idx in problem_indices:
            period = period_list[idx]
            group = group_list[idx]
            adj = adjacency_list[idx]
            if period == 1:
                if adj > 1:
                    new_z_list[idx] = 10 - adj
            else:
                if adj < 5:
                    if period == 2:
                        new_z_list[idx] = 10 - adj
                    else:
                        new_z_list[idx] = 18 - adj
                else:
                    new_z_list[idx] = 10 + adj #replace with one higher period element with proper valency
        #Construct new Molecule
        virtual_molecule = Molecule([new_z_list, adj_matrix, None, None])

        #Construct BO and Chg
        new_bo_matrix = process.get_bo_matrix_from_adj_matrix(virtual_molecule, virtual_chg)
        new_bo_sum = np.sum(new_bo_matrix, axis=0)

        new_period_list, new_group_list = virtual_molecule.get_period_group_list()
        # Modify new_z for atoms containing unpaired electron ...
        for i in range(n):
            group = new_group_list[i]
            bo = new_bo_sum[i]
            if group % 2 != bo % 2:
                new_group_list[i] += 1
                virtual_molecule.atom_list[i].atomic_number += 1

        for i in range(n):
            group = new_group_list[i]
            bo = new_bo_sum[i]
            parity = 0
            # new_z: Modified z when overvalence is observed. = 10-bo for period=2, = 18 - bo (bo <=4), = bo + 10 (bo>4)
            if new_period_list[i] == 1:
                octet = 1 # For checking validity 
                new_z = 10 - bo
            elif new_period_list[i] == 2:
                octet = min(group,4)
                new_z = 10 - bo
            else:
                octet = group
                if bo > 4:
                    new_z = bo + 10 # Reconsider valence expansion
                else:
                    new_z = 18 - bo # Just same with F, O, N substitution, but with higher order
            if not process.check_atom_validity(group,bo,0,octet): # Set every charge equal to zero 
                virtual_molecule.atom_list[i].set_atomic_number(int(new_z))
        virtual_molecule.set_bo_matrix(new_bo_matrix)
        virtual_molecule.atom_feature['chg'] = np.zeros((n)) # Set charge zero
        new_z_list = virtual_molecule.get_z_list()
        return virtual_molecule


    def make_3d_coordinate(self,library='rdkit',atom_stereos = dict(),bond_stereos = dict()):
        return self.make_3d_coordinates(1,library)[0]


    def make_3d_coordinates(self,num_conformer = 1,library = 'rdkit',atom_stereos = dict(),bond_stereos = dict()):
        """
        Returns possible 3d molecular geometry using other libraries, mainly 'babel' and 'rdkit'

        :param library(str):
            Either 'babel' or 'rdkit' are possible. You can add your own found libraries for generating 3d structure 
        
        :return coordinate_list(list of tuples with size 3(float)):
            3d geometry of molecule generated by other libraries

        """
        from rdkit.Chem import AllChem
        from rdkit import Chem
        mol = None
        coordinates = []
        changed_indices = []
        n = len(self.atom_list)
        if library == 'rdkit':
            params = Chem.rdDistGeom.srETKDGv3()
            #params.pruneRmsThresh = 0.25
            try:
                mol = self.get_rd_mol(atom_stereos = atom_stereos, bond_stereos = bond_stereos)
                Chem.SanitizeMol(mol)
                mol = Chem.AddHs(mol)
            except:
                virtual_molecule,changed_indices = self.get_valid_molecule()
                mol = virtual_molecule.get_rd_mol(atom_stereos = atom_stereos, bond_stereos = bond_stereos)
                Chem.SanitizeMol(mol)
                mol = Chem.AddHs(mol)
                print("in making 3d coord", Chem.MolToSmiles(mol))
                print("SMILES for Hypothetical TS:\t", Chem.MolToSmiles(mol))
            if mol is None:
                print ('Impossible embedding')
                return []
            conformer_id_list = Chem.rdDistGeom.EmbedMultipleConfs(mol, num_conformer, params)
            #print (conformer_id_list)
            conformer_energy_list = dict()
            converged_conformer_id_list = []
            for conformer_id in conformer_id_list:
                converged = not AllChem.UFFOptimizeMolecule(mol, confId=conformer_id)
                #if converged: converged_conformer_id_list.append(conformer_id)
                if converged: converged_conformer_id_list.append(conformer_id)
                #print(conformer_id, "CONVERGED?", converged)
                conformer_energy_list[conformer_id] = AllChem.UFFGetMoleculeForceField(mol,confId=conformer_id).CalcEnergy()

            conformers = mol.GetConformers()
            #print("Energy List", conformer_energy_list)
            #print(f"{len(converged_conformer_id_list)} generated TS Conformers CONVERGED out of {len(conformer_id_list)}")
            #print(f"{converged_conformer_id_list} CONVERGED")

            #for conformer_id in converged_conformer_id_list:
            print ('Conformer XYZs ...')
            for conformer_id in conformer_id_list:
                conformer = conformers[conformer_id]
                #energy = conformer_energy_list[conformer_id]
                #Chem.rdmolfiles.MolToXYZFile(mol, f"{conformer_id}.xyz", confId=conformer_id)
                coordinate_list = []
                for i in range(n): 
                    position = conformer.GetAtomPosition(i) 
                    coordinate_list.append((position[0],position[1],position[2]))
                   # print(virtual_molecule.atom_list[i].get_element(),position[0], position[1], position[2])
                if len(coordinate_list) > 0:
                    coordinate_list = np.array(coordinate_list)
                    coordinates.append(coordinate_list)
        elif library == 'babel': 
            from openbabel import pybel
            from openbabel import openbabel 
            import os
            #### pybel method
            try:
                ob_mol = self.get_ob_mol()
            except:
                virtual_molecule,changed_indices = self.get_valid_molecule()
                ob_mol = virtual_molecule.get_ob_mol()
            pybel_mol = pybel.Molecule(ob_mol)
            for i in range(num_conformer):
                coordinate_list = []
                pybel_mol.make3D()
                pybel_mol.localopt('uff',1000)
                pybel_atom_list = pybel_mol.atoms
                #print ('bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb')
                for atom in pybel_atom_list:
                    position = atom.coords
                    coordinate_list.append((position[0],position[1],position[2]))
                    #print (atom.atomicnum,position[0],position[1],position[2])
                if len(coordinate_list) > 0:
                    coordinate_list = np.array(coordinate_list)
                    coordinates.append(coordinate_list)
        else:
            print ('Give us algorithm for generating 3D!!!')
            print ('You can try your own algorithm here!!!')
            ######### Your own algorithm here #########
            return None
        # (TODO): Need to come up with better algorithm, best one is to simply reoptimize with uff
        # Reupdate molecule
        scale = 1.05 # Make distance to bond length, value between 1.05 ~ 1.2
        for i in range(len(coordinates)):
            coordinate_list = coordinates[i]
            if len(changed_indices) > 0:
                internal_coordinates = self.get_bond_list(False)
                q_updates = dict()
                radius_list = self.get_radius_list()
                for bond in internal_coordinates:
                    start, end = bond
                    if start in changed_indices or end in changed_indices:
                        delta_d = scale * (radius_list[start] + radius_list[end]) - ic.get_distance(coordinate_list,start,end)
                        if np.abs(delta_d) > 0.15:
                            q_updates[bond] = delta_d
                        else:
                            q_updates[bond] = 0.0
                    else:
                        q_updates[bond] = 0.0
                #print ('fixed 3d generation', q_updates)
                updateGood, trj = ic.update_xyz(coordinate_list,q_updates)
                #print("is use_ic_update GOOD?", updateGood)
                #print(trj[-1])
            coordinates[i] = coordinate_list
        return coordinates


    def sample_conformers(self,n_conformer = 20,library = 'rdkit',atom_stereos = dict(),bond_stereos = dict()):
        conformer_list = []
        coordinates = self.make_3d_coordinates(n_conformer,library)
        chg = self.get_chg()
        multiplicity = self.get_multiplicity()
        for coordinate_list in coordinates:
            new_atom_list = [atom.copy() for atom in self.atom_list]
            conformer = Molecule((self.get_z_list(),self.get_adj_matrix(),None,self.get_chg_list()))
            process.locate_molecule(conformer,coordinate_list)
            conformer_list.append(conformer)
        return conformer_list


    def get_coordinate_list(self):
        coordinate_list = [[atom.x,atom.y,atom.z] for atom in self.atom_list]
        return np.array(coordinate_list)
    

    def print_coordinate_list(self,option = 'element'):
        print (self.get_content(option))
        

    def get_content(self,option='element',criteria = 1e-4):
        atom_list = self.atom_list
        content = ''
        for atom in atom_list:
            content += atom.get_content(option,criteria)
        return content.strip()


    def write_geometry(self,file_directory, option='element',criteria = 1e-4):
        """
        Writes xyz file that contains the 3d molecular geometry
    
        :param file_directory(str):
            Directory for saving 3d geometry xyz file

        :return:
        """
        atom_list = self.atom_list
        n = len(atom_list)
        f = open(file_directory, 'w')
        #f = open(file_directory, 'a')
        if True: # If inappropriate geometry condition is determined, it will be added
            content = str(n)+'\n'            
            if self.energy is not None:
                content = content + str(self.energy)+'\n'
            else:
                content = content + '\n'
            f.write(content)
            for atom in atom_list:
                f.write(atom.get_content(option,criteria))
            f.close()
        else:
            print ('Wrong geometry!!!')

    def get_distance_between_atoms(self,idx1,idx2):
        """
        Returns the distance between chosen two atoms

        :param idx1,idx2(int):
            indices of chosen two atoms. 

        :return distance(float):
            Distance between selected two atoms

        """
        coordinate_list = self.get_coordinate_list()
        return ic.get_distance(coordinate_list,idx1,idx2)

    def get_angle_between_atoms(self,idx1,idx2,idx3,unit='rad'):
        """
        Returns the distance between chosen two atoms

        :param idx1,idx2(int):
            indices of chosen two atoms. 

        :return distance(float):
            Distance between selected two atoms

        """
        coordinate_list = self.get_coordinate_list()
        angle = ic.get_angle(coordinate_list,idx1,idx2,idx3)
        if unit == 'degree':
            angle *= 180/np.pi
        return angle


    def get_dihedral_angle_between_atoms(self,idx1,idx2,idx3,idx4,unit='rad'):
        coordinate_list = self.get_coordinate_list()
        angle = ic.get_dihedral_angle(coordinate_list,idx1,idx2,idx3,idx4)
        if unit == 'degree':
            angle *= 180/np.pi
        return angle




class Intermediate(Molecule):
    """
    :class Intermediate:
        class Intermediate mainly contains atom_list, atom_feature, bo_matrix, adj_matrix, energy, smiles
        which are the same thing in class Molecule
        class Intermediate is the set of molecules. Thus, it also has two more variables: molecule_list, a set of molecule and atom_indices_for_each_molecule
        which corresponds to the actual indices within Intermediate for each atom_list in molecule_list
        For example, if intermediate contains two H2 molecules, then one example for atom_indices_for_each_molecule is
        "[[0,1],[2,3]]" or "[[0,2],[1,3]]"
        number_list is **** (need change)

    :param data(various types are possible):
        data that represents Intermediate is used
        Possible input for data: 'CCC.CCC', list of pyclass Molecule

    : param data_type(string):
         the data type for generating class Intermediate. Mainly 'smiles','molecule_list' are possible input
         when 'smiles' is used, it generates molecule using RDkit within class Molecule
         when 'molecule_list' is used, it directly generates Intermediate, by appending molecule in molecule_list
         when None is used, Intermediate with an empty data is generated.
    """
    def __init__(self,data = None):
        self.molecule_list = None
        if type(data) == str:
            super().__init__(data)
        elif type(data) is list:
            try:
                len(data[0])
                super().__init__(data)
            except:
                super().__init__(None)
                self.molecule_list = data
        elif data is None:
            super().__init__(None)            
        else:
            print ('Unknown input type!!!')
            pass
        self.number_list = []
        self.atom_indices_for_each_molecule = None
        if self.molecule_list is None: # Build molecule_list 
            adj_matrix = self.get_adj_matrix() 
            molecule_list = []
            if adj_matrix is not None and data is not None:
                molecule_list = self.get_molecule_list()
                self.molecule_list = molecule_list
        else:
            # First, obtain atom indices
            atom_list = []
            atom_indices_for_each_molecule = []
            cnt = 0
            atom_feature = dict()
            molecule_list = self.molecule_list
            # From molecules obtain information: atom indices (atom_list), atom_feature,
            for molecule in molecule_list:
                atom_indices = []
                molecule_atom_list = molecule.atom_list
                molecule_atom_feature = molecule.atom_feature
                m = len(molecule_atom_list)
                # Generate atom list
                for i in range(m):
                    atom_list.append(molecule_atom_list[i])
                    atom_indices.append(i+cnt)
                # Generate key feature
                if molecule_atom_feature is not None:
                    for key_feature in molecule_atom_feature:
                        if key_feature in atom_feature:
                            atom_feature[key_feature] = np.concatenate((atom_feature[key_feature],molecule_atom_feature[key_feature]))
                        else:
                            atom_feature[key_feature] = molecule_atom_feature[key_feature]
                atom_indices_for_each_molecule.append(atom_indices)
                cnt += m
            total_bo_matrix = np.zeros((cnt,cnt))
            total_adj_matrix = np.zeros((cnt,cnt))
            cnt = 0
            # Now, generate molecule adj/bo matrix using slice
            update_bo_matrix = True
            update_adj_matrix = True
            for molecule in molecule_list:
                m = len(molecule.atom_list)
                bo_matrix = molecule.bo_matrix
                if bo_matrix is None:
                    update_bo_matrix = False
                    # Check adj
                    adj_matrix = molecule.adj_matrix
                    if adj_matrix is None:
                        print ('Molecule connectivity is not given!!!!')
                        update_adj_matrix = False
                    else:
                        total_adj_matrix[cnt:cnt+m,cnt:cnt+m] = adj_matrix[:,:]
                else:
                    total_bo_matrix[cnt:cnt+m,cnt:cnt+m] = bo_matrix[:,:]
                    total_adj_matrix = np.where(total_bo_matrix>0,1,0)
                if molecule.chg is None:
                    molecule.chg = 0
                cnt += m 
            # Give initial values
            self.atom_indices_for_each_molecule = atom_indices_for_each_molecule
            self.atom_list = atom_list
            self.atom_feature = atom_feature
            
            if update_bo_matrix:
                self.bo_matrix = total_bo_matrix
            if update_adj_matrix:
                self.adj_matrix = total_adj_matrix
        
        if self.chg is None:
            if self.molecule_list is not None:
                chg = 0
                for molecule in self.molecule_list:
                    mol_chg = molecule.get_chg()
                    if mol_chg is None:
                        chg = None
                        break
                    chg += mol_chg
                self.chg = chg

    def setCeig(self,Ceig):
        self.Ceig = Ceig

    def set_name(self,name): # Function for ACE-Reaction
        self.name = name

    def get_atom_indices_for_each_molecule(self):
        """
        Computes the atom_index for each atom within molecules. The example of atom_indices_for_each_molecule is provided
        at __init__ section.
        
        :param:
            
        :return atom_indices_for_each_molecule(list of list):
            list of atom_indices for pyclass Intermediate for each molecule (Example shown in __init__)
        """ 
        if self.atom_indices_for_each_molecule is not None:
            return self.atom_indices_for_each_molecule
        index = 0
        n = len(self.atom_list)
        atom_indices_for_each_molecule = []
        if n > 1:
            adj_matrix = self.get_matrix('adj')
            atom_indices_for_each_molecule = process.group_molecules(adj_matrix)
            return atom_indices_for_each_molecule
        else:
            return [[0]]

    def get_molecule_list(self,update=True):
        """
        Returns the molecule_list stored in pyclass Intermediate
        
        :param:
            
        :return molecule_list(list of pyclass 'Molecule'):
            list of class 'Molecule'. If the 'Intermediate' does not contain molecule_list as variable, the method
            finds molecule_list if the atom_list and adjacency_matrix of Intermediate are provided
        """
        if self.molecule_list is not None:
            return self.molecule_list
        else:
            molecule_list = []
            atom_indices_for_each_molecule = self.atom_indices_for_each_molecule
            if atom_indices_for_each_molecule is None:
                atom_indices_for_each_molecule = self.get_atom_indices_for_each_molecule()
            for atom_indices in atom_indices_for_each_molecule:
                molecule = self.get_molecule_from_indices(atom_indices)
                molecule_list.append(molecule)
            if update:
                self.molecule_list = molecule_list
            return molecule_list
    
    def get_molecule_from_indices(self,indices):
        """
        Returns pyclass 'Molecule' for corresponding indices within pyclass 'Intermediate'
        Ex. [1,2,3,4] is given, with atom_list([N atom,H atom, H atom, H atom, O atom, H atom, H atom]), it returns 
        class 'Molecule', corresponding to NH3, 
        
        :param indices(list of int):
        indices of atoms within atom_list of class 'Intermediate' 
        
        :return molecule(pyclass 'Molecule'):
        A pyclass 'Molecule' whose variables that can be copied from class 'Intermediate' are copied
        """
        atom_list = self.atom_list
        adj_matrix = self.adj_matrix
        bo_matrix = self.bo_matrix
        atom_feature = self.atom_feature
        #### Get molecule feature
        n = len(indices)
        molecule_atom_list = []
        molecule_atom_feature = dict()
        molecule = Molecule()
        # Generate atom_list
        for indice in indices:
            molecule_atom_list.append(atom_list[indice])
            molecule.atom_list = molecule_atom_list
        # Generate atom feature
        for key_feature in atom_feature:
            if atom_feature[key_feature] is not None:
                molecule_atom_feature[key_feature] = atom_feature[key_feature][indices]
            else:
                molecule_atom_feature = None
            molecule.atom_feature = molecule_atom_feature
        # Generate bo_matrix, if exists
        if bo_matrix is None:
            # Check adj matrix
            if adj_matrix is None:
                print ('We need connectivity to get information of desired molecule!!!')
            else:
                reduced_adj_matrix = adj_matrix[:,indices]
                reduced_adj_matrix = reduced_adj_matrix[indices,:]
                molecule.adj_matrix = reduced_adj_matrix
        else:
            reduced_bo_matrix = bo_matrix[:,indices]
            reduced_bo_matrix = reduced_bo_matrix[indices,:]
            molecule.bo_matrix = reduced_bo_matrix
        return molecule


    def __eq__(self,intermediate):
        return self.is_same_intermediate(intermediate,True)

    def is_same_intermediate(self,intermediate,option = False):
        """
        Checks whether the two intermediates are same or not. 
        The method checks whether the molecules within molecule_list are the same
                 
        :param intermediate(pyclass 'Intermediate'):
            molecule_list of intermediates are not necessarily needed to be specified. The method automatically 
            checks molecule_list and computes molecule_list if the molecule_list is not provided.

        :return True/False (boolean):
            If True, they are the same. Otherwise, different
        """
        if len(self.atom_list) != len(intermediate.atom_list):
            return False
        molecule_list1 = self.get_molecule_list()
        molecule_list2 = intermediate.get_molecule_list()
        if molecule_list1 == None or molecule_list2 == None:
            print ('intermediate is not well prepared!')
            return False
        elif len(molecule_list1) != len(molecule_list2):
            return False
        n = len(molecule_list1)
        molecule_indices1 = list(range(n))
        molecule_indices2 = list(range(n))
        cnt = n
        while cnt > 0:
            molecule = molecule_list1[molecule_indices1[0]]
            found = False
            for i in range(cnt):
                index = molecule_indices2[i]
                molecule_prime = molecule_list2[index]
                if molecule.is_same_molecule(molecule_prime,option):
                    del(molecule_indices1[0])
                    del(molecule_indices2[i])
                    cnt -= 1
                    found = True
                    break
            if not found:
                break
        return (cnt == 0)

    def copy(self,copy_all = False):
        new_intermediate = Intermediate()
        atom_list = self.atom_list
        # First copy atoms
        new_atom_list = []
        for atom in atom_list:
            new_atom_list.append(atom.copy())
        new_intermediate.atom_list = new_atom_list
        # Copy connectivity information
        bo_matrix = self.get_matrix('bo')
        if bo_matrix is not None:
            new_intermediate.bo_matrix = np.copy(bo_matrix)
        else:
            adj_matrix = self.get_matrix('adj')
            if adj_matrix is not None:
                new_intermediate.adj_matrix = np.copy(adj_matrix)
            else:
                print ('Warning: Connectivity information is not included in the molecule!!!')
        # Finally, copy charge
        chg = 0
        for molecule in self.get_molecule_list():
            mol_chg = molecule.get_chg()
            if mol_chg is not None:
                chg += mol_chg
            else:
                chg = None
                break
        new_intermediate.chg = chg
 
        if 'chg' in self.atom_feature:
            new_intermediate.atom_feature['chg'] = np.copy(self.atom_feature['chg'])
        # Above things are essential for copy
        if copy_all:
            # Copy other attributes
            new_intermediate.energy = self.energy
            new_intermediate.smiles = self.smiles
            new_intermediate.c_eig_list = self.c_eig_list
        return new_intermediate


    def initialize(self):
        atom_indices_for_each_molecule = self.get_atom_indices_for_each_molecule()
        molecule_list = []
        for atom_indices in atom_indices_for_each_molecule:
            molecule = self.get_molecule_from_indices(atom_indices)
            molecule_list.append(molecule)
        self.molecule_list = molecule_list

    def get_smiles(self,method = 'none',find_stereocenter = 'N'):
        """ 
        Returns the smiles of molecules connected with .
        For example, if two methane molecules are given in the intermeidate, 
        it returns 'C.C'
        
        :param:
            
        :return smiles(str):
            Total smiles where each molecule is split by using '.' 
        """
        if self.smiles != None:
            return self.smiles
        else:
            smiles = ''
            bo_matrix = self.get_matrix('bo')
            fc_list = self.get_chg_list()
            bond_list = self.get_bond_list(False)
            if bo_matrix is None:
                adj_matrix = self.get_matrix('adj')
                chg = self.get_chg()
                if adj_matrix is None or chg is None:
                    print ('We need to know both adjacency and charge!!!')
                    return None
                else:
                    bo_matrix = process.get_bo_matrix_from_adj_matrix(self,chg)
                    fc_list = process.get_chg_list_from_bo_matrix(self,chg,bo_matrix)
            elif fc_list is None:
                fc_list = process.get_chg_list_from_bo_matrix(self,chg,bo_matrix)
            molecule_list = self.get_molecule_list()
            atom_indices_for_each_molecule = self.get_atom_indices_for_each_molecule()
            smiles_list = []
            for i in range(len(molecule_list)):
                molecule = molecule_list[i]
                old_chg = molecule.chg
                if old_chg is None:
                    molecule.chg = int(np.sum(fc_list[atom_indices_for_each_molecule[i]]))
                smiles_list.append(molecule.get_smiles(method,find_stereocenter))
                molecule.chg = old_chg
            if smiles_list[0] is None:
                print ('SMILES is not well created !!!')
                smiles = None
            else:
                smiles = '.'.join(smiles_list)
            self.smiles = smiles
            return smiles

    def get_energy(self):
        if self.energy is None:
            energy = 0
            molecule_list = self.get_molecule_list()
            for molecule in molecule_list:
                if molecule.energy is None:
                    return None
                else:
                    energy += molecule.energy
            return energy
        else:
            return self.energy


    def make_3d_coordinates(self,num_conformer = 1,library = 'rdkit'):
        grid = 2.0
        coordinates = super().make_3d_coordinates(num_conformer,library)
        atom_indices_for_each_molecule = self.get_atom_indices_for_each_molecule()
        if atom_indices_for_each_molecule is None:
            return []
        geometry_info = dict()
        atom_list = self.atom_list
        n = len(atom_indices_for_each_molecule)
        translation_vector = dict()
        # Make sphere for each molecule
        for coordinate_list in coordinates:
            for atom_indices in atom_indices_for_each_molecule:
                center = np.mean(coordinate_list[atom_indices],axis=0)
                radius = 0
                for atom_index in atom_indices:
                    r = np.linalg.norm(center - coordinate_list[atom_index]) + atom_list[atom_index].get_radius() * 2
                    if r > radius:
                        radius = r
                geometry_info[tuple(atom_indices)] = [center,radius]
            if len(geometry_info) == 1:
                continue
            # Scatter molecules
            index = 0
            while index < n:
                # Grid search for molecule
                overlap = True
                l = 0
                current_atom_indices = tuple(atom_indices_for_each_molecule[index])
                current_center,current_radius = geometry_info[current_atom_indices]
                good_combination = None
                while overlap:
                    # For every grid, check overlap
                    combinations = []
                    for x in range(l+1):
                        for y in range(l+1-x):
                            combinations.append((x,y,l-x-y))
                    for combination in combinations:
                        x,y,z = combination
                        good_combination = combination
                        vector = np.array([grid*x,grid*y,grid*z])
                        # Check overlap
                        for i in range(index):
                            atom_indices = atom_indices_for_each_molecule[i]
                            center,radius = geometry_info[tuple(atom_indices)]
                            if np.linalg.norm(vector + current_center - center) < radius + current_radius:
                                good_combination = None
                                break
                        if good_combination is not None:
                            overlap = False
                            break
                    if good_combination is not None:
                        break
                    l += 1
                # Move current center
                x,y,z = good_combination
                vector = np.array([grid*x,grid*y,grid*z])
                for atom_index in list(current_atom_indices):
                    coordinate_list[atom_index] += vector
                geometry_info[current_atom_indices][0] += vector
                index += 1
        return coordinates
             

class Conformer:

    def __init__(self,data = None):
        energy = None
        if type(data) == str: # Read directory
            f = open(data,'r')
            atom_list = []
            energy = None
            if data[-4:] == '.xyz':
                try:
                    atom_num = int(f.readline())
                except:
                    print ('Wrong format! Should start with number of atoms!')
                try:
                    energy = float(f.readline())
                    self.energy = energy
                except:
                    a = 1
                for i in range(atom_num):
                    try:
                        content = f.readline().strip()
                        atom_line = content.split()
                        #atomic_number = int(atom_line[0])
                        element_symbol = atom_line[0]
                        x = float(atom_line[1]) 
                        y = float(atom_line[2]) 
                        z = float(atom_line[3])
                        new_atom = Atom(element_symbol)
                        new_atom.x = x
                        new_atom.y = y
                        new_atom.z = z
                        atom_list.append(new_atom)
                    except:
                        print ('Error found in:',content)
                        print ('Check the file again:',data)
                chg = None
                multiplicity = None
            else:
                try:
                    info = f.readline().strip().split()
                    chg = int(info[0])
                    multiplicity = int(info[1])
                except:
                    print ('Wrong format! Should start with number of atoms!')
                while True:
                    try:
                        content = f.readline().strip()
                        atom_line = content.split()
                        #atomic_number = int(atom_line[0])
                        element_symbol = atom_line[0]
                        x = float(atom_line[1]) 
                        y = float(atom_line[2]) 
                        z = float(atom_line[3])
                        new_atom = Atom(element_symbol)
                        new_atom.x = x
                        new_atom.y = y
                        new_atom.z = z
                        atom_list.append(new_atom)
                    except:
                        break
             
        elif data is not None:
            atom_list, chg, multiplicity = data
        self.atom_list = atom_list
        self.chg = chg
        if chg is not None and multiplicity is None:
            z_sum = sum([atom.get_atomic_number() for atom in atom_list])
            multiplicity = (z_sum - chg) % 2
        self.multiplicity = multiplicity
        self.energy = energy

    def __str__(self):
        header = f'{self.chg} {self.multiplicity}\n'
        content = self.get_content()
        return header + content


    def get_coordinate_list(self):
        coordinate_list = [[atom.x,atom.y,atom.z] for atom in self.atom_list]
        return np.array(coordinate_list)
    

    def print_coordinate_list(self,option = 'element'):
        print (self.get_content(option))
        

    def get_content(self,option='element',criteria = 1e-4):
        atom_list = self.atom_list
        content = ''
        for atom in atom_list:
            content += atom.get_content(option,criteria)
        return content.strip()



    def get_adj_matrix(self,coeff=1.2):
        atom_list = self.atom_list
        n = len(atom_list)
        radius_list = [atom.get_radius() for atom in atom_list]
        radius_matrix_flatten = np.repeat(radius_list,n)
        radius_matrix = radius_matrix_flatten.reshape((n,n))
        radius_matrix_transpose = np.transpose(radius_matrix)
        distance_matrix = self.get_distance_matrix()
        criteria_matrix = coeff * (radius_matrix + radius_matrix_transpose)
        adj_matrix = np.where(distance_matrix<criteria_matrix,1,0)
        np.fill_diagonal(adj_matrix,0)
        return adj_matrix

    def get_distance_between_atoms(self,idx1,idx2):
        """
        Returns the distance between chosen two atoms

        :param idx1,idx2(int):
            indices of chosen two atoms. 

        :return distance(float):
            Distance between selected two atoms

        """
        coordinate_list = self.get_coordinate_list()
        return ic.get_distance(coordinate_list,idx1,idx2)

    def get_angle_between_atoms(self,idx1,idx2,idx3,unit='rad'):
        """
        Returns the distance between chosen two atoms

        :param idx1,idx2(int):
            indices of chosen two atoms. 

        :return distance(float):
            Distance between selected two atoms

        """
        coordinate_list = self.get_coordinate_list()
        angle = ic.get_angle(coordinate_list,idx1,idx2,idx3)
        if unit == 'degree':
            angle *= 180/np.pi
        return angle


    def get_dihedral_angle_between_atoms(self,idx1,idx2,idx3,idx4,unit='rad'):
        coordinate_list = self.get_coordinate_list()
        angle = ic.get_dihedral_angle(coordinate_list,idx1,idx2,idx3,idx4)
        if unit == 'degree':
            angle *= 180/np.pi
        return angle



    def get_bo_matrix(self,coeff=1.2):
        pass


    def get_distance_matrix(self):
        coordinate_list = self.get_coordinate_list()
        return spatial.distance_matrix(coordinate_list,coordinate_list)


    def get_ace_mol(self):
        pass


