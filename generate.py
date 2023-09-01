### python provided modules ###
import argparse
import datetime
import os
import pickle
import subprocess
import sys
import time
from copy import deepcopy
from itertools import combinations

### extra common libraries ###
import numpy as np
from rdkit import Chem

### ace-reaction libraries ###
from autoCG import chem
from autoCG.Calculator import gaussian, orca
from autoCG.stereochemistry import chiral
from autoCG.utils import arranger, ic, process


class GuessGenerator:
    def __init__(
        self,
        calculator=None,
        ts_scale=1.33,
        form_scale=1.2,
        break_scale=1.8,
        step_size=0.15,
        qc_step_size=0.2,
        num_relaxation=3,
    ):
        self.ts_scale = ts_scale
        self.form_scale = form_scale
        self.break_scale = break_scale
        self.step_size = step_size
        self.qc_step_size = qc_step_size
        self.num_step = 30  # May not change ...
        self.use_gurobi = True
        self.num_conformer = 1
        self.num_relaxation = num_relaxation
        self.unit_converter = 627.5095
        self.energy_criteria = 10.0 / self.unit_converter  # 10.0 kcal/mol
        self.library = "rdkit"
        self.save_directory = None
        # Check calculator
        if calculator is None:  # Default as Gaussian
            #'''
            command = "g09"
            calculator = gaussian.Gaussian(command)  # Run with pm6
            # calculator = gaussian.FastGaussian(command,working_directory) # Run with pm6
            #'''
            # command = 'orca'
            # calculator = orca.Orca(command,working_directory)
        else:
            if not os.path.exists(calculator.working_directory):
                calculator.change_working_directory(
                    os.getcwd()
                )  # Use the final working directory (may be save_directory) as default working directory, it's safe!
        self.calculator = calculator
        self.considerRS = False
        self.considerEndoExo = False
        self.considerEZ = False

    def read_option(self, input_directory):
        if not os.path.exists(input_directory):
            print("No inputs are found!!!")
            return
        fd = open(os.path.join(input_directory, "option"), "r")
        for i in range(1000):
            line = fd.readline()
            if line == "":
                break
            line = line.split()
            if "scale" in line:
                idx = line.index("scale")
                try:
                    self.scale = float(line[idx + 2])
                except:
                    print("input scale is uncertain!!")
            if "b_scale" in line:
                idx = line.index("b_scale")
                try:
                    self.break_scale = float(line[idx + 2])
                except:
                    print("Cannot read b_scale value!!!")
            if "f_scale" in line:
                idx = line.index("f_scale")
                try:
                    self.form_scale = float(line[idx + 2])
                except:
                    print("Cannot read f_scale value!!!")
            if "ts_scale" in line:
                idx = line.index("ts_scale")
                try:
                    self.ts_scale = float(line[idx + 2])
                except:
                    print("Cannot read ts_scale uncertain!!")
            if "guess_step_size" in line:
                idx = line.index("guess_step_size")
                try:
                    self.step_size = float(line[idx + 2])
                except:
                    print("Cannot read step size for generarting guess!!")
            if "guess_num_conformer" in line:
                idx = line.index("guess_num_conformer")
                try:
                    self.num_conformer = int(line[idx + 2])
                except:
                    print("Cannot read the number of ts guess!!")
            if "guess_num_relaxation" in line:
                idx = line.index("guess_num_relaxation")
                try:
                    self.num_relaxation = int(line[idx + 2])
                except:
                    print(
                        "Cannot read the number of optimization for guess generation!!"
                    )
        fd.close()

    def read_qc_input(self, input_directory, file_name="qc_input"):
        qc_directory = os.path.join(input_directory, file_name)
        if os.path.exists(qc_directory):
            self.write_log(f"Read qc_input from {qc_directory}\n\n")
            content = ""
            with open(qc_directory) as f:
                for line in f:
                    content = content + line
            self.calcualtor.content = content
        else:
            self.write_log(
                f"Input template for quantum calculation is not found! Path will be found with default parameters!\n\n"
            )

    def change_save_directory(self, save_directory):
        if not os.path.exists(save_directory):
            os.system(f"mkdir {save_directory}")
            if not os.path.exists(save_directory):
                print("Directory not found! Cannot change save directory!")
            else:
                self.save_directory = save_directory
        else:
            self.save_directory = save_directory

    def change_working_directory(self, working_directory):
        if self.calculator is not None:
            self.calculator.change_working_directory(working_directory)
        else:
            print("Calculator does not exist!!! Define a calcualtor first!")

    def set_energy_criteria(self, criteria):
        self.energy_criteria = criteria / self.unit_converter

    def write_log(self, content, mode="a", file_name="guess"):
        if self.save_directory is None:
            # print ('No log directory!!! Need to give log directory!!!')
            return
        log_directory = os.path.join(self.save_directory, f"{file_name}.log")
        with open(log_directory, mode) as f:
            f.write(content)
            f.flush()

    def get_geometry(self, molecule):
        content = ""
        for atom in molecule.atom_list:
            content = content + atom.get_content()
        return content

    def write_geometry(self, molecule, name="test", extension="xyz"):
        if self.save_directory is None:
            return
        file_name = name + "." + extension
        content = self.get_geometry(molecule)
        with open(os.path.join(self.save_directory, file_name), "w") as f:
            if extension == "xyz":
                f.write(f"{len(molecule.atom_list)}\n\n")
            else:
                try:
                    f.write(f"{molecule.get_chg()} {molecule.get_multiplicity()}\n")
                except:
                    f.write("0 1\n")
            f.write(content + "\n")
            f.flush()

    def write_optimization(self, molecule, file_name=None, mode="a"):
        if self.save_directory is None:
            return
        if file_name is None:
            return
        content = self.get_geometry(molecule)
        e = molecule.get_energy()
        with open(os.path.join(self.save_directory, file_name), mode) as f:
            f.write(f"{len(molecule.atom_list)}\n")
            f.write(f"{e}\n")
            f.write(content)
            f.flush()

    def save_geometries(self, molecules, name="test", extension="xyz"):
        for i, molecule in enumerate(molecules):
            self.save_geometry(molecule, i, name)

    def save_result(self, reaction_info, mapping_info, matching_results):
        if self.save_directory is None:
            return
        with open(os.path.join(self.save_directory, "info.pkl"), "wb") as f:
            pickle.dump(reaction_info, f)
            pickle.dump(mapping_info, f)
            pickle.dump(matching_results, f)

    def get_reaction_info(self, reactant, product):
        reaction_info = dict()
        if self.use_gurobi:
            broken_bonds, formed_bonds = arranger.get_bond_change(reactant, product)
            reaction_info["b"] = broken_bonds
            reaction_info["f"] = formed_bonds
            return reaction_info
        else:  # Can be implemented later ...
            return None

    def make_constraint(self, molecule, bond, scale):
        start = bond[0]
        end = bond[1]
        atom_list = molecule.atom_list
        r1 = atom_list[start].get_radius()
        r2 = atom_list[end].get_radius()
        d = scale * (r1 + r2)
        return d

    def scale_bonds(self, molecule, constraints):
        step_size = self.step_size
        max_num_step = 0
        coordinate_list = np.array(molecule.get_coordinate_list())
        update_q = dict()
        displacement = 0
        for constraint in constraints:
            start = constraint[0]
            end = constraint[1]
            delta_d = constraints[constraint] - ic.get_distance(
                coordinate_list, start, end
            )
            update_q[constraint] = delta_d
            displacement += abs(delta_d)
            if abs(delta_d) / step_size > max_num_step - 1:
                max_num_step = int(abs(delta_d) / step_size) + 1
        # print (update_q)
        n = len(coordinate_list)
        for constraint in update_q:
            update_q[constraint] /= max_num_step
        displacement /= max_num_step
        # Now use internal coordinate
        for i in range(max_num_step):
            converged = ic.update_geometry(molecule, update_q)
            if not converged:
                print("scaling update not converged!!!")
            # Relax geometry

    def find_spectator_molecules(self, reactant, reaction_info):
        # Find participating molecules (Non-participating is removed ...)
        reactant_atom_indices_for_each_molecule = (
            reactant.get_atom_indices_for_each_molecule()
        )
        molecule_indices = []
        # print (reactant_atom_indices_for_each_molecule)
        for bond in reaction_info["f"] + reaction_info["b"]:
            for i, atom_indices in enumerate(reactant_atom_indices_for_each_molecule):
                if bond[0] in atom_indices or bond[1] in atom_indices:
                    if i not in molecule_indices:
                        molecule_indices.append(i)
        # print (molecule_indices)
        molecule_list = reactant.get_molecule_list()
        not_used_molecule_indices = list(
            set(range(len(reactant_atom_indices_for_each_molecule)))
            - set(molecule_indices)
        )
        spectating_molecules = [molecule_list[i] for i in not_used_molecule_indices]

        used_atom_indices = []
        for i in molecule_indices:
            used_atom_indices = (
                used_atom_indices + reactant_atom_indices_for_each_molecule[i]
            )
        used_atom_indices.sort()

        # Make mapping that reduces the total_atom_indices to used_atom_indices
        mapping_info = dict()
        for i, index in enumerate(used_atom_indices):
            mapping_info[index] = i
        for i in range(len(reaction_info["f"])):
            start, end = reaction_info["f"][i]
            reaction_info["f"][i] = (mapping_info[start], mapping_info[end])
        for i in range(len(reaction_info["b"])):
            start, end = reaction_info["b"][i]
            reaction_info["b"][i] = (mapping_info[start], mapping_info[end])

        return reaction_info, mapping_info, spectating_molecules

    def get_ts_like_molecules(self, reactant, reaction_info, n_conformer=1):
        # Get reduced molecule and modify reaction_info with reduced indices
        z_list = np.copy(reactant.get_z_list())
        adj_matrix = np.copy(reactant.get_matrix("adj"))
        for bond in reaction_info["f"]:
            start, end = bond
            adj_matrix[start][end] += 1
            adj_matrix[end][start] += 1
        ts_molecule = chem.Molecule([z_list, adj_matrix, None, None])
        # Modify bo matrix and make hypothetical ts molecule, also sample geometries
        # print ('aklsdjfalkdjfas',reduced_chg_list)

        # Reset charge if possible ...
        ts_molecule.chg = reactant.get_chg()
        # print ('chgchgchg',ts_molecule.chg,ts_molecule.get_chg_list())
        ts_molecules = ts_molecule.sample_conformers(
            self.num_conformer, library=self.library
        )

        print("hestetsetset")
        if self.considerRS or self.considerEndoExo:
            reactive_atoms = np.unique(np.array(reaction_info["f"]).flatten())
        
            print("reactive_atoms: ", reactive_atoms)
            self.write_log("###### Sample Stereoisomers ######\n")
            self.write_log(f"Stereochemistry including following atoms will be considered: {reactive_atoms}\n")
            

            if self.considerRS:
                container = []
                self.write_log("Chirality option is ACTIVE\n")
                for ts_conformer in ts_molecules:
                    container += chiral.sampleRS(ts_conformer)
                ts_molecules += container

            if self.considerEndoExo:
                container = []
                self.write_log("Endo/Exo option is ACTIVE\n")
                for ts_conformer in ts_molecules:
                    stereo = deepcopy(ts_conformer)
                    for bond in reaction_info["f"]:
                        start, end = bond
                        stereo.adj_matrix[start][end] -= 1
                        stereo.adj_matrix[end][start] -= 1
                    chiral.changeEndoExo(stereo, reactive_atoms)
                    for bond in reaction_info["f"]:
                        start, end = bond
                        stereo.adj_matrix[start][end] += 1
                        stereo.adj_matrix[end][start] += 1

                    container.append(ts_conformer)
                    container.append(stereo)
                ts_molecules += container
            
            if self.considerEZ:
                container = []
                self.write_log("Double Bond E/Z option is ACTIVE\n")
                for ts_conformer in ts_molecules:
                    container += chiral.sampleEZ(ts_conformer)
                ts_molecules += container


        return ts_molecules  # constraint_info for ts_molecule, mapping_info: intermediate to ts_molecule (reduced_idx)
    
        """
        def sample_diastereomers(self, 
                                ts_conformer: chem.Molecule, 
                                reactive_atoms):
            
            if self.considerRS:
                return chiral.sampleRS(ts_conformer)
                #return chiral.sampleRS(ts_conformer, scope=reactive_atoms)
            
            if self.considerEndoExo:
                stereo = deepcopy(ts_conformer)
                chiral.changeEndoExo(stereo, tuple(reactive_atoms))
                return [ts_conformer, stereo]
        """
        """
        z_list = ts_conformer.get_z_list()
        import queue
        q = queue.SimpleQueue()
        def prior(adj_matrix:np.ndarray, queue:queue.SimpleQueue, depth:int):
            center = queue.get()

            if depth > 2:
                return center

            neighbors = np.where(adj_matrix[center] > 1, 
                                 adj_matrix[center], 
                                 -1)

            if np.all(neighbors) == 0:
                return center
            
            p = z_list[neighbors]
            if not np.all(p == p[0]):
                return np.argmax(p)
            else:
                for neighbor in neighbors:
                    queue.put(x)
                return prior(adj_matrix, queue, depth+1)
        """


    def move_to_reactant(self, molecule, formed_bonds, broken_bonds, file_name=None):
        step_size = self.step_size
        coordinate_list = molecule.get_coordinate_list()
        n = len(molecule.atom_list)
        energy_list = [molecule.energy]
        trajectory = []
        # total_bonds = broken_bonds + formed_bonds # Must fix broken bonds
        total_bonds = molecule.get_bond_list(False)
        # print ('total bonds: ',total_bonds)
        update_q = dict()
        check_increase = False
        undesired_bonds = []
        self.write_optimization(molecule, file_name, "w")
        for i in range(self.num_step):
            update_q = dict()
            for bond in total_bonds:
                if bond in formed_bonds:
                    start, end = bond
                    d = molecule.get_distance_between_atoms(start, end)
                    r1 = molecule.atom_list[start].get_radius()
                    r2 = molecule.atom_list[end].get_radius()
                    delta = d - (r1 + r2) * self.form_scale
                    if delta > 0:
                        update_q[bond] = -min(step_size, delta)
                    else:
                        update_q[bond] = 0.0
                elif bond in broken_bonds:
                    start, end = bond
                    d = molecule.get_distance_between_atoms(start, end)
                    r1 = molecule.atom_list[start].get_radius()
                    r2 = molecule.atom_list[end].get_radius()
                    delta = (r1 + r2) * self.break_scale - d
                    if delta > 0:
                        update_q[bond] = min(step_size, delta)
                    else:
                        update_q[bond] = 0.0
                else:
                    update_q[bond] = 0.0
            update = sum([abs(value) for value in list(update_q.values())])
            if update < 0.0001:  # If nothing to update ...
                print("No more coordinate to update ...")
                break
            # print (f'{i+1}th geometry ...')
            # molecule.print_coordinate_list()
            converged = ic.update_geometry(molecule, update_q)
            if n * (n - 1) / 2 > len(update_q):
                self.calculator.relax_geometry(
                    molecule,
                    update_q,
                    None,
                    None,
                    "test",
                    self.num_relaxation,
                    self.step_size,
                )
            self.write_optimization(molecule, file_name, "a")

    def relax_geometry(self, reactant, constraints, file_name=None):
        energy_list = [self.calculator.get_energy(reactant)]
        self.write_optimization(reactant, file_name, "w")
        for i in range(self.num_step):
            self.calculator.relax_geometry(
                reactant,
                constraints,
                None,
                None,
                "test",
                self.num_relaxation,
                self.qc_step_size,
            )
            self.write_optimization(reactant, file_name, "a")
            energy_list.append(reactant.energy)
            if energy_list[-1] - energy_list[-2] > -self.energy_criteria:
                break
        # print ('#### Optimized geometry ####')
        # reactant.print_coordinate_list()

    def get_oriented_RPs(
        self,
        reactant,
        product,
        reaction_info=None,
        chg=None,
        multiplicity=None,
        save_directory=None,
    ):
        working_directory = None
        # Check save_directory
        if type(save_directory) is str:
            if not os.path.exists(save_directory):
                os.system(f"mkdir {save_directory}")
                if not os.path.exists(save_directory):
                    print("Path does not exist! Guess is generated without log file!")
                    save_directory = None
                else:
                    working_directory = save_directory
            else:
                working_directory = save_directory
        if working_directory is not None:
            if self.calculator.working_directory == os.getcwd():
                self.calculator.change_working_directory(working_directory)
        self.save_directory = save_directory

        starttime = datetime.datetime.now()
        self.write_log(f"##### Guess generator info #####\n", "w")
        # Parameter directly related to precomplexes
        self.write_log(f"library: {self.library}\n")
        self.write_log(f"ts_scale: {self.ts_scale}\n")
        self.write_log(f"form_scale: {self.form_scale}\n")
        self.write_log(f"broken_scale: {self.break_scale}\n")
        self.write_log(f"num_conformer: {self.num_conformer}\n")
        self.write_log(
            f"E_criteria(Delta E) = {self.energy_criteria * self.unit_converter}kcal/mol\n"
        )

        # Fine tuning parameters for delicate generation
        self.write_log(
            f"step_size (change in coordinates for bond participating coordinates): {self.step_size}\n"
        )
        self.write_log(f"Maximal displacement change: {self.qc_step_size}\n")
        self.write_log(f"num_relaxation: {self.num_relaxation}\n")
        self.write_log(f"Calculator: {self.calculator.command}\n\n")

        self.write_log(f"Starting time: {starttime}\n\n")
        self.write_log(f"All guesses will be saved in {self.save_directory}\n")

        if reaction_info is None:
            self.write_log("Finding reaction information with gurobi!\n")
            reaction_info = self.get_reaction_info(reactant, product)
            if reaction_info is None:
                self.write_log(
                    "Internal error occured for finding reaction information! Guess generation terminated ..."
                )
                exit()
        else:
            self.write_log("Using the provided reaction!\n")

        reactant_copy = reactant.copy()
        if chg is not None:
            reactant_copy.chg = chg
        if multiplicity is not None:
            reactant_copy.multiplicity = multiplicity
        # Need to write reaction information here!
        self.write_log(f"####### Final reaction information #######\n")
        bond_form = reaction_info["f"]
        bond_break = reaction_info["b"]
        print("reaction info:", reaction_info)
        try:
            self.write_log(
                f'reaction SMILES: {reactant.get_smiles("ace")}>>{product.get_smiles("ace")}\n'
            )
        except:
            self.write_log("reaction SMILES cannot be generated !!!\n")

        # self.write_log(f'reaction SMILES: {reactant.get_smiles("ace")}>>{product.get_smiles("ace")}\n')
        self.write_log(f"bond form: {bond_form}\n")
        self.write_log(f"bond dissociation: {bond_break}\n\n")
        # Mapped reaction_info (reduced), to reconstruct original, use mapping_info!
        self.write_log("Start to get ts-like molecules.\n")
        st = datetime.datetime.now()
        (
            reaction_info,
            mapping_info,
            spectating_molecules,
        ) = self.find_spectator_molecules(reactant, reaction_info)
        if len(spectating_molecules) == 0:
            reduced_reactant = reactant.copy()
            if product is not None:
                reduced_product = product.copy()
            else:
                reduced_product = None
        else:
            reduced_reactant = process.get_reduced_intermediate(reactant, mapping_info)
            if product is not None:
                reduced_product = process.get_reduced_intermediate(
                    product, mapping_info
                )
            else:
                reduced_product = None

        # Check total charge of reduced_reactant
        ts_molecules = self.get_ts_like_molecules(reduced_reactant, reaction_info)
        #ts_molecules = ts_molecules[1:] #### MUST BE DELETED!!!!
        # exit()
        et = datetime.datetime.now()
        self.write_log(
            f"[{et}] {len(ts_molecules)} ts conformers generated... Taken time: {et-st}\n\n"
        )
        # Note that ts_molecules only contain atoms of molecules that are participating in the reaction
        # Relax ts_molecules with constraints
        energy_list = []
        ts_constraints = dict()
        if len(ts_molecules) == 0:
            print("No conformer generated !!!")
            return [], reaction_info, mapping_info, []

        molecule = ts_molecules[0]
        # Build constraints for ts optimization
        atom_list = ts_molecules[0].atom_list
        n = len(ts_molecules)
        bond_list = molecule.get_bond_list(False)
        for bond in bond_list:
            if bond in reaction_info["f"]:
                ts_constraints[bond] = self.make_constraint(
                    molecule, bond, self.ts_scale
                )
            elif bond in reaction_info["b"]:
                ts_constraints[bond] = self.make_constraint(
                    molecule, bond, self.ts_scale
                )
            else:
                ts_constraints[bond] = self.make_constraint(molecule, bond, 1.0)
        self.write_log("Set up done!!!\n")
        print("n:", n)
        if chg is None:
            chg = molecule.chg
            if chg is None:
                print("Charge is not given !!!")
        if multiplicity is None:
            multiplicity = (np.sum(reduced_reactant.get_z_list()) - chg) % 2 + 1
        reactant_copy = reduced_reactant.copy()
        if chg is not None:
            reactant_copy.chg = chg
        if multiplicity is not None:
            reactant_copy.multiplicity = multiplicity

        # current_directory = os.getcwd()
        # self.write_log('Relaxing TS structures ...\n\n')
        for i in range(n):
            # Write guess log for each directory
            print("final relaxation ...")
            # print (ts_constraints)
            ts_molecule = ts_molecules[i]
            self.scale_bonds(ts_molecule, ts_constraints)
            # self.save_geometry(ts_molecules[i],i,'relaxed_ts','xyz')
            st = datetime.datetime.now()
            save_directory = self.save_directory
            # Log changed into the specific folder directory
            save_directory = self.save_directory
            if save_directory is not None:
                folder_directory = os.path.join(save_directory, str(i + 1))
                name = f"guess_{i+1}"
                os.system(f"mkdir -p {folder_directory}")
                self.save_directory = folder_directory

            folder_directory = os.path.join(save_directory, str(i + 1))
            name = f"guess_{i+1}"
            os.system(f"mkdir -p {folder_directory}")
            self.write_log(
                f"Start geometry optimization for {i+1}th guess of TS structure.\n"
            )
            # Log changed into the specific folder directory
            self.save_directory = folder_directory
            self.write_log(
                f"{i+1}th Initial guess generation\n\nInitial\n{self.get_geometry(ts_molecules[i])}",
                file_name=name,
                mode="w",
            )
            # Repeat relaxation and separating bonds
            undesired_bonds = []
            formed_bonds = reaction_info["f"]
            broken_bonds = reaction_info["b"]
            self.write_geometry(ts_molecule, "initial_ts", "xyz")
            self.write_log(
                f"\nFinal\n{self.get_geometry(ts_molecule)}\n", file_name=name
            )
            et = datetime.datetime.now()
            # Move to original guess geometry
            self.save_directory = save_directory
            energy_list.append(ts_molecules[i].energy)
            self.calculator.clean_scratch()
        # os.chdir(current_directory)
        self.write_log("TS generation finished!!! \n\n")
        #print("Energy: ", energy_list)
        # May filter some ts_molecules using the energy information
        reactant_product_pairs = []
        matching_results = []
        n = len(ts_molecules)
        for i in range(n):
            self.write_log(f"Finding RP pair for {i+1}th TS structure.\n")
            # Move save_directory
            save_directory = self.save_directory
            if save_directory is not None:
                folder_directory = os.path.join(save_directory, str(i + 1))
                self.save_directory = folder_directory
            name = f"guess_{i+1}"
            ts_molecule = ts_molecules[i]
            # Set new constraints from the reaction information
            reactant_molecules = ts_molecule.copy()
            reactant_molecules.energy = ts_molecule.energy
            product_molecules = ts_molecule.copy()
            product_molecules.energy = ts_molecule.energy
            # print ('bf',reactant_molecules.chg,product_molecules.chg,ts_molecule.chg)
            if chg is not None:
                reactant_molecules.chg = chg
                product_molecules.chg = chg
            if multiplicity is not None:
                reactant_molecules.multiplicity = multiplicity
                product_molecules.multiplicity = multiplicity
            self.write_log(f"{i+1}th structure moving to R ...\n", file_name=name)
            # print ('moving to R ...')
            self.move_to_reactant(
                reactant_molecules, broken_bonds, formed_bonds, "TS_to_R.xyz"
            )  # New module for finding reactant complex
            undesired_bonds = []
            # print ('R geometry ...')
            # reactant_molecules.print_coordinate_list()
            # Modify connectivity
            for bond in formed_bonds:
                start, end = bond
                reactant_molecules.adj_matrix[start][end] -= 1
                reactant_molecules.adj_matrix[end][start] -= 1
            print("###### R optimization #####")
            self.relax_geometry(reactant_molecules, [], "opt_R.xyz")
            print("R optimization finished !!!")
            # reactant_molecules.print_coordinate_list()
            self.write_geometry(reactant_molecules, "R", "xyz")
            self.write_geometry(reactant_molecules, "R", "com")
            self.write_log(f"{i+1}th structure moving to P ...\n", file_name=name)
            self.calculator.clean_scratch()
            # print ('moving to P ...')
            # print (formed_bonds,broken_bonds)
            self.move_to_reactant(
                product_molecules, formed_bonds, broken_bonds, "TS_to_P.xyz"
            )
            # self.calculator.optimize_geometry(molecule=product_molecules,constraints=reaction_info['b'],chg=chg,multiplicity=multiplicity,file_name='product',extra='(modredundant,maxcycles=10) Symmetry=None')
            #'''
            undesired_bonds = []
            # print ('P geometry ...')
            # product_molecules.print_coordinate_list()
            for bond in broken_bonds:
                start, end = bond
                product_molecules.adj_matrix[start][end] -= 1
                product_molecules.adj_matrix[end][start] -= 1
            print("###### P optimization #####")
            # self.separated_reactant_optimization(product_molecules,
            self.relax_geometry(product_molecules, [], "opt_P.xyz")
            print("P optimization finished !!!")
            self.write_geometry(product_molecules, "P", "xyz")
            self.write_geometry(product_molecules, "P", "com")

            reactant_molecules.adj_matrix = None
            reactant_molecules.bo_matrix = None
            product_molecules.adj_matrix = None
            product_molecules.bo_matrix = None

            RP_pair = (reactant_molecules, product_molecules)
            reactant_product_pairs.append(RP_pair)
            self.calculator.clean_scratch()

            # Back to original ...
            self.save_directory = save_directory
            self.write_log(
                f"[{datetime.datetime.now()}] {i+1}th reactant-product pair well created !!!\n"
            )
            matching_result = [True, True]
            # print ('Checking R connectivity ... ')
            # print (reduced_reactant.get_chg())
            # print (reduced_reactant.get_adj_matrix())
            matching_result[0], coeff = process.is_same_connectivity(
                reduced_reactant, reactant_molecules
            )
            # print (reduced_reactant.get_smiles('ace'))
            if reduced_product is not None:
                # print ('Checking P connectivity ...')
                matching_result[1], coeff = process.is_same_connectivity(
                    reduced_product, product_molecules
                )
            if matching_result[0]:
                # print ('Well matching R generated !!!')
                content = f"{i+1}th reactant well matches with the original one!\n"
            else:
                # print ('Well matching R not found !!!')
                content = f"{i+1}th reactant does not match with the original one!\n"
            self.write_log(content)
            if matching_result[1]:
                # print ('Well matching P generated !!!')
                content = f"{i+1}th product well matches with the original one!\n\n"
            else:
                # print ('Well matching P not found !!!')
                content = f"{i+1}th product does not match with the original one!\n\n"
            self.write_log(content)
            matching_results.append(matching_result)
            # print (matching_result,matching_results)
            print("af", reactant_molecules.chg, product_molecules.chg)
        # Save info in the save directory
        print(matching_results)
        self.save_result(reaction_info, mapping_info, matching_results)

        # All molecules here, are in ts connectivity matrix!!! To reset, need to remove bonds using reaction_info
        endtime = datetime.datetime.now()
        self.write_log(
            f"[{endtime}] Total {endtime-starttime} was taken for making {n} RP pairs.\n"
        )
        self.save_directory = None

        return (
            reactant_product_pairs,
            reaction_info,
            mapping_info,
            matching_results,
        )  # Reduced indices for reaction_info

    def get_oriented_RPs_from_smiles(
        self, reaction_smiles, chg=None, multiplicity=None, save_directory=None
    ):
        smiles_list = reaction_smiles.split(">>")
        reactant_smiles = smiles_list[0]
        product_smiles = smiles_list[1]
        reactant = chem.Intermediate(reactant_smiles)
        product = chem.Intermediate(product_smiles)
        return self.get_oriented_RPs(
            reactant,
            product,
            None,
            chg=chg,
            multiplicity=multiplicity,
            save_directory=save_directory,
        )

    def get_oriented_RPs_from_geometry(
        self, input_directory, chg=None, multiplicity=None, save_directory=None
    ):
        # Input directory: name.com, coordinates
        reaction_info = dict()
        reaction_info["f"] = []
        reaction_info["b"] = []
        atom_list = []
        reactant = chem.Intermediate()
        with open(input_directory) as f:
            # Read chg, multiplicity
            chg, multiplicity = f.readline().strip().split()
            chg = int(chg)
            multiplicity = int(multiplicity)
            # Read molecule geometry
            content = f.readline()
            while content != "\n":
                element, x, y, z = content.split()
                x = float(x)
                y = float(y)
                z = float(z)
                new_atom = chem.Atom(element)
                new_atom.x = x
                new_atom.y = y
                new_atom.z = z
                atom_list.append(new_atom)
                content = f.readline()
            # Read coordinates
            while True:
                line = f.readline()
                try:
                    coordinate_info = line.strip().split(" ")
                    start = int(coordinate_info[0]) - 1
                    end = int(coordinate_info[1]) - 1
                    reaction_type = coordinate_info[2].lower()  # b or f
                    bond = (start, end)
                    if start > end:
                        bond = (end, start)
                    reaction_info[reaction_type].append(bond)
                except:
                    break
        reactant.atom_list = atom_list
        reactant.adj_matrix = process.get_adj_matrix_from_distance(reactant, 1.2)
        return self.get_oriented_RPs(
            reactant,
            None,
            reaction_info,
            chg=chg,
            multiplicity=multiplicity,
            save_directory=save_directory,
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", "-i", type=str, help="Input directory for reaction")
    parser.add_argument(
        "--save_directory",
        "-sd",
        type=str,
        help="save directory for precomplex",
        default=None,
    )
    parser.add_argument(
        "--library",
        "-l",
        type=str,
        help="Library for generating conformers of rough TS guess",
        default="rdkit",
    )
    parser.add_argument(
        "--calculator",
        "-ca",
        type=str,
        help="Using calculator for constraint optimization",
        default="orca",
    )
    parser.add_argument(
        "--num_conformer",
        "-nc",
        type=int,
        help="Number of precomplex conformations",
        default=1,
    )
    parser.add_argument(
        "--num_relaxation",
        "-nr",
        type=int,
        help="Number of relaxation during each scan",
        default=3,
    )
    parser.add_argument(
        "--ts_scale",
        "-ts",
        type=float,
        help="Scaling factor for bonds participating reactions",
        default=1.33,
    )
    parser.add_argument(
        "--form_scale",
        "-fs",
        type=float,
        help="Scaling factor for bond formation",
        default=1.2,
    )
    parser.add_argument(
        "--break_scale",
        "-bs",
        type=float,
        help="Scaling factor for bond dissociation",
        default=1.8,
    )
    parser.add_argument(
        "--step_size",
        "-s",
        type=float,
        help="Step size for constraint optimization",
        default=0.15,
    )
    parser.add_argument(
        "--qc_step_size",
        "-qc",
        type=float,
        help="Maximal displacement change during QC calculation",
        default=0.25,
    )
    parser.add_argument(
        "--energy_criteria",
        "-ec",
        type=float,
        help="Energy criteria for precomplex optimization",
        default=10.0,
    )
    parser.add_argument("--chg", type=int, help="Total charge of the system", default=0)
    parser.add_argument(
        "--mult", type=int, help="Total multiplicity of the system", default=1
    )
    parser.add_argument("--RS",
                        "-rs",
                        action="store_true",
                        help="sample R/S conformers when reaction center is chiral center"
    )
    parser.add_argument("--EndoExo",
                        "-nx",
                        action="store_true",
                        help="sample Endo/Exo conformers when ring is formed during the reaction"
    )

    args = parser.parse_args()
    save_directory = args.save_directory
    print(args)
    if save_directory is None:
        input_directory = args.input
        names = input_directory.split("/")
        save_directory = "/".join(names[:-1])
    if args.calculator == "gaussian":
        calculator = gaussian.Gaussian()
    else:
        calculator = orca.Orca()
    generator = GuessGenerator(
        calculator,
        args.ts_scale,
        args.form_scale,
        args.break_scale,
        args.step_size,
        args.qc_step_size,
        args.num_relaxation,
    )
    generator.library = args.library
    generator.num_conformer = args.num_conformer
    generator.set_energy_criteria(args.energy_criteria)

    generator.considerRS = args.RS
    generator.considerEndoExo = args.EndoExo

    if ">>" in args.input:
        (
            RP_pairs,
            reaction_info,
            mapping_info,
            matching_results,
        ) = generator.get_oriented_RPs_from_smiles(
            args.input, save_directory=save_directory
        )
    else:
        (
            RP_pairs,
            reaction_info,
            mapping_info,
            matching_results,
        ) = generator.get_oriented_RPs_from_geometry(
            args.input, save_directory=save_directory
        )
