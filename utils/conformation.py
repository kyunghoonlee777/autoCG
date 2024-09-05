import numpy as np
import os

from scipy.spatial.transform import Rotation

from autoCG import chem 
from autoCG.utils import process


def sample_from_crest(molecule,constraints=[],working_directory = None,num_process = 1):
    if working_directory == None:
        working_directory = os.getcwd()
    current_directory = os.getcwd()
    os.chdir(working_directory)

    # Make input R.xyz
    molecule.write_geometry('R.xyz')

    # Make constraints.inp file if exist ...
    command = f'crest R.xyz -T {num_process} --noreftopo'
    force = 0.25

    if len(constraints) > 0:
        # Make constraints.inp file ...
        with open('constraints.inp','w') as f:
            f.write('$constrain\n')
            for bond in constraints:
                start, end = bond
                d = molecule.get_distance_between_atoms(start,end)
                f.write(f'  force constant={force}\n')
                f.write(f'  distance: {start+1}, {end+1}, {d}\n')
            f.write('$end') 
        command = command + ' --cinp constraints.inp'
    #print (command)
    # Next, run crest
    os.system(command)
    #exit() 
    # Finally, read conformers
    conformers = process.read_geometries('crest_conformers.xyz')
    #os.system('rm crest_conformers.xyz R.xyz')

    # Move to original directory ...
    os.chdir(current_directory)
    for conformer in conformers:
        conformer.bo_matrix = molecule.bo_matrix
        conformer.chg = molecule.get_chg()
        conformer.multiplicity = molecule.get_multiplicity()
    return conformers


if __name__ == '__main__':
    import sys
    molecule = chem.Molecule(sys.argv[1])
    conformer = molecule.sample_conformers(1)[0]
    try:
        working_directory = sys.argv[2]
    except:
        working_directory = None
    conformers = sample_from_crest(conformer,working_directory)
    for i,conformer in enumerate(conformers):
        print (f'{i+1}th conformer')
        conformer.bo_matrix = molecule.bo_matrix
        conformer.print_coordinate_list()


