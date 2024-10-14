# Implementation instruction

## Overview of Implementation

To interface autoCG with other QC packages, you need to implement new python scripts. (We may extend some other QC packages if requested, but it will take some time because we need to get used to those programs) Here, we explain how users can implement a module for interfacing with their available QC package in detail.

Basically, you need to add your own module in the ‘Calculator’ folder for interfacing with a new QC package. In the folder, you will already see three scripts: gaussian.py, orca.py, and template.py. 

In the template.py, you can see that there are three functions that require implementation: **get_energy**, **get_force**, and **relax_geometry** (also, **get_hessian**, if users want to use the original one!): Each function is used for obtaining the energy, the gradient, and the optimized geometry of a given molecular structure. To do so, of course, **three things** should be well implemented: the code should **create the input file** for running a QC package, **run the QC package** with the input file, and **extract the desired result** from the output file. 

Here, we will explain how we implemented gaussian.py, which can perform those three things. Basically, the three functions are almost the same. They just need to slightly change the input creation and read different information depending on their purpose. Therefore, we’ll simply explain the **get_energy** function.

## Example with gaussian.py

First, we switch the current directory to the working directory, where the input file and output file of Gaussian are created. Here, `file_name` and `extra` are not necessary variables! 

```python
def get_energy(self,molecule,chg=None,multiplicity=None,file_name='test',extra=''):
        current_directory = os.getcwd()
        os.chdir(self.working_directory)
```

Then, it writes an input file. Here, we implemented the “make_input”, since this job is repeated. Depending on your coding style, you can implement it differently. 

In the “make_input”, it first writes the content in qc_input file saved in `self.content` and some additional content given as extra. 

```python
def make_input(self,molecule,chg=None,multiplicity=None,option='sp',file_name='test',constraints={},extra=''):
        f = open(f'{file_name}.com','w')
        if chg is None or multiplicity is None:
            chg, multiplicity = self.get_default_mol_params(molecule)
        content = self.get_content()
        content = content + option + f'{extra}\n\ntest\n\n'
        f.write(content)
```

Then, it writes the molecule information (charge, multiplicity, geometry) through write_molecule_info (we strongly recommend using this method). (We omit the explanation of write_molecule_info, since it’s very straightforward)

```python
self.write_molecule_info(f,molecule,chg,multiplicity)
```

Finally, it writes freezing coordinates stored in the `constraints`. `constraints` are dicts where keys are the atom indices of active coordinates and values are the target values.

```python
    self.write_molecule_info(f,molecule,chg,multiplicity)
### Write constraints
        for constraint in constraints:
            constraint_info = []
            if len(constraint) == 2:
                constraint_info.append('B')
            elif len(constraint) == 3:
                constraint_info.append('A')
            else:
                constraint_info.append('D')
            for index in constraint:
                constraint_info.append(str(index+1))
            constraint_info.append('F')
            f.write(' '.join(constraint_info)+'\n')
        f.write('\n') # Additional empty line required
        self.write_basis(f)
        f.write('\n')
        f.close()
```

If input files differ a lot depending on the calculation option, we suggest making input files by directly writing contents with f.write(). 

Back to the get_energy, after writing input, it directly runs Gaussian for the created input, using the `self.command` stored in Gaussian. (Either g09 or g16, depending on your available Gaussian software)

```python
self.make_input(molecule,chg,multiplicity,option='sp',file_name=file_name,extra=extra)
os.system(f'{self.command} {file_name}.com')
```

Then, it reads the output of logfile using cclib, and reads the energy stored in `data.scfenergies` Since, cclib reads the energy as eV unit, we introduce converter to converter, 

```python
p = cclib.parser.Gaussian('{file_name}.log')
data = p.parse()
converter = 1
os.chdir(current_directory)
if self.energy_unit == 'kcal':
        converter = 23.0605506577
if self.energy_unit == 'Hartree':
        converter = 0.036749326681
os.chdir(current_directory) # Swiwtching working_directory to current_directory
return converter*data.scfenergies[-1]
```

Others are almost the same and here, we’ll briefly mention some caution for implementing them.

## Additional notes on Implementation

For force calculation, Gaussian automatically reorients the molecules. Therefore, the obtained force does not correspond to the force applied to the current geometry of the molecules. To remove this problem, we add **“Symmetry=None”** in ****`extra`, so that Gaussian does not try reorientation. If that keyword is written in the input, Gaussian calculates the force of a given geometry without reorientation. 

For relax_geometry, the important part is to introduce constraints within the optimizations. Users should fully understand using QC package that they are trying to run. In case of Gaussian, those options can be given as maxcycles=num_relaxation, maxstep=displacement, and notrust option, written as follows:

```python
if maximal_displacement < 100:
        max_step = int(maximal_displacement*100) + 1
    extra = f'(modredundant,loose,maxcycles={num_relaxation},maxstep={max_step},notrust) Symmetry=None'
```

