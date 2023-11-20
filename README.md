# Energy analysis exercise 2023-24

To be able to execute the code you will need the package biobb_structure_checking that can be installed through \$pip install "biobb_structure_checking\>=3.13.4"

#### [Necessary Modules](./modules_classes.py)

This file contains the modules needed for step 2, these are force fields module, that is used to manage the force field parameters and the residue library module, that is used to manage the residue library.

This both modules are later used in step 2.

#### [data](./assignment_data)

In this directory is located all the necessary data that is used in the codes, like the 6m0j.pdb or the parameters_vanderw.txt.

#### Further Dependencies

```
Bio.PDB.NeighborSearch (from BioPython)
Bio.PDB.PDBParser (from BioPython)
```
Both this imports are necessary for the correct functioning of the code.


####[code](scrips.py) execution

In this file there are all the steps together so the execution of this file should be enough
