# Energy analysis exercise 2023-24

To be able to execute the code you will need the package biobb_structure_checking that can be installed through \$pip install "biobb_structure_checking\>=3.13.4"

#### [Necessary Modules](./modules_classes.py)

This file contains the modules needed for step 2, these are force fields module, that is used to manage the force field parameters and the residue library module, that is used to manage the residue library.

This both modules are later used in step 2.

We also have an additional directory, called [deprecated code](Deprecated_code) where it is stored all the versions of the functions used in the code before merging it all in the script.py, most of this code won't work as it is not in the directory where it is designed to launch from.

#### [data](./assignment_data)

In this directory is located all the necessary data that is used in the codes, like the 6m0j.pdb or the parameters_vanderw.txt.

#### Further Dependencies

```
Bio.PDB.NeighborSearch (from BioPython)
Bio.PDB.PDBParser (from BioPython)
```
Both this imports are necessary for the correct functioning of the code.
The better course of action would be
```
from Bio.PDB import *
```

In order for it to work properly you have to change the path of the NACCESS to your own.


#### Code execution
Our code can be found in the jupyter notebook named biophysics_project.ipynb
Another option is to execute the python script named scripts.py
