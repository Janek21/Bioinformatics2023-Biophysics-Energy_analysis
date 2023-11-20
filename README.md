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

## Processes code description

Up nexts are the explanations for every step broken down of our [code](scripts.py)
#### Preparation steps

In the preparation step the structure of 6m0j is retrieved from PDB by using the ```parser.get_structure()``` imported before as stated in further dependencies.

Then, using ```args``` it is checked at PDB the composition of a biological unit and all chains but the involved in the biological unit are removed.

After this, the next step is the cleansing and checking of the structure, this is done with the ```StructureChecking``` and the ```args```. What we do with this step is remove all heteroatoms, fix the amides, the chirality and the backbone, detect and rebuild if there are missing protein chains and add hydrogens.

After all this, the desired structure is correctly imported, stored and refined.

#### Step 1

#### Step 2

#### Step 3

#### Step 4


