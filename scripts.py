
#EXECUTE IN TERMINAL, NOT VISUAL
#OTEHRWISE FOLLOW STEPS IN LINE 18

# Getting basic libraries:
import os, sys, math
import numpy as np
from Bio.PDB import *

# Libraries for structure checking: --> 
import biobb_structure_checking
import biobb_structure_checking.constants as cts
from biobb_structure_checking.structure_checking import StructureChecking
base_dir_path=biobb_structure_checking.__path__[0]
args = cts.set_defaults(base_dir_path,{'notebook':True})
###Preparation
    #1 Obtain the required structure from the PDB
f=open("./assignment_data/6m0j.pdb", "r")#for executing in visual use absolute path
parser = PDBParser()
structure = parser.get_structure("6m0j",f)
    #2 Check at PDB which is the composition of a “Biological unit”.
# Remove all chains but those involved in the biological unit, if necessary
base_path = './assignment_data/'
args['output_format'] = "pdb"
args['keep_canonical'] = False
args['input_structure_path'] = base_path + '6m0j.cif'
args['output_structure_path'] = base_path + '6m0j_fixed.pdb'
args['output_structure_path_charges'] = base_path + '6m0j_fixed.pdbqt'
args['time_limit'] = False
args['nocache'] = False
args['copy_input'] = False
args['build_warnings'] = False
args['debug'] = False
args['verbose'] = False
args['coords_only'] = False
args['overwrite'] = False
#Intializing the checking engine, loading the structure and showing info
st_c = StructureChecking(base_dir_path, args)
    #3 Remove all heteroatoms
#Remove Hydrogen
st_c.rem_hydrogen()
#Remove water
st_c.water("yes")
#Remove metals
st_c.metals("All")
#Remove ligands
st_c.ligands("All")
    #4
biobb_structure_checking

####STEP 1


