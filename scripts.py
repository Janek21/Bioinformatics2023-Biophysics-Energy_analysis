
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
#Fix amides
st_c.amide("All")
#fix the chirality of some aa
st_c.chiral("All") 
#fix the backbone
st_c.backbone('--fix_atoms All --fix_chain none --add_caps none')
#detects and rebuilds missing protein side chains
st_c.fixside("All")
#add hydrogens
st_c.add_hydrogen("auto")

#you could also do everything with st_c.checkall() but we did it manually so its clearer what we do
st_c._save_structure(args['output_structure_path'])

#st_c.rem_hydrogen('yes')
#st_c.add_hydrogen('--add_charges --add_mode auto')
#st_c._save_structure(args['output_structure_path_charges'])


####STEP 1
from step1_2 import *
dt=5 #quin es el distance treshold?
chain_atoms(structure, dt)


####STEP 2

# Importing libraries:
import argparse
import sys
import os
import math

from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.NACCESS import NACCESS_atomic
from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.PDBIO import PDBIO, Select

# Classes / functions needed to calculate VanderWaals parameters or the residue library
#create file with classes for easy loading
#L ARXIU INCLU LES LIBRARIES; PODEM TREURE LO DE SOBRE

from modules_classes import *

# loading residue library from data/aaLib.lib
res_lib = ResiduesDataLib('./assignment_data/parameters_vanderw.txt')
# loading VdW parameters
ff_params = VdwParamset('./assignment_data/parameters_step2.txt')

#PQ NO VA?

# set the pdb_path and load the structure
pdb_path = "./assignment_data/6m0j_fixed.pdb"
# Setting the Bio.PDB.Parser object
parser = PDBParser(PERMISSIVE=1)
# Loading structure
st = parser.get_structure('st', pdb_path)



#Add atom parameters

def add_atom_parameters(st, res_lib, ff_params):
    for at in st.get_atoms():
        resname = at.get_parent().get_resname() #Finds residues 
        params = res_lib.get_params(resname, at.id) #Gets the parameters of the atoms
        if not params:  # If it's not in an AA get the atom without charge
            print(at)
            at.xtra['atom_type'] = at.element
            at.xtra['charge'] = 0
        else:           #If it's in an AA assignt type and charge to the atom
            at.xtra['atom_type'] = params.at_type
            at.xtra['charge'] = params.charge
        at.xtra['vdw'] = ff_params.at_types[at.xtra['atom_type']]


add_atom_parameters(st, res_lib, ff_params)





f.close()
