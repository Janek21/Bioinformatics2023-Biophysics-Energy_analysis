
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



####STEP 2

# Importing libraryes:
import argparse
import sys
import os
import math

from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.NACCESS import NACCESS_atomic
from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.PDBIO import PDBIO, Select

# Classes / functions needed to calculate VanderWaals parameters or the residue library

class ResiduesDataLib(): # Creating the class

    def __init__(self, fname): # initialitation of the class with closed variables
        self.residue_data = {}

        try:
            fh = open(fname, "r") # Oppening the file for reading

        except OSError: # Searching for errors and if found fiving msg error and quitting
            print("#ERROR while loading library file (", fname, ")")
            sys.exit(2)

        for line in fh: # Reading the file line per line

            if line[0] == '#': # Checking for the header and continuing
                continue

            data = line.split() # splitting the data
            r = Residue(data) # calling the class Residue to get the residue

            self.residue_data[r.id] = r # Assigning the residue to a local var of the class

        self.nres = len(self.residue_data) # Getting the lenght of the variable and assigining it to a local variable

    def get_params(self, resid, atid): # Defining function to get the parameters
        atom_id = resid + ':' + atid # getting the atom id

        if atom_id in self.residue_data: # Checking if the atom id is inside the residue data and if so returning it
            return self.residue_data[atom_id]
        
        else: # If not in the data getting an error
            print("WARNING: atom not found in library (", atom_id, ')')
            return None

class Residue(): # Defining the classs
    
    def __init__(self,data): # Initializing the class
        self.id     = data[0]+':'+data[1] # getting id
        self.at_type = data[2] # 
        self.charge  = float(data[3]) # gettig charge
        
class AtomType(): # Defining class

    def __init__(self, data): # Initializing the class
        # Getting some info and storing it in local variables
        self.id   = data[0]
        self.eps  = float(data[1])
        self.sig  = float(data[2])
        self.mass = float(data[3])
        self.fsrf = float(data[4])
        self.rvdw = self.sig * 0.5612
        
class VdwParamset(): #extracted from GELPI's github
    #parameters for the VdW
    def __init__ (self, file_name): # initializing the class
        self.at_types = {}

        try: # oppening the file for reading
            fh = open(file_name, "r")

        except OSError: # checking for errors and quitting
            print ("#ERROR while loading parameter file (", file_name, ")")
            sys.exit(2)

        for line in fh: # Reading the file

            if line[0] == '#': # checking for the headder
                continue

            data = line.split() # splitting the line
            self.at_types[data[0]] = AtomType(data) # getting the atom type by calling the class AtomType

        self.ntypes = len(self.at_types) # assigning the length of the variable to a local var
        fh.close() # clossing




































#Add atom parameters

def add_atom_parameters(st, res_lib, ff_params):
    for at in st.get_atoms():
        resname = at.get_parent().get_resname() #Finds residues 
        params = res_lib.get_params(resname, at.id)
        if not params:
        params = res_lib.get_params(resname, at.id) #Gets the parameters of the atoms
        if not params:  # If it's not in an AA get the atom without charge
            print(at)
            at.xtra['atom_type'] = at.element
            at.xtra['charge'] = 0
        else:
        else:           #If it's in an AA assignt type and charge to the atom
            at.xtra['atom_type'] = params.at_type
            at.xtra['charge'] = params.charge
        at.xtra['vdw'] = ff_params.at_types[at.xtra['atom_type']]


add_atom_parameters(st, residue_library, ff_params)
    self.ntypes = len(self.at_types) # getting the length of at_types and storing it in a local variable
    fh.close() # closing file
add_atom_parameters()

