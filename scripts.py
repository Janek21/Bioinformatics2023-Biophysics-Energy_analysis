
#EXECUTE IN TERMINAL, NOT VISUAL
#OTEHRWISE FOLLOW STEPS IN LINE 18

# Getting basic libraries:
import os, sys, math
from symbol import try_stmt
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
st_c.backbone()
#biobb_structure_checking


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
            r = Residue(data)

            self.residue_data[r.id] = r

        self.nres = len(self.residue_data)

    def get_params(self, resid, atid):
        atom_id = resid + ':' + atid

        if atom_id in self.residue_data:
            return self.residue_data[atom_id]
        
        else:
            print("WARNING: atom not found in library (", atom_id, ')')
            return None

class Residue():
    
    def __init__(self,data):
        self.id     = data[0]+':'+data[1]
        self.at_type = data[2]
        self.charge  = float(data[3])
        
class AtomType():

    def __init__(self, data):
        self.id   = data[0]
        self.eps  = float(data[1])
        self.sig  = float(data[2])
        self.mass = float(data[3])
        self.fsrf = float(data[4])
        self.rvdw = self.sig * 0.5612
        
class VdwParamset(): #extracted from GELPI's github
    #parameters for the VdW
    def __init__ (self, file_name):
        self.at_types = {}

        try:
            fh = open(file_name, "r")

        except OSError:
            print ("#ERROR while loading parameter file (", file_name, ")")
            sys.exit(2)

        for line in fh:

            if line[0] == '#':
                continue

            data = line.split()
            self.at_types[data[0]] = AtomType(data)

        self.ntypes = len(self.at_types)
        fh.close()