
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
#finding out atoms in both chains
def chain_atoms(structure, dt):
	model = structure[0]
	t = 0
	chain1_atoms = []
	chain2_atoms = []
	for chain in model:
		if t==0:
			t=1
			for residual in chain.get_residues():
				for atom1 in residual:
					chain1_atoms.append(atom1)
		else:
			for residual in chain.get_residues():
				for atom2 in residual:
					chain2_atoms.append(atom2)
	return chain_comparison(chain1_atoms, chain2_atoms, dt)

#comparing atoms of both chains
def chain_comparison(chain1_atoms, chain2_atoms, dt):
    interface_residues=set()
    for a1 in chain1_atoms:
        for a2 in chain2_atoms:
            if a1.get_parent().id==a1.get_parent().id and a2.get_parent().id == a2.get_parent().id:
                distance=a1-a2
            if distance<=dt:#dt is the treshold of the distance
                interface_residues.add(a1.get_parent().id[1])
                interface_residues.add(a2.get_parent().id[1])
    return interface_residues

dt=5 #distance treshold
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

# ----------------------------------------------------------------

# <----> \ | / <---> | STEP 3 | <---> \ | / <---->

# 
"""
	Initial setup a structure for Energy evaluation
"""
import argparse
import os

import argparse
import os

from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.NACCESS import NACCESS_atomic
from Bio.PDB.PDBIO import PDBIO, Select

from modules_classes import ResiduesDataLib
from modules_classes import VdwParamset
import step2m2_energies as en

# Loading Libraries
# loading residue library from data/aaLib.lib
dir = os.getcwd()

residue_library = ResiduesDataLib(dir+'/assignment_data/parameters_step2.lib')

# loading VdW parameters
ff_params = VdwParamset(dir+'/assignment_data/parameters_vanderw.txt')

# Important variables:
pdb_file = dir+'/assignment_data/6m0j_fixed.pdb'
cutoff_dist = 5.0
NACCESS_BINARY = dir+'/soft/NACCESS/naccess'


parser = PDBParser(PERMISSIVE=1)
print('Parsing', pdb_file)

# load structure from PDB file of PDB ifle handler
st = parser.get_structure('STR', pdb_file)

# assign data types, and charges from libraries
# We will use the xtra attribute in Bio.PDB.Atom to hold the new data
# Possible errors on N-term and C-Term atoms
# Possible errors on HIS alternative forms

en.add_atom_parameters(st, residue_library, ff_params)

# Calculating surfaces
# The specific PATH to naccess script (in soft) is needed
# ASA goes to .xtra field directly

srf = NACCESS_atomic(st[0], naccess_binary=NACCESS_BINARY)

# Prepare surfaces for the separate chains
# Alternatively the twp PDB files can be prepared outside and parsed here

io = PDBIO()
st_chains = {}
# Using BioIO trick (see tutorial) to select chains
class SelectChain(Select):
	def __init__(self, chid):
		self.id = chid

	def accept_chain(self, chain):
		if chain.id == self.id:
			return 1
		else:
			return 0

for ch in st[0]:
	io.set_structure(st)
	io.save('tmp.pdb', SelectChain(ch.id))
	st_chains[ch.id] = parser.get_structure('stA', 'tmp.pdb')
	en.add_atom_parameters(st_chains[ch.id], residue_library, ff_params)
	srfA = NACCESS_atomic(st_chains[ch.id][0], naccess_binary=NACCESS_BINARY)
os.remove('tmp.pdb')

## Interface residues
if cutoff_dist > 0.:
	interface = en.get_interface(st, cutoff_dist)

## Initiatlize Energy aggregates
elec = {}
elec_ala = {}

vdw = {}
vdw_ala = {}

solvAB = {}
solvAB_ala = {}

solvA = {}
solvA_ala = {}

totalIntElec = 0.
totalIntVdw = 0.
totalSolv = 0.
totalSolvMon = {}

## We get the chsin ids,not always they are A and B
chids = []
for ch in st[0]:
	chids.append(ch.id)
	totalSolvMon[ch.id] = 0

total = 0.


print(f'\nInteraction energy based in interface residues ONLY')

with open("inter_en_res.csv", "w") as file:

	print(
		'D#{:11}  {:11s} {:11s} {:11s} {:11s} | {:11s} {:11s} {:11s} {:11s}'.format(
			'res_id',
			'elec_res', 'vdw_res', 'solv_AB_res', 'solv_A_res',
			'elec_ala', 'vdw_ala', 'solv_AB_ala', 'solv_A_ala'))

	file.write(
		'D#{:11},{:11s},{:11s},{:11s},{:11s}, - ,{:11s},{:11s},{:11s},{:11s}\n'.format(
			'res_id',
			'elec_res', 'vdw_res', 'solv_AB_res', 'solv_A_res',
			'elec_ala', 'vdw_ala', 'solv_AB_ala', 'solv_A_ala'))

	for ch in st[0]:
		for res in ch.get_residues():
			if cutoff_dist > 0 and res not in interface[ch.id]:
				continue
			elec[res], elec_ala[res], vdw[res], vdw_ala[res] = en.calc_int_energies(st[0], res)
			solvAB[res], solvAB_ala[res] = en.calc_solvation(st[0], res)
			solvA[res], solvA_ala[res] = en.calc_solvation(
				st_chains[ch.id],
				st_chains[ch.id][0][ch.id][res.id[1]]
			)
			totalIntElec += elec[res]
			totalIntVdw += vdw[res]
			totalSolv += solvAB[res]
			totalSolvMon[ch.id] += solvA[res]
			total += elec[res] + vdw[res] + solvAB[res] - solvA[res]

			print(
				'D#{:11}  {:11.4f} {:11.4f} {:11.4f} {:11.4f} | {:11.4f} {:11.4f} {:11.4f} {:11.4f}'.format(
					en.residue_id(res),
					elec[res], vdw[res], solvAB[res], solvA[res],
					elec_ala[res], vdw_ala[res], solvAB_ala[res], solvA_ala[res]))
			
			file.write(
				'D#{:11},{:11.4f},{:11.4f},{:11.4f},{:11.4f}, - ,{:11.4f},{:11.4f},{:11.4f},{:11.4f}\n'.format(
					en.residue_id(res),
					elec[res], vdw[res], solvAB[res], solvA[res],
					elec_ala[res], vdw_ala[res], solvAB_ala[res], solvA_ala[res]))

print(f'\nTOTAL ENERGIES of interaction energy based in interface residues ONLY')
print('{:20}: {:11.4f}'.format('Total Elec Int.', totalIntElec))
print('{:20}: {:11.4f}'.format('Total Vdw Int.', totalIntVdw))
print('{:20}: {:11.4f}'.format('Total Solv AB', totalSolv))
print('{:19}{}: {:11.4f}'.format('Total Solv ', chids[0], totalSolvMon[chids[0]]))
print('{:19}{}: {:11.4f}'.format('Total Solv ', chids[1], totalSolvMon[chids[1]]))
print('{:20}: {:11.4f}'.format('DGintAB-A-B', total))
print('\nWritting done in file: inter_en_res.tsv\n\n')


print(f'\nAla Scanning: DDGs for X->Ala mutations on interface residues')
with open("ala_scaning.csv", "w") as file:
	file.write(
		'{:11},{:11s},{:11s},{:11s},{:11s},{:11s}\n'.format(
		'res_id',
		'elec',
		'vdw',
		'solvAB',
		'solv',
		'total'))

	print(
		'{:11}  {:11s} {:11s} {:11s} {:11s} {:11s}'.format(
			'res_id',
			'elec',
			'vdw',
			'solvAB',
			'solv',
			'total'))
	
	for ch in st[0]:
		for res in ch.get_residues():
			if cutoff_dist > 0 and res not in interface[ch.id]:
				continue

			print(
				'{:11}  {:11.4f} {:11.4f} {:11.4f} {:11.4f} {:11.4f}'.format(
					en.residue_id(res),
					- elec[res] + elec_ala[res],
					- vdw[res] + vdw_ala[res],
					- solvAB[res] + solvAB_ala[res],
					- solvA[res] + solvA_ala[res],
					- elec[res] + elec_ala[res] - vdw[res] + vdw_ala[res] -solvAB[res] +\
					solvAB_ala[res] -solvA[res] + solvA_ala[res]))

			file.write(
				'{:11},{:11.4f},{:11.4f},{:11.4f},{:11.4f},{:11.4f}\n'.format(
					en.residue_id(res),
					- elec[res] + elec_ala[res],
		   		 - vdw[res] + vdw_ala[res],
		  		  - solvAB[res] + solvAB_ala[res],
		  		  - solvA[res] + solvA_ala[res],
		  		  - elec[res] + elec_ala[res] - vdw[res] + vdw_ala[res] -solvAB[res] +\
		  			  solvAB_ala[res] -solvA[res] + solvA_ala[res]))
print('Writing  done in file: ala_scaning.tsv')
