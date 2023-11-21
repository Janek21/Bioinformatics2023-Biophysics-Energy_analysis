
import argparse
import os

from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.NACCESS import NACCESS_atomic
from Bio.PDB.PDBIO import PDBIO, Select

from modules_classes import ResiduesDataLib
from modules_classes import VdwParamset
import step2_energies as en

#n_input=str(input("Input path to soft/NACCESS(including this folder)"))
NACCESS_BINARY = '/home/jj/Desktop/Bioinformatics/Github/Bioinformatics_p/Biophysics/Biophysics_A1/soft/NACCESS'

parse_cmd = argparse.ArgumentParser(
	prog='binding',
	description='binding energy calculation'
)

parse_cmd.add_argument(
	'--rlib',
	action='store',
	dest='reslib_file',
	default='/home/jj/Desktop/Bioinformatics/Github/Bioinformatics_p/Biophysics/Biophysics_A1/assignment_data/parameters_step2.lib',
	help='Residue Library'
)
parse_cmd.add_argument(
	'--vdw',
	action='store',
	dest='vdwprm_file',
	default='/home/jj/Desktop/Bioinformatics/Github/Bioinformatics_p/Biophysics/Biophysics_A1/assignment_data/parameters_vanderw.txt',
	help='Vdw parameters'
)

parse_cmd.add_argument(
	'--dist',
	action='store',
	dest='cutoff_dist',
	default=8.0,
	type=float,
	help='Cutoff distance for determining the interface (0: use all residues):'
)
parse_cmd.add_argument('pdb_file', help='Input PDB', type=open)

args = parse_cmd.parse_args()
print(args)
print("PDB.filename:", args.pdb_file.name)
print("Residue Lib.:", args.reslib_file)
print("PDB.filename:", args.vdwprm_file)
print("Distance:", args.cutoff_dist)

# Loading Libraries
# loading residue library from data/aaLib.lib
residue_library = ResiduesDataLib(args.reslib_file)

# loading VdW parameters
ff_params = VdwParamset(args.vdwprm_file)

parser = PDBParser(PERMISSIVE=1)
print('Parsing', args.pdb_file)
# load structure from PDB file of PDB ifle handler
st = parser.get_structure('STR', args.pdb_file.name)
print(args.pdb_file.name)
print(st)

'''
PDB.filename: ./assignment_data/6m0j_fixed.pdb
Residue Lib.: ./assignment_data/parameters_step2.lib
PDB.filename: ./assignment_data/parameters_vanderw.txt
Distance: 8.0
'''
'''
res_lib=residue_library


def add_atom_parameters(st, res_lib, ff_params):
     	#Adds parameters from libraries to atom .xtra field
        #For not recognized atoms, issues a warning and put default parameters
    for at in st.get_atoms():
        resname = at.get_parent().get_resname()
        params = res_lib.get_params(resname, at.id)
        if not params:
            #print("WARNING: residue/atom pair not in library ("+atom_id(at) + ')')
            at.xtra['atom_type'] = at.element
            at.xtra['charge'] = 0
        else:
            at.xtra['atom_type'] = params.at_type
            at.xtra['charge'] = params.charge
        at.xtra['vdw'] = ff_params.at_types[at.xtra['atom_type']]
        
add_atom_parameters(st, residue_library, ff_params)

srf = NACCESS_atomic(st[0], naccess_binary=NACCESS_BINARY)
print(srf)
'''