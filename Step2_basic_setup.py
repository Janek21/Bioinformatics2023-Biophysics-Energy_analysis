#! /usr/bin/python3

"""
    Initial setup a structure for Energy evaluation
    To execute the script --> python Step2_basic_setup.py --rlib assignment_data/parameters_step2.lib --vdw assignment_data/parameters_vanderw.txt assignment_data/6m0j.pdb
"""
import argparse
import sys

from Bio.PDB.PDBParser import PDBParser
from modules_classes import ResiduesDataLib
from modules_classes import VdwParamset
from Bio.PDB.NACCESS import NACCESS_atomic


parser = argparse.ArgumentParser(
    prog='structure_setup',
    description='basic structure setup'
)

parser.add_argument(
    '--rlib',
    action='store',
    dest='reslib_file',
    default='assignment_data/parameters_step2.lib',
    help='Residue Library'
)

parser.add_argument(
    '--vdw',
    action='store',
    dest='vdwprm_file',
    default='assignment_data/parameters_vanderw.txt',
    help='Vdw parameters'
)

parser.add_argument('pdb_file', help='Input PDB', type=open)

args = parser.parse_args()

print("PDB.filename:", args.pdb_file.name)
print("Residue Lib.:", args.reslib_file)
print("PDB.filename:", args.vdwprm_file)

# Loading Libraries
# loading residue library from data/aaLib.lib
residue_library = ResiduesDataLib(args.reslib_file)

# loading VdW parameters
ff_params = VdwParamset(args.vdwprm_file)

parser = PDBParser(PERMISSIVE=1)
print('Parsing', args.pdb_file)
# load structure from PDB file of PDB ifle handler
st = parser.get_structure('STR', args.pdb_file.name)

# assign data types, and charges from libraries
# We will use the xtra attribute in Bio.PDB.Atom to hold the new data
# Possible errors on N-term and C-Term atoms
# Possible errors on HIS alternative forms

for at in st.get_atoms():
    resname = at.get_parent().get_resname()
    params = residue_library.get_params(resname, at.id)
    if not params:
        sys.exit("ERROR: residue/atom pair not in library (" + resname + ' ' + at.id + ')')
    at.xtra['atom_type'] = params.at_type
    at.xtra['charge'] = params.charge
    at.xtra['vdw'] = ff_params.at_types[at.xtra['atom_type']]


# Calculating surfaces
# The specific PATH to naccess script (in soft) is needed
# Srf goes to .xtra field directly
'''
NACCESS_BINARY = '/home/jaume/Desktop/Year_2/Biophysics/Project:Energy_Analysis/Biophysics_A1/soft/NACCESS'

!sudo srf = NACCESS_atomic(st[0], naccess_binary = NACCESS_BINARY)
'''

import subprocess

NACCESS_BINARY = '/home/jaume/Desktop/Year_2/Biophysics/Project:Energy_Analysis/Biophysics_A1/soft/NACCESS'

# Assuming st[0] is the argument you want to pass to your command
command = f'sudo {NACCESS_BINARY} {st[0]}'

# Run the command
subprocess.run(command, shell=True)