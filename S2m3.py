from Bio.PDB import *
from modules_classes import *
from S2m3_module2 import *

pdb_path = "/home/jj/Desktop/Bioinformatics/Github/Bioinformatics_p/Biophysics/Biophysics_A1/assignment_data/6m0j_fixed.pdb"
parser = PDBParser(PERMISSIVE=1)
structure = parser.get_structure('st', pdb_path)
dt=5
st=structure

residue_library = ResiduesDataLib('/home/jj/Desktop/Bioinformatics/Github/Bioinformatics_p/Biophysics/Biophysics_A1/assignment_data/parameters_step2.lib')
ff_params = VdwParamset('/home/jj/Desktop/Bioinformatics/Github/Bioinformatics_p/Biophysics/Biophysics_A1/assignment_data/parameters_vanderw.txt')
NACCESS_BINARY = '/home/jj/Desktop/Bioinformatics/Github/Bioinformatics_p/Biophysics/Biophysics_A1/soft/NACCESS/naccess'
srfA = NACCESS_atomic(st[0], naccess_binary=NACCESS_BINARY)

add_atom_parameters(st, residue_library,ff_params)

def chain_atoms(structure, dt):
    model =structure[0]
    t=0
    chain1_atoms=set()
    chain2_atoms=set()
    for chain in model:
        if t==0:
            t=1
            chain1=chain
        else: chain2=chain
    NeighborSearch_chain2 = NeighborSearch(list(chain2.get_atoms()))
    for residual1 in chain1.get_residues():
        for atom1 in residual1:
            Neighbor_atom_chain2 = NeighborSearch_chain2.search(atom1.coord, dt)
            for atom2 in Neighbor_atom_chain2:
                residual2=atom2.get_parent()
                chain1_atoms.add(residual1)
                chain2_atoms.add(residual2) 
    return calc(chain1_atoms, chain2_atoms)

def calc(chain1_atoms, chain2_atoms):
    s=sA=sE=e=v=0
    for res in chain1_atoms:
        s += calc_solvation(st[0], res)
        sA += calc_solvation(st[0]['A'], res)
        E, V = calc_int_energies(st, res)
        e+=E
        v+=V
    for res in chain2_atoms:
        s += calc_solvation(st[0], res)
        sE += calc_solvation(st[0]['E'], res)
        E, V = calc_int_energies(st, res)
        e+=E
        v+=V
    G = e + v + s - sA - sE
    return(G)

print(chain_atoms(structure, dt))
