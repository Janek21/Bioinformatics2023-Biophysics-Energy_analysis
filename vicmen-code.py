import argparse
import sys
import os
import math

from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.NACCESS import NACCESS_atomic
from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.PDBIO import PDBIO, Select



def residue_id(res):
    '''Returns readable residue id'''
    return '{} {}{}'.format(res.get_resname(), res.get_parent().id, res.id[1])

def atom_id(at):
    '''Returns readable atom id'''
    return '{}.{}'.format(residue_id(at.get_parent()), at.id)


class ResiduesDataLib():
    def __init__(self, fname):
        self.residue_data = {}
        try:
            fh = open(fname, "r")
        except OSError:
            print("#ERROR while loading library file (", fname, ")")
            sys.exit(2)
        for line in fh:
            if line[0] == '#':
                continue
            data = line.split()
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


def get_interface_residues(st, id1, id2, distance):
    interface_chain1 = set()
    interface_chain2 = set()
    chain1 = st[0][id1]    ## Get the desired chains from the protein
    chain2 = st[0][id2]
    NeighborSearch_chain2 = NeighborSearch(list(chain2.get_atoms()))    ## Prepare the Neighbour search for chain E
    for res_chain1 in chain1:
        for atom_chain1 in res_chain1:    ## Iterate over all the atoms from chain A
            Neighbor_atom_chain2 = NeighborSearch_chain2.search(atom_chain1.coord, distance)    ## Look for atoms in chain E within the distance from the atom in chain A
            for atom_chain2 in Neighbor_atom_chain2:    ## Itarate over the atoms we know they are within the distance
                res_chain2 = atom_chain2.get_parent()   ## Get the residues to which the atom belongs to
                interface_chain1.add(res_chain1)        ## Use a set to only save each residue once
                interface_chain2.add(res_chain2)
    return list(interface_chain1), list(interface_chain2)    ## Return the list for each chain



def add_atom_parameters(st, res_lib, ff_params):
    ''' Adds parameters from libraries to atom .xtra field
        For not recognized atoms, issues a warning and put default parameters
    '''
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
def MH_diel(r):
    '''Mehler-Solmajer dielectric'''
    return 86.9525 / (1 - 7.7839 * math.exp(-0.3153 * r)) - 8.5525

def elec_int(at1, at2, r):
    '''Electrostatic interaction energy between two atoms at r distance'''
    print(at1.xtra['charge'])
    return 332.16 * at1.xtra['charge'] * at2.xtra['charge'] / MH_diel(r) / r

def vdw_int(at1, at2, r):
    '''Vdw interaction energy between two atoms'''
    eps12 = math.sqrt(at1.xtra['vdw'].eps * at2.xtra['vdw'].eps)
    sig12_2 = at1.xtra['vdw'].sig * at2.xtra['vdw'].sig
    return 4 * eps12 * (sig12_2**6/r**12 - sig12_2**3/r**6)
    
def calc_solvation_ala(st, res):
    '''Solvation energy based on ASA'''
    solv = 0.
    solv_ala = 0.
    for at in res.get_atoms():
        if 'EXP_NACCESS' not in at.xtra:
            continue
        s = float(at.xtra['EXP_NACCESS'])* at.xtra['vdw'].fsrf
        solv += s
        if at.id in ala_atoms:
            solv_ala += s
    return solv, solv_ala

def calc_int_energies_ala(st, res):
    '''Returns interaction energies (residue against other chains)
        for all atoms and for Ala atoms
    '''
    elec = 0.
    elec_ala = 0.
    vdw = 0.
    vdw_ala = 0.

    for at1 in res.get_atoms():
        for at2 in st.get_atoms():
        # skip same chain atom pairs
            if at2.get_parent().get_parent() != res.get_parent():
                r = at1 - at2
                e = elec_int(at1, at2, r)
                elec += e
                if at1.id in ala_atoms: #GLY are included implicitly
                    elec_ala += e
                e = vdw_int(at1, at2, r)
                vdw += e
                if at1.id in ala_atoms: #GLY are included implicitly
                    vdw_ala += e
    return elec, elec_ala, vdw, vdw_ala

def calc_solvation(st, res):
    '''Solvation energy based on ASA'''
    solv = 0.
    for at in res.get_atoms():
        if 'EXP_NACCESS' not in at.xtra:
            continue
        s = float(at.xtra['EXP_NACCESS'])* at.xtra['vdw'].fsrf
        solv += s
    return solv

def calc_int_energies(st, res):
    elec = 0.
    vdw = 0.

    for at1 in res.get_atoms():
        for at2 in st.get_atoms():
        # skip same chain atom pairs
            if at2.get_parent().get_parent() != res.get_parent():
                r = at1 - at2
                e = elec_int(at1, at2, r)
                elec += e
                e = vdw_int(at1, at2, r)
                vdw += e
    return elec, vdw


pdb_path = "/home/jj/Desktop/Bioinformatics/Github/Bioinformatics_p/Biophysics/Biophysics_A1/assignment_data/6m0j_fixed.pdb"
parser = PDBParser(PERMISSIVE=1)
st = parser.get_structure('st', pdb_path)

residue_library = ResiduesDataLib('/home/jj/Desktop/Bioinformatics/Github/Bioinformatics_p/Biophysics/Biophysics_A1/assignment_data/parameters_step2.lib')
ff_params = VdwParamset('/home/jj/Desktop/Bioinformatics/Github/Bioinformatics_p/Biophysics/Biophysics_A1/assignment_data/parameters_vanderw.txt')
NACCESS_BINARY = '/home/jj/Desktop/Bioinformatics/Github/Bioinformatics_p/Biophysics/Biophysics_A1/soft/NACCESS/naccess'
srfA = NACCESS_atomic(st[0], naccess_binary=NACCESS_BINARY)

add_atom_parameters(st, residue_library,ff_params)

interface_residues_chainA, interface_residues_chainE = get_interface_residues(st, 'A', 'E', 4)    ## Call the function with the desired arguments

'''
print(f"Interface residues in chain A: {interface_residues_chainA}")    ## List of residues in chain A within the distance from chain E
print(f"Interface residues in chain E: {interface_residues_chainE}")    ## List of residues in chain E within the distance from chain A
'''

s=sA=sE=e=v=0
for res in interface_residues_chainA:
    s += calc_solvation(st[0], res)
    sA += calc_solvation(st[0]['A'], res)
    E, V = calc_int_energies(st, res)
    e+=E
    v+=V

for res in interface_residues_chainE:
    s += calc_solvation(st[0], res)
    sE += calc_solvation(st[0]['E'], res)
    E, V = calc_int_energies(st, res)
    e+=E
    v+=V

G = e + v + s - sA - sE
print(G)

ala_atoms = {'N', 'H', 'CA', 'HA', 'C', 'O', 'CB', 'HB', 'HB1', 'HB2', 'HB3', 'HA1', 'HA2', 'HA3'}

