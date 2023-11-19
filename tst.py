from Bio.PDB import *
parser=PDBParser()
structure= parser.get_structure("6m0j", "./assignment_data/6m0j_fixed.pdb")
model =structure[0]
t=0
chain1_atoms=[]
chain2_atoms=[]
distance_treshold=7
            
def chain_atoms(structure):
    model =structure[0]
    t=0
    chain1_atoms=[]
    chain2_atoms=[]
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
    return chain_comparison(chain1_atoms, chain2_atoms)
    
def chain_comparison(chain1_atoms, chain2_atoms):
    interface_residues=set()
    for a1 in chain1_atoms:
        for a2 in chain2_atoms:
            if a1.get_parent().id==a1.get_parent().id and a2.get_parent().id == a2.get_parent().id:
                distance=a1-a2
                #NO DISTANCES, WHY EMPTY
            if distance<=distance_treshold:
                interface_residues.add(a1.get_parent().id[1])
                interface_residues.add(a2.get_parent().id[1])
    return interface_residues

int_re=chain_atoms(structure)
print(f"Interface residues : {int_re}")