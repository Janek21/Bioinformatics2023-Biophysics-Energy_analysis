from Bio.PDB import *
parser=PDBParser()

#DISCLAIMER ONLY WORKS FOR PROTEINS WITH 2 CHAINS
# IF YOU WANT IT TO WORK FOR MORE CHANGE if t==0 and make it a triple nested for loop (or search your solution)      
def chain_atoms(structure, dt):
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
    return chain_comparison(chain1_atoms, chain2_atoms, dt)
    

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

if __name__=='__main__':
    total_interface_residues=chain_atoms(structure)
    print(f"The interface residues are: {total_interface_residues}")
