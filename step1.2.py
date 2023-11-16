from Bio.PDB import *
parser=PDBParser()
structure= parser.get_structure("6m0j", "./assignment_data/6m0j_fixed.pdb")
model =structure[0]
t=0
chain1_atoms=[]
chain2_atoms=[]
#same resseq and/or name in both chains
for chain in model:
    if t==0:
        t=1
        for residual in chain.get_residues():
            for atom in residual:
                chain1_atoms.append(residual)
    else:
        for residual in chain.get_residues():
            for atom in residual:
                chain2_atoms.append(residual)
            i=str(i)
            i=i.split()
            i=i[1]
            ls.append(i)
