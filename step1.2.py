from Bio.PDB import *
parser=PDBParser()
structure= parser.get_structure("6m0j", "./assignment_data/6m0j_fixed.pdb")
model =structure[0]
#chain=model["A"]
ls=[]

#same resseq and/or name in both chains
#or same position? NO
#make dict name:resseq
for c in model:
    print("For chain ", c)
    for i in c.get_residues():
        print(i)
        i=str(i)
        i=i.split()
        i=i[1]
        ls.append(i)
