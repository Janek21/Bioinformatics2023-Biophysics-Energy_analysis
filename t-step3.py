# Calculate the  effect of replacing each residue by an alanine
with open("ala_scanning.txt", "a") as alanine:
    for ch in st[0]:
        for res in ch.get_residues():
            if max_dist > 0 and res not in interface[ch.id]:
                continue
            print('{:1}, {:1.4}'.format(residue_id(res), 
                    - elec[res] + elec_ala[res] - vdw[res] + vdw_ala[res] - solvAB[res] +\
                        solvAB_ala[res] -solvA[res] + solvA_ala[res]), file = alanine)
