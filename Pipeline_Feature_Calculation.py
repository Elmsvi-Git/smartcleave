import RNA
# compute RNA/DNA duplex ΔG at 37 degrees (Binding Energy)
def Calculate_Binding_Energy(seq):
    di=[]
    e=[]
    energy = {'AA':-0.7, 'AC':-1.5, 'AG':-1.3, 'AU':-0.4, 'CA':-1.2, 'CC':-1.7, 'CG':-1.4, 'CU':-0.4, 'GA':-1.5, 'GC':-2, 'GG':-2.3, 'GU':-1.4, 'UA':-0.5, 'UC':-1.4, 'UG':-1.6, 'UU':0.2}
    for i in range (len(seq)):
        di.append(''.join(seq[i:i+2]))
    di.pop(-1)
    for i in range (len(di)):
        e.append(energy[di[i]])
    a=sum(e)
    if di[0][0] in ['G','C'] or di[-1][1] in ['G','C']:
        a+=2
    else:
        a+=2.6
    return a

# compute RNA/DNA duplex ΔG at 23 degrees (Binding Energy)
def delta_G_banerjee(seq):
    di=[]
    e=[]
    energy = {'AA':-1.02, 'AC':-2.02, 'AG':-1.64, 'AU':-0.62, 'CA':-1.67, 'CC':-2.16, 'CG':-1.8, 'CU':-0.8, 'GA':-1.9, 'GC':-2.4, 'GG':-2.73, 'GU':-1.75, 'UA':-0.77, 'UC':-1.67, 'UG':-2, 'UU':0.2}
    for i in range (len(seq)):
        di.append(''.join(seq[i:i+2]))
    di.pop(-1)
    for i in range (len(di)):
        e.append(energy[di[i]])
    a=sum(e)
    if di[0][0] in ['G','C'] or di[-1][1] in ['G','C']:
        a+=1.45
    else:
        a+=2.07
    return a

# compute minimum free energy (MFE) and corresponding structure of Dz (Dz internal energy)
def Calculate_Dz_Minimum_Free_Energy(seq):
    RNA.cvar.temperature = 37.0
    RNA.params_load_DNA_Mathews2004()
    fc = RNA.fold_compound(seq)
    
    structure,mfe = fc.mfe()
    return mfe

# compute minimum free energy (MFE) and corresponding structure of Dz homodimer
def Calculate_Homodimer_Energy(seq):
    RNA.params_load_DNA_Mathews2004()
    RNA.cvar.temperature = 37.0
    dimer_sequence = str(seq) + "&" + str(seq)
    (structure, energy) = RNA.cofold(dimer_sequence) 
    return energy

# Calculate the local free energy cost of keeping a specified region unpaired in an RNA sequence.
def Calculate_RNA_Cost_Energy(sequence, region_start, region_end):
    RNA.cvar.temperature = 37.0
    fc = RNA.fold_compound(sequence)
    (ss_pf, ensemble_energy) = fc.pf()
    constraint = ["."] * len(sequence)
    for i in range(region_start, region_end+1):
        constraint[i] = "x"  
    constraint_str = "".join(constraint)
    fc.hc_add_from_db(constraint_str)
    (ss_pf_constrained, ensemble_energy_constrained) = fc.pf()
    RNA_Cost_Energy = ensemble_energy_constrained - ensemble_energy
    return RNA_Cost_Energy