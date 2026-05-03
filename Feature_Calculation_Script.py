# compute RNA/DNA duplex ΔG (Binding Energy)
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

# compute minimum free energy (MFE) and corresponding structure of Dz (Dz internal energy)
def Calculate_Dz_Minimum_Free_Energy(seq):
    import RNA
    RNA.cvar.temperature = 37.0
    RNA.params_load_DNA_Mathews2004()
    fc = RNA.fold_compound(seq)
    
    structure,mfe = fc.mfe()
    return mfe

# compute minimum free energy (MFE) and corresponding structure of Dz homodimer
def Calculate_Homodimer_Energy(seq):
    import RNA
    RNA.params_load_DNA_Mathews2004()
    RNA.cvar.temperature = 37.0
    dimer_sequence = seq + "&" + seq
    (structure, energy) = RNA.cofold(dimer_sequence) 
    return energy

# Calculate the local free energy cost of keeping a specified region unpaired in an RNA sequence.
def Calculate_RNA_Cost_Energy(sequence, region_start, region_end):
    import RNA
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

# Calculate the Total_Energy
def Calculate_Total_Energy(Binding_Energy, Internal_Energy,RNA_Cost_Energy):
    Total_Energy= Binding_Energy-Internal_Energy+RNA_Cost_Energy
    return Total_Energy

# Defining Binding_region
def Defining_Binding_Region(left_arm, right_arm):
    from Bio.Seq import Seq
    binding_region=str(Seq(left_arm).reverse_complement_rna()+Seq(right_arm).reverse_complement_rna())
    return binding_region

# defining Binding_region_Start_point
def Finding_Binding_Region_Start_point(rna_sequence,binding_region):
    start_point=rna_sequence.find(binding_region)
    return start_point

# defining Binding_region_End_point
def Finding_Binding_Region_End_point(start_point,binding_region):
    end_point=start_point+len(binding_region)+1
    return end_point	

def Get_full_dnazyme_results(params):
    import requests
    from bs4 import BeautifulSoup
    import pandas as pd
    from io import BytesIO

    base_url = "https://iimcb.genesilico.pl/DNAzymeBuilder/dnazymebuilder"

    r = requests.get(base_url, params=params)

    soup = BeautifulSoup(r.text, "lxml")

    csv_link = soup.find("a", href=lambda x: x and "/api/download/csv/" in x)

    if not csv_link:
        raise Exception("CSV download link not found")

    csv_url = "https://iimcb.genesilico.pl" + csv_link["href"]

    r_csv = requests.get(csv_url)

    return pd.read_csv(BytesIO(r_csv.content), encoding="utf-8")

def analyze_dnazyme_pipeline(params):
    """
    Full analysis pipeline:
    - Submits RNA to DNAzymeBuilder
    - Downloads full results
    - Computes thermodynamic features
    - Returns enriched DataFrame
    """

    rna_sequence = params["sequence"]

    df = Get_full_dnazyme_results(params)
    df['rna_sequence']=rna_sequence
    df['binding_region'] = df.apply(lambda row:Defining_Binding_Region (row['left_arm'], row['right_arm']), axis=1)
    df['region_start'] = df.apply(lambda row:Finding_Binding_Region_Start_point (row['rna_sequence'], row['binding_region']), axis=1)
    df['region_end'] = df.apply(lambda row:Finding_Binding_Region_End_point (row['region_start'], row['binding_region']), axis=1)
    

    Binding_Energy_list= []
    Internal_Energy_list= []
    Homodimer_Energy_list = []
    RNA_Cost_Energy_list = []
    Total_Energy_list = []

    
    for _, row in df.iterrows():

        dz_seq=row["seq"]  # adjust column name if needed
        right_arm=row["right_arm"]
        left_arm=row["left_arm"]
        cleaving_pos= int(row["cleaving_pos"])
        binding_region=row['binding_region']
        region_start=row['region_start']
        region_end=row['region_end']

        # Binding Energy
        Binding_Energy = Calculate_Binding_Energy(binding_region)

        # DNAzyme minimum free energy
        Internal_Energy = Calculate_Dz_Minimum_Free_Energy(dz_seq)

        # Homodimer energy
        Homodimer_Energy = Calculate_Homodimer_Energy(dz_seq)

        # RNA accessibility cost (region around cleavage)
        region_start = max(0, cleaving_pos - 3)
        region_end = min(len(rna_sequence), cleaving_pos + 3)


        RNA_Cost_Energy = Calculate_RNA_Cost_Energy(
                rna_sequence,
                region_start,
                region_end
        )

        # total_energy
        Total_Energy = Calculate_Total_Energy(Binding_Energy,Internal_Energy,RNA_Cost_Energy) 


        Binding_Energy_list.append(Binding_Energy)
        Internal_Energy_list.append(Internal_Energy)
        Homodimer_Energy_list.append(Homodimer_Energy)
        RNA_Cost_Energy_list.append(RNA_Cost_Energy)
        
        
        Total_Energy_list.append(Total_Energy)

    # 4️⃣ Add computed columns
    df["Binding_Energy"] = Binding_Energy_list
    df["Internal_Energy"] = Internal_Energy_list
    df["Homodimer_Energy"] = Homodimer_Energy_list
    df["RNA_Cost_Energy"] = RNA_Cost_Energy_list
    
    df['Total_Energy'] = Total_Energy_list
    
    return df



