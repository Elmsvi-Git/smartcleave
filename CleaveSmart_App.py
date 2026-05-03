
import streamlit as st
import pandas as pd
import numpy as np
import joblib
from Bio.Seq import Seq


# --- 1. PAGE CONFIGURATION ---
st.set_page_config(page_title="CleaveSmart Web App", layout="wide")


# --- 2. MODEL LOADING ---
@st.cache_resource
def load_prediction_model():
    try:
        return joblib.load("classifier.pkl")
    except:
        return None

model = load_prediction_model()

# --- 3. CORE FUNCTIONS 
def generate_10_23_dnazymes(target_rna, left_arm_len, right_arm_len):
    catalytic_core = "GGCTAGCTACAACGA"
    candidates = []
    for i in range(right_arm_len, len(target_rna) - left_arm_len):
        if target_rna[i] in ['A', 'G'] and target_rna[i+1]== 'U':
            target_segment_L = target_rna[i - right_arm_len  : i ]
            target_segment_R = target_rna[i + 1 : i + left_arm_len+1]
            arm_L = str(Seq(target_segment_R).reverse_complement())
            arm_R = str(Seq(target_segment_L).reverse_complement())
            full_dz_sequence = arm_L + catalytic_core + arm_R
            candidates.append({
                "Position": i + 1, "Junction": f"{target_rna[i]}{target_rna[i+1]}",
                "Local_Context": f"...{target_rna[max(0,i-4):i+1]}[↓]{target_rna[i+1:i+5]}...",
                "DNAzyme_Sequence": full_dz_sequence,
   

            })
    return pd.DataFrame(candidates)

def calculate_features(dz_seq, target_rna, left_arm_len, right_arm_len):
    try:
        from Pipeline_Feature_Calculation import (Calculate_Binding_Energy, Calculate_Dz_Minimum_Free_Energy, 
                                       Calculate_Homodimer_Energy, Calculate_RNA_Cost_Energy)
    except ImportError:
        return {"Internal_Energy": 0, "Binding_Energy": 0, "Total_Energy": 0, "RNA_Cost_Energy": 0, "Homodimer_Energy": 0}

    catalytic_core = "GGCTAGCTACAACGA"
    arm_L_seq = dz_seq.split(catalytic_core)[0]
    arm_R_seq = dz_seq.split(catalytic_core)[1]

    rc_right = str(Seq(arm_R_seq).reverse_complement()).replace('T', 'U')
    rc_left = str(Seq(arm_L_seq).reverse_complement()).replace('T', 'U')
    binding_substrate = rc_right + rc_left 

    start_point = target_rna.find(rc_right)
    end_point = start_point + len(binding_substrate) if start_point != -1 else 0

    
    internal_energy = Calculate_Dz_Minimum_Free_Energy(str(dz_seq))
    homodimer_energy = Calculate_Homodimer_Energy(str(dz_seq))
    binding_energy = Calculate_Binding_Energy(binding_substrate)
    rna_cost_energy = Calculate_RNA_Cost_Energy(str(target_rna), start_point, end_point)
   

    return {
        "Internal_Energy": internal_energy,
        "Homodimer_Energy": homodimer_energy,
        "Binding_Energy": binding_energy,
        "RNA_Cost_Energy": rna_cost_energy,
        "Total_Energy": binding_energy - internal_energy + rna_cost_energy,
        
    }
st.markdown("""
    <style>
    html, body, [class*="st-"], button {
        font-family: "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif;
    }

    button[data-baseweb="tab"] {
        font-size: 18px !important;
        font-weight: 500 !important;
        color: #454545 !important;
        padding: 10px 20px !important;
        border-radius: 5px 5px 0 0 !important;
    }

    button[data-baseweb="tab"][aria-selected="true"] {
        color: #000000 !important;
        background-color: #f0f2f6 !important;
        border-bottom: 3px solid #262730 !important;
    }

    .justified-text {
        font-size: 16px;
        line-height: 1.7;
        color: #1f1f1f;
        text-align: justify;
    }

    h1, h2, h3 {
        color: #131b23;
        font-weight: 1000;
    }
    </style>
    """, unsafe_allow_html=True)
# --- 4. TABS SETUP ---
tab1, tab2, tab3, tab4 = st.tabs(["🏠 Home", "🧬 Analyze & Predict", "⚙️ Algorithm", "👥 About Us"])

# --- TAB 1: HOME ---
with tab1:

   
    col_logo, col_text = st.columns([2, 5]) 

    with col_logo:
        
        st.image('app_logo.png', use_container_width=True) 

    with col_text:
        st.title("Welcome to CleaveSmart Web App!")
        st.subheader("AI-Powered Selection of RNA-Cleaving DNAzyme 10-23 ")

    st.divider() 
    st.markdown("""
    ### ✔️ Overview
    🤖 This platform uses AI to predict RNA-cleavage efficiency of DNAzyme 10-23 . It all starts with your input. 
    Based on the RNA sequence (substrate) and the specific arm lengths you provide, the app scans for optimal target sites. 
    Specifically, it identifies GU and AU positions—the ***sweet spots*** where enzymes are known to be most effective. 

    The platform subsequently generates a library of candidate DNAzymes specifically optimized for these target sites. For each generated sequence, the system computes a comprehensive suite of thermodynamic and structural descriptors. 
    We focus on parameters that accurately predict the enzyme's cleavage efficiency within a 60-minute timeframe under standard conditions.
    The model classifies enzymes into two distinct groups based on their performance:
     Class 1 (High-efficiency): Predicted cleavage > 40%.
     *Class 0 (Low-efficiency): Predicted cleavage ≤ 40%.

    Beyond simple labels, the app calculates a **confidence score** (probability %). The higher the percentage, the more robust the prediction of its catalytic activity.
    """)

   

   
    st.markdown("### 📋 What you will receive:")

    col_feat1, col_feat2, col_feat3 = st.columns(3)

    with col_feat1:
        st.markdown("🎯 The Best Match")
        st.caption("Highlighting the cleavage site with the highest probability of success.")

    with col_feat2:
        st.markdown("📜 Comprehensive List")
        st.caption("A full list of all candidate enzymes and their predicted scores.")

    with col_feat3:
        st.markdown("💡 Interpretable AI")
        st.caption("Full transparency! We display the calculated features for every prediction.")

    st.success("Proceed to the **Analyze and Predict** tab to start selecting your enzymes!")
# --- TAB 2: ANALYZE & PREDICT 
with tab2:
    
    col_logo, col_text = st.columns([2, 5]) 

    with col_logo:
        
        st.image('app_logo.png', use_container_width=True) 

    with col_text:
        st.title("Welcome to CleaveSmart Pipeline:")
        st.subheader("From Modeling to DNAzyme Application")

    st.divider() 

    st.title("🧬 Analyze & Predict")
    st.markdown("#### Run AI-Driven DNAzyme Efficacy Prediction")

    default_seq = """AGATGGAGAGCCTTGTCCCTGGTTTCAACGAGAAAACACACGTCCAACTCAGTTTGCCTGTTTTACAGG
    TTCGCGACGTGCTCGTACGTGGCTTTGGAGACTCCGTGGAGGAGGTCTTATCAGAGGCACGTCAACATC
    TTAAAGATGGCACTTGTGGCTTAGTAGAAGTTGAAAAAGGCGTTTTGCCTCAACTTGAACAGCCCTATG
    TGTTCATCAAACGTTCGGATGCTCGAACTGCACCTCATGGTCATGTTATGGTTGAGCTGGTAGCAGAAC
    TCGAAGGCATTCAGTACGGTCGTAGTGGTGAGACACTTGGTGTCCTTGTCCCTCATGTGGGCGAAATAC
    CAGTGGCTTACCGCAAGGTTCTTCTTCGTAAGAACGGTAATAAAGGAGCTGGTGGCCATAGTTACGGCG
    CCGATCTAAAGTCATTTGACTTAGGCGACGAGCTTGGCACTGATCCTTATGAAGATTTTCAAGAAAACT
    GGAACACTAAACATAGCAGTGGTGTTACCCGTGAACTCATGCGTGAGCTTAACGGAGGG"""
    
    raw_input = st.text_area(
        "Enter Target RNA Sequence (5'->3'), max_chars=550",
        value=default_seq,
        height=250,
        max_chars=550
    )
    
    # Clean + normalize
    target_seq = (
        raw_input
        .replace("\n", "")
        .replace(" ", "")
        .upper()
        .replace("T", "U")
        .strip()
    )    
    
    import re

    target_seq = re.sub(r"\s+", "", raw_input).upper().replace("T", "U")
    # Validate
    invalid_chars = re.findall(r"[^AUCG]", target_seq)
    
    if invalid_chars:
        unique_invalid = sorted(set(invalid_chars))
        st.error(
            f"Invalid sequence detected \n\n"
            f"Invalid characters found: {', '.join(unique_invalid)}\n\n"
            f"Only A, T, C, G nucleotides are allowed."
        )
        st.stop()   # This prevents any code below from running
    
    # If we reach here → sequence is valid
    st.success("Sequence is valid ")
    st.write(f"Sequence length: {len(target_seq)} nt")

    st.write("Running downstream analysis...")
    l_arm = st.slider("Enter Left Arm Length (nt)", 5, 25, 16)
    r_arm = st.slider("Enter Right Arm Length (nt)", 5, 25, 8)

    
    

    if st.button("🚀 Run AI Analysis"):
        if model is None:
            st.error("Error: Prediction model (classifier.pkl) not found.")
        else:
            with st.spinner('Processing...\n This may take a while. Please be patient.'):
                df_candidates = generate_10_23_dnazymes(target_seq, l_arm, r_arm)
                if not df_candidates.empty:
                  
                    features_list = [calculate_features(s, target_seq, l_arm, r_arm) for s in df_candidates['DNAzyme_Sequence']]
                    df_features = pd.DataFrame(features_list)
                    X_input = df_features[["Internal_Energy", "Binding_Energy", "Homodimer_Energy", "RNA_Cost_Energy", "Total_Energy"]]

                    if X_input.empty:
                        st.error("No valid candidates remained after processing.")
                        st.stop()

                    df_candidates['Activity_Probability'] = model.predict_proba(X_input)[:, 1]
                    df_results = pd.concat([df_candidates, df_features], axis=1).sort_values(by='Activity_Probability', ascending=False)

                    
                    top = df_results.iloc[0]
                    st.success(f"Successfully evaluated {len(df_results)} cleavage sites!")
                    c1, c2, c3 = st.columns(3)
                    c1.metric("Target Position", int(top['Position']))
                    c2.metric("Best Junction", top['Junction'])
                    c3.metric("Top Activity Score", f"{top['Activity_Probability']:.2%}")
                    

                    st.dataframe(df_results.style.background_gradient(subset=['Activity_Probability'], cmap='RdYlGn'), use_container_width=True)
                else:
                    st.warning("No valid 10-23 cleavage sites (AU/GU) were found.")
    st.markdown("<br><hr>", unsafe_allow_html=True)
    st.caption("CleaveSmart AI | Developed for Molecular Biology Research")

        
# --- TAB 3: ALGORITHM ---
with tab3:
    st.header("⚙️ The CleaveSmart Algorithm")
    st.subheader("Predicting Activity through Four Intelligence Layers")

    st.markdown("""
    Our platform transforms raw RNA sequences into interpretable mechanistical insights. 
    The workflow is divided into four sophisticated stages:
    """)
    
    col_a, col_b = st.columns(2)

    with col_a:
        st.markdown("""
        #### 1️⃣ Precision Structural Configuration
        The process begins with an intelligent scan of your RNA substrate. By analyzing the 
        sequence against your specified arm lengths, CleaveSmart pinpoint-targets the 
        GU and AU junctions. These regions are the ***sweet spots*** where 
        DNAzyme 10-23 exhibits peak catalytic activity. The system then architecturally 
        engineers bespoke enzyme candidates for every viable site.

        #### 2️⃣ Molecular Feature Extraction
        For each candidate, the algorithm extracts a comprehensive set of "Digital Fingerprints". 
        These include unique thermodynamic and structural descriptors. 
        We focus on high-impact parameters that dictate cleavage kinetics within a 
        standardized 60-minute window.
        """)

    with col_b:
        st.markdown("""
        #### 3️⃣ Machine Learning Classification
        The extracted features are processed by our core ML engine. The model performs 
        a binary classification based on a 40% efficiency threshold:
         **Class 1 (High-Efficiency)**: Predicted cleavage > 40%.
         **Class 0 (Low-Efficiency)**: Predicted cleavage ≤ 40%.

        Beyond simple labels, the system assigns a confidence score (probability %), 
        quantifying the likelihood of an enzyme’s success in a laboratory setting.

        #### 4️⃣ Interpretable AI
        We move beyond opaque computations to embrace traceable, insight-driven AI. 
        CleaveSmart provides a clear rationale for every prediction:

         **Primary Candidate**: Immediate identification of the top-performing enzyme.
         **Global Candidate Map**: A full-scale library of all potential enzymes.
         **Logical Traceability**: We expose the underlying logic—such as free energy levels—allowing you to validate the model's decision-making process.
      """)



    st.divider()
    st.info("💡 Our algorithm is designed to bridge the gap between computational prediction and experimental validation.*")

# --- TAB 4: ABOUT US ---
with tab4:
     
    col1, col2, col3= st.columns([0.2, 1, 0.2])

   
 
    
    with col2:
        
        st.image("Logos.png",use_container_width=True)

    
      
        

    
    st.divider()
    st.header("Our Collaboration & Team")

    st.markdown("""
This project is the result of a multidisciplinary scientific collaboration aimed at bridging the gap between 
computational intelligence and molecular biology. This research was supported by the Iran National Science Foundation (INSF) under project No. 4037511.
""")

    st.write("The development and validation of this platform have been conducted through the joint efforts of:")
    st.markdown("""
- Isfahan University of Medical Sciences (School of Pharmacy and Pharmaceutical Sciences)
- University of Isfahan (UINA Lab: University of Isfahan Nucleic Acid Lab)
""")

    st.write("""
By combining expertise in biology, bioinformatics, and machine learning, our team strives to 
provide researchers with reliable tools for precision enzyme selection.
""")



    st.divider() 

# --- Section: Contact Us ---
    st.header("Contact Us")
    st.write("We value your feedback and are open to potential scientific collaborations.")

    

    st.markdown("##### Principal Investigators & Academic Correspondence")
    st.markdown("""
- Dr. Vajihe Akbari: v_akbari@pharm.mui.ac.ir, akbarivajihe@yahoo.com
- Dr. Elahe Mousavi: El.msvi@gmail.com
- Dr. Fatemeh Javadi-Zarnaghi: Fa.javadi@sci.ui.ac.ir, ui.biochemistry@gmail.com
""")
    st.markdown("UINA Lab Official Website: https://uina.ui.ac.ir/page-UINA1/en/216/form/Cl4572/link")
    st.caption("Designed and developed by Fatemeh Rahimpour")

   
    
