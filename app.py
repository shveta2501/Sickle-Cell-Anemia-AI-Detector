import streamlit as st
from Bio import SeqIO
from io import StringIO

# 1. Page Configuration
st.set_page_config(page_title="Sickle Cell Detector", layout="centered")
st.title("🧬 AI Genetic Disease Detection Platform")
st.write("Variant Interpretation Dashboard – Developer: Shveta")

st.markdown("---")

# Reference Flanking Sequences for HBB Codon 6 (Normal vs Sickle)
# Normal (HbA):   CCT GAG GAG  (Pro-Glu-Glu) -> Flanked by CTG GGC AGG TGG TCT C...
# Mutant (HbS):   CCT GTG GAG  (Pro-Val-Glu)
HBA_TARGET = "CCTGAGGAG"
HBS_TARGET = "CCTGTGGAG"


# 2. File Uploader
uploaded_file = st.file_uploader("Upload Patient FASTA File", type=["fasta", "fa"])

if uploaded_file is not None:
    st.success("File uploaded successfully!")
    
    # 3. Stream to Text Conversion
    raw_text = uploaded_file.getvalue().decode("utf-8")
    data_stream = StringIO(raw_text)
    
    # 4. Parsing Engine
    for record in SeqIO.parse(data_stream, "fasta"):
        patient_id = record.id
        # Remove whitespaces/newlines for clean string matching
        dna_sequence = "".join(str(record.seq).upper().split())
        
        # Display Patient Summary
        st.subheader("📊 Patient Genomic Analysis Summary")
        col1, col2 = st.columns(2)
        with col1:
            st.metric(label="Patient/Record ID", value=patient_id)
        with col2:
            st.metric(label="Sequence Length", value=f"{len(dna_sequence)} bp")
            
        st.text_area("Parsed Raw DNA Sequence (Preview)", dna_sequence, height=120)
        
        st.markdown("---")
        
        # 5. Diagnostic Prediction (Targeted Matching)
        st.subheader("🔬 Diagnostic Mutation Analysis")
        
        # Check specifically for the HbS point mutation context
        if HBS_TARGET in dna_sequence:
            st.error("🚨 TARGET MUTATION DETECTED: Sickle Cell Positive (HbS)")
            
            # Find location to show the user
            idx = dna_sequence.find(HBS_TARGET)
            st.warning(f"**Mutation Found at Position:** base {idx+1} to {idx+11}")
            st.markdown(f"Sequence Context: `...{dna_sequence[max(0, idx-10):idx]}**{HBS_TARGET}**{dna_sequence[idx+11:idx+21]}...`")
            st.info("**Clinical Note:** Identified the $A \\rightarrow T$ point mutation at codon 6, causing a Glutamic Acid to Valine substitution.")
            
        elif HBA_TARGET in dna_sequence:
            st.success("✅ NO MUTATION FOUND: Normal Beta-Globin (HbA)")
            idx = dna_sequence.find(HBA_TARGET)
            st.markdown(f"Sequence Context: `...{dna_sequence[max(0, idx-10):idx]}**{HBA_TARGET}**{dna_sequence[idx+11:idx+21]}...`")
            st.info("**Clinical Note:** Sequence matches wild-type normal HBB sequence (GAG at codon 6).")
            
        else:
            st.warning("⚠️ INCONCLUSIVE REFERENCE MATCH")
            st.info("**Clinical Note:** The sequence could not be confidently mapped to the HBB Codon 6 region. Please verify if the uploaded FASTA contains the Beta-Globin gene region.")
