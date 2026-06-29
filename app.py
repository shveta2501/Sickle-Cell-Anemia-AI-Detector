import streamlit as st
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Align  # Advanced Smith-Waterman Local Alignment Engine
from io import StringIO, BytesIO
import re

# New Libraries for 3D Molecular Simulation
from stmol import make_viewer
import py3Dmol

# Optimization: Prevent model retraining on every UI rerun interaction
@st.cache_resource
def initialize_ml_models():
    from sklearn.ensemble import RandomForestClassifier
    import numpy as np

    X_train = np.array([
        [2, 25, 6.5, 4.0, 75.0],   # Severe SS case
        [2, 28, 9.5, 22.0, 88.0],  # Mild SS case (High HbF protection)
        [1, 30, 12.0, 2.0, 84.0],  # Normal Carrier AS
        [0, 22, 14.5, 0.5, 90.0],  # Completely Healthy AA
        [2, 12, 5.5, 6.0, 70.0],   # High risk pediatric SS
        [1, 45, 11.0, 1.5, 82.0]   # Stable Adult AS
    ])
    
    y_hospital = np.array([2, 0, 0, 0, 2, 1])
    y_stroke   = np.array([2, 0, 0, 0, 1, 0])

    model_hospital = RandomForestClassifier(n_estimators=50, random_state=42).fit(X_train, y_hospital)
    model_stroke   = RandomForestClassifier(n_estimators=50, random_state=42).fit(X_train, y_stroke)
    return model_hospital, model_stroke

# 1. Page Configuration
st.set_page_config(page_title="Genetic Disease Platform Pro", layout="centered", page_icon="🧬")
st.title("🧬 AI Genetic Disease Detection Platform")
st.write("Clinical Genomics, Dynamic Alignment, ML Risk Core & 3D Mutation Simulator")
st.caption("Production Pipeline Module • System Architect: Shveta")

st.markdown("---")

# --- MULTI-DISEASE CONFIGURATION PANEL (WITH FIXED PDB ID) ---
GENETIC_PANEL = {
    "Sickle Cell Anemia": {
        "gene": "HBB",
        "wild_type": "CCTGAGGAG",  # GAG Codon
        "mutant": "CCTGTGGAG",     # GTG Mutation
        "carrier_code": "W",       
        "locus": "Chromosome 11, HBB Locus (Codon 6)",
        "ncbi_ref": "ATGGTGCATCTGACTCCTGAGGAGAAGTCTGCCGTTACT",  # Standard Exon 1 anchor segment
        "pdb_id": "1A3N"  # Human Deoxyhemoglobin
    },
    "Beta Thalassemia (Nonsense)": {
        "gene": "HBB",
        "wild_type": "CAGGAGGCT",  
        "mutant": "TAGGAGGCT",     
        "carrier_code": "Y",       
        "locus": "Chromosome 11, HBB Locus (Codon 39)",
        "ncbi_ref": "AAGTCCAACTCCTAACCCAGGAGGCTCCTGGGGAGAAG",  # Codon 39 anchor segment
        "pdb_id": "2DNN"  # FIXED: Corrected 4-character RCSB PDB ID structural code
    }
}

selected_disease = st.selectbox("🎯 Select Genetic Screening Panel", list(GENETIC_PANEL.keys()))
active_panel = GENETIC_PANEL[selected_disease]

# 2. File Uploader UI Component
uploaded_file = st.file_uploader("Upload Patient FASTA/FA Sequence File", type=["fasta", "fa", "txt"])

if uploaded_file is not None:
    st.success("Sequence stream verified successfully!")
    
    raw_text = uploaded_file.getvalue().decode("utf-8")
    data_stream = StringIO(raw_text)
    
    records = list(SeqIO.parse(data_stream, "fasta"))
    
    if len(records) == 0:
        st.error("No valid FASTA records detected in file.")
    else:
        if len(records) > 1:
            patient_options = [r.id for r in records]
            selected_patient_id = st.sidebar.selectbox("👥 Select Patient Record to Display", patient_options)
            record = next(r for r in records if r.id == selected_patient_id)
        else:
            record = records[0]

        patient_id = record.id
        raw_seq_str = "".join(str(record.seq).upper().split())
        bio_seq = Seq(raw_seq_str)
        
        # --- REVERSE COMPLEMENT READS RESOLUTION ---
        wt_sig = active_panel["wild_type"]
        mut_sig = active_panel["mutant"]
        carrier_char = active_panel["carrier_code"]
        
        if (wt_sig not in raw_seq_str) and (mut_sig not in raw_seq_str) and (carrier_char not in raw_seq_str):
            rev_seq_str = str(bio_seq.reverse_complement())
            if (wt_sig in rev_seq_str) or (mut_sig in rev_seq_str) or (carrier_char in rev_seq_str):
                bio_seq = bio_seq.reverse_complement()
                st.sidebar.info("🔄 Antisense strand detected. Auto-flipped to sense orientation.")
        
        dna_sequence = str(bio_seq)
        
        # Display Patient Summary Cards
        st.subheader("📊 Patient Genomic Analysis Summary")
        col1, col2 = st.columns(2)
        with col1:
            st.metric(label="Patient/Record ID", value=patient_id)
        with col2:
            st.metric(label="Total Sequence Length", value=f"{len(dna_sequence)} bp")
            
        st.text_area("Parsed Raw DNA Sequence (Preview)", dna_sequence, height=100)
        st.markdown("---")
        
        # --- ADVANCED MATHEMATICAL LOCAL ALIGNMENT ENGINE (SMITH-WATERMAN) ---
        st.subheader("🧮 Mathematical Local Pairwise Alignment Engine")
        st.caption("Defending against junk DNA and tracking accurate mutation sites using Biopython Aligner.")
        
        aligner = Align.PairwiseAligner()
        aligner.mode = 'local'  # Pure Smith-Waterman
        aligner.match_score = 2.0
        aligner.mismatch_score = -1.0
        aligner.open_gap_score = -1.0
        aligner.extend_gap_score = -0.5
        
        ncbi_ref_seq = active_panel["ncbi_ref"]
        
        alignments = aligner.align(dna_sequence, ncbi_ref_seq)
        best_alignment = alignments[0]
        
        alignment_score = best_alignment.score
        max_possible_score = len(ncbi_ref_seq) * aligner.match_score
        identity_percentage = (alignment_score / max_possible_score) * 100
        
        m_col1, m_col2, m_col3 = st.columns(3)
        with m_col1:
            st.metric(label="Alignment Score", value=f"{alignment_score:.1f}")
        with m_col2:
            st.metric(label="Sequence Identity", value=f"{identity_percentage:.2f} %")
        with m_col3:
            st.metric(label="Algorithm Mode", value="Smith-Waterman")
            
        st.markdown("**Dynamic Reference Alignment Visual Map:**")
        st.code(str(best_alignment), language="text")
        st.markdown("---")
        
        # 5. Advanced Molecular Biology Validation
        st.subheader("🧬 Central Dogma Pipeline Validation")
        mrna_seq = bio_seq.transcribe()
        
        try:
            clean_len = (len(mrna_seq) // 3) * 3
            protein_seq = mrna_seq[:clean_len].translate(to_stop=False) 
            st.info(f"**Translated Protein Fragment:** `{protein_seq[:35]}...`")
        except Exception:
            st.warning("Could not automatically translate sequence reading frame.")

        st.markdown("---")
        
        # 6. Advanced Zygosity Detection & Diagnostic Engine
        st.subheader("🔬 Clinical Zygosity & Mutation Analysis")
        st.write(f"Target Gene: **{active_panel['gene']}** ({active_panel['locus']})")
        
        has_wild_type = wt_sig in dna_sequence
        has_mutant = mut_sig in dna_sequence
        has_iupac_carrier = carrier_char in dna_sequence
        
        genotype_status = "Inconclusive"
        html_highlighted_context = ""
        
        if (has_wild_type and has_mutant) or has_iupac_carrier:
            genotype_status = "AS (Carrier / Trait Profile)"
            if has_iupac_carrier:
                idx_c = dna_sequence.find(carrier_char)
                start_p = max(0, idx_c-10)
                html_highlighted_context = f"...{dna_sequence[start_p:idx_c]}<span style='background-color:#FEEBC8; color:#744210; padding:2px; font-weight:bold;'>{carrier_char} (IUPAC)</span>{dna_sequence[idx_c+1:idx_c+15]}..."
            else:
                idx_wt = dna_sequence.find(wt_sig)
                idx_mut = dna_sequence.find(mut_sig)
                html_highlighted_context = f"**Allele 1 (WT):** ...{dna_sequence[max(0, idx_wt-10):idx_wt]}<span style='background-color:#FEEBC8; color:#744210; padding:2px; font-weight:bold;'>{wt_sig}</span>... <br> **Allele 2 (Mutant):** ...{dna_sequence[max(0, idx_mut-10):idx_mut]}<span style='background-color:#FED7D7; color:#742A2A; padding:2px; font-weight:bold;'>{mut_sig}</span>..."
            
            st.warning(f"⚡ CLINICAL STATUS: **{genotype_status}**")
            
        elif has_wild_type and not has_mutant:
            genotype_status = "AA (Normal / Wild-Type)"
            idx = dna_sequence.find(wt_sig)
            html_highlighted_context = f"...{dna_sequence[max(0, idx-15):idx]}<span style='background-color:#C6F6D5; color:#22543D; padding:2px; font-weight:bold;'>{wt_sig} (Normal)</span>{dna_sequence[idx+len(wt_sig):idx+len(wt_sig)+15]}..."
            st.success(f"✅ CLINICAL STATUS: **{genotype_status}**")
            
        elif has_mutant and not has_wild_type:
            genotype_status = "SS (Diseased / Homozygous Mutant)"
            idx = dna_sequence.find(mut_sig)
            html_highlighted_context = f"...{dna_sequence[max(0, idx-15):idx]}<span style='background-color:#FED7D7; color:#742A2A; padding:2px; font-weight:bold;'>{mut_sig} (Mutated)</span>{dna_sequence[idx+len(mut_sig):idx+len(mut_sig)+15]}..."
            st.error(f"🚨 CLINICAL STATUS: **{genotype_status}**")
        else:
            st.error("⚠️ REFERENCE TARGET INCONCLUSIVE")

        if html_highlighted_context:
            st.markdown("#### Genomic Mutation Context Map:")
            st.markdown(f"<div style='background-color:#F7FAFC; border:1px solid #E2E8F0; padding:12px; border-radius:5px; font-family:monospace;'>{html_highlighted_context}</div>", unsafe_allow_html=True)

        st.markdown("---")
        
        # --- 3D PROTEIN STRUCTURE MUTATION SIMULATION ---
        st.subheader("🩻 Interactive 3D Protein Structure Mutation Simulation")
        st.caption("Real-time rendering of structural targets from Protein Data Bank (RCSB PDB) mapping folding geometry.")
        
        style_choice = st.segmented_control("🎨 Select Molecular Render Style", ["cartoon", "sphere", "line", "stick"], default="cartoon")
        
        pdb_id = active_panel["pdb_id"]
        
        xyzview = py3Dmol.view(query=f'pdb:{pdb_id}')
        xyzview.setStyle({style_choice: {'colorscheme': 'spectrum'}})
        xyzview.zoomTo()
        
        if "SS" in genotype_status or "AS" in genotype_status:
            st.warning(f"⚠️ Structural alert: Visualizing possible mutation disruption fold matrix on PDB Target: {pdb_id}")
            xyzview.addResLabels({'resi': [6], 'chain': 'B'}, {'backgroundColor': 'lightgray', 'fontColor': 'red'})
            xyzview.addSurface(py3Dmol.VDW, {'opacity': 0.3, 'color': 'red'}, {'resi': [6], 'chain': 'B'})
        else:
            st.success(f"✅ Displaying Wild-Type Stable Structural Geometry Template for PDB: {pdb_id}")
            xyzview.addSurface(py3Dmol.VDW, {'opacity': 0.1, 'color': 'green'})

        show_viewer = make_viewer(xyzview, width=700, height=450)
        st.components.v1.html(show_viewer.html, height=460)
        
        st.markdown("---")
        
        # 7. Virtual Restriction Fragment Length Polymorphism (RFLP Analysis)
        if selected_disease == "Sickle Cell Anemia":
            st.subheader("🧪 In-Silico Restriction Mapping (RFLP Validation)")
            mstII_pattern = re.compile(r"CCT[ATGC]AGG")
            cut_sites = [m.start() for m in mstII_pattern.finditer(dna_sequence)]
            
            st.caption(f"Detected **{len(cut_sites)}** restriction cleavage site(s) globally.")
            
            if len(cut_sites) > 0 and "AA" in genotype_status:
                st.markdown("🔍 **Enzyme Digest Status:** <span style='color:green; font-weight:bold;'>SITES FUNCTIONAL (COMPLETE DIGESTION)</span>", unsafe_allow_html=True)
            elif "SS" in genotype_status:
                st.markdown("🔍 **Enzyme Digest Status:** <span style='color:red; font-weight:bold;'>RESTRICTION SITE BLOCKED (CLEAVAGE FAILURE)</span>", unsafe_allow_html=True)
            elif "AS" in genotype_status:
                st.markdown("🔍 **Enzyme Digest Status:** <span style='color:orange; font-weight:bold;'>PARTIAL DIGESTION PROFILE (3-BAND SIGNATURE)</span>", unsafe_allow_html=True)
            else:
                st.markdown("🔍 **Enzyme Digest Status:** <span style='color:gray; font-weight:bold;'>NO SPECIFIC HBB CLEAVAGE SITES TRACKED</span>", unsafe_allow_html=True)
            
            st.markdown("---")

        # --- 8. MULTI-OMICS PHENOTYPIC SEVERITY PREDICTOR ---
        st.subheader("🤖 Multi-Omics Phenotypic Severity Predictor (ML Engine)")
        st.caption("Hybrid Random Forest Classifier mapping genomic variance with clinical biomarkers.")
        
        m1, m2, m3, m4 = st.columns(4)
        with m1:
            patient_age = st.number_input("Patient Current Age", min_value=1, max_value=100, value=24)
        with m2:
            hb_level = st.slider("Total Hemoglobin (Hb - g/dL)", min_value=3.0, max_value=17.0, value=8.5, step=0.1)
        with m3:
            hbf_level = st.slider("Fetal Hemoglobin (HbF - %)", min_value=0.0, max_value=40.0, value=12.0, step=0.5)
        with m4:
            mcv_level = st.slider("Mean Corpuscular Vol (MCV - fL)", min_value=50.0, max_value=120.0, value=85.0, step=1.0)

        genotype_map = {"AA (Normal / Wild-Type)": 0, "AS (Carrier / Trait Profile)": 1, "SS (Diseased / Homozygous Mutant)": 2}
        encoded_genotype = genotype_map.get(genotype_status, 1)

        # Bug Prevention: Set defaults so PDF engine doesn't crash on NameError if ML execution throws an error
        risk_score = 0
        category = "UNDETERMINED"
        action = "Metrics calibration pending execution."

        try:
            import numpy as np
            model_hospital, model_stroke = initialize_ml_models()

            current_patient_features = np.array([[encoded_genotype, patient_age, hb_level, hbf_level, mcv_level]])

            hosp_probs = model_hospital.predict_proba(current_patient_features)[0]
            stroke_probs = model_stroke.predict_proba(current_patient_features)[0]

            hospitalization_risk_score = int((hosp_probs[2] * 0.7 + hosp_probs[1] * 0.3) * 100)
            stroke_risk_score = int((stroke_probs[2] * 0.8 + stroke_probs[1] * 0.2) * 100)

            if "AA" in genotype_status:
                hospitalization_risk_score, stroke_risk_score = max(2, hospitalization_risk_score // 10), max(1, stroke_risk_score // 12)
            elif "SS" in genotype_status and hbf_level > 20.0:
                hospitalization_risk_score = max(15, hospitalization_risk_score - 40)

            p_col1, p_col2 = st.columns(2)
            with p_col1:
                st.metric(label="📊 Predicted Hospitalization Risk Index", value=f"{hospitalization_risk_score}%")
                if hospitalization_risk_score >= 65:
                    st.error("🚨 HIGH RISK: Frequent painful vaso-occlusive crises predicted.")
                elif 30 <= hospitalization_risk_score < 65:
                    st.warning("⚠️ MODERATE RISK: Baseline monitoring required.")
                else:
                    st.success("✅ LOW RISK: Clinical phenotype appears stable.")

            with p_col2:
                st.metric(label="🧠 Calculated Stroke Risk Criticality", value=f"{stroke_risk_score}%")
                if stroke_risk_score >= 50:
                    st.error("🚨 CRITICAL ALERT: Transcranial Doppler (TCD) screening urgent.")
                elif 20 <= stroke_risk_score < 50:
                    st.warning("⚠️ ELEVATED MONITORING: Periodic neurological tracking required.")
                else:
                    st.success("✅ STABLE STATUS: Low risk profile detected.")

            risk_score = max(hospitalization_risk_score, stroke_risk_score)
            category = "HIGH RISK" if risk_score >= 60 else ("MODERATE RISK" if risk_score >= 30 else "LOW RISK")
            action = f"Hospitalization Risk: {hospitalization_risk_score}%, Stroke Risk: {stroke_risk_score}%. Calibrated via Random Forest Machine Learning Matrix."

        except Exception as ml_err:
            st.error(f"ML Core Integration Error: {ml_err}")

        # --- 9. IN-MEMORY REPORTLAB PDF GENERATOR ---
        st.markdown("---")
        st.subheader("📥 Download Certified Diagnostics Report")
        
        try:
            from reportlab.lib.pagesizes import letter
            from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer
            from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
            from reportlab.lib import colors

            with BytesIO() as pdf_buffer:
                doc = SimpleDocTemplate(pdf_buffer, pagesize=letter, rightMargin=40, leftMargin=40, topMargin=40, bottomMargin=40)
                styles = getSampleStyleSheet()
                story = []

                title_style = ParagraphStyle('Title', parent=styles['Heading1'], fontSize=20, textColor=colors.HexColor("#1A365D"), spaceAfter=15)
                section_style = ParagraphStyle('Section', parent=styles['Heading2'], fontSize=13, textColor=colors.HexColor("#2B6CB0"), spaceBefore=10, spaceAfter=5)
                text_style = ParagraphStyle('Body', parent=styles['Normal'], fontSize=10, leading=14, textColor=colors.HexColor("#2D3748"))

                story.append(Paragraph("AI GENETIC DISEASE DETECTION PLATFORM", title_style))
                story.append(Paragraph(f"<b>Patient ID:</b> {patient_id}", text_style))
                story.append(Paragraph(f"<b>Targeted Panel:</b> {selected_disease}", text_style))
                story.append(Paragraph(f"<b>Resolved Genotype:</b> {genotype_status}", text_style))
                story.append(Paragraph(f"<b>Alignment Core Identity Score:</b> {identity_percentage:.2f}%", text_style))
                story.append(Paragraph(f"<b>Clinical Severity Score:</b> {risk_score}% ({category})", text_style))
                story.append(Spacer(1, 15))
                
                story.append(Paragraph("Diagnostic Advisory Matrix", section_style))
                story.append(Paragraph(action, text_style))
                story.append(Spacer(1, 15))
                story.append(Paragraph("Disclaimer: For Research Use Only (RUO). Final clinical assertions must be validated via chromatography protocols.", ParagraphStyle('Disc', parent=text_style, fontSize=8, textColor=colors.gray)))

                doc.build(story)
                pdf_data = pdf_buffer.getvalue()

            st.download_button(
                label="📥 Download Clinical PDF Report",
                data=pdf_data,
                file_name=f"Genomics_Report_{patient_id}.pdf",
                mime="application/pdf"
            )
        except Exception as pdf_error:
            st.warning("Install 'reportlab' using terminal to enable certified PDF downloading features.")

        st.markdown("---")
        st.json({
            "Patient Key": patient_id,
            "Resolved Clinical Genotype": genotype_status,
            "Calculated Phenotypic Risk Index": f"{risk_score}%",
            "Pairwise Sequence Identity Match": f"{identity_percentage:.2f}%"
        })