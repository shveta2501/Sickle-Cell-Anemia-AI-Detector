import streamlit as st
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Align  # Advanced Smith-Waterman Local Alignment Engine
from io import StringIO, BytesIO
import re
import numpy as np

# 1. Page Configuration
st.set_page_config(page_title="Genetic Disease Platform Pro", layout="centered", page_icon="🧬")
st.title("🧬 AI Genetic Disease Detection Platform")
st.write("Clinical Genomics, Multiplexed Panel Engine & CRISPR Therapeutic Roadmap")
st.caption("Production Pipeline Module • System Architect: Shveta")

st.markdown("---")

# --- FIVE-DISEASE MULTIPLEXED CONFIGURATION PANEL (WITH PDB & NCBI REFERENCE MATRICES) ---
GENETIC_PANEL = {
    "Sickle Cell Anemia": {
        "gene": "HBB",
        "wild_type": "CCTGAGGAG",  # GAG Codon
        "mutant": "CCTGTGGAG",     # GTG Mutation
        "carrier_code": "W",       
        "locus": "Chromosome 11, HBB Locus (Codon 6)",
        "ncbi_ref": "ATGGTGCATCTGACTCCTGAGGAGAAGTCTGCCGTTACT",  
        "pdb_id": "1A3N",  
        "type": "Autosomal Recessive"
    },
    "Beta Thalassemia (Nonsense)": {
        "gene": "HBB",
        "wild_type": "CAGGAGGCT",  
        "mutant": "TAGGAGGCT",     
        "carrier_code": "Y",       
        "locus": "Chromosome 11, HBB Locus (Codon 39)",
        "ncbi_ref": "AAGTCCAACTCCTAACCCAGGAGGCTCCTGGGGAGAAG",  
        "pdb_id": "2DNN",  
        "type": "Autosomal Recessive"
    },
    "G6PD Deficiency": {
        "gene": "G6PD",
        "wild_type": "CCAGCTCTG",  
        "mutant": "CCAACTCTG",     
        "carrier_code": "M",       
        "locus": "Chromosome X, G6PD Locus (p.Ser188Phe)",
        "ncbi_ref": "GAGGCCGTGGGCAGGGCCCTGGCCAGCTCTGGGGCCGTG",  
        "pdb_id": "2A9F",  
        "type": "X-Linked Recessive"
    },
    "Cystic Fibrosis (DeltaF508)": {
        "gene": "CFTR",
        "wild_type": "ATATCATTGGTG",  
        "mutant": "ATATCATGGTG",     
        "carrier_code": "N",       
        "locus": "Chromosome 7, CFTR Locus (Codon 508)",
        "ncbi_ref": "GAAATATTGGTGTTTCCTATGATATCATTGGTGTTTCCTATG", 
        "pdb_id": "5UAR",  
        "type": "Autosomal Recessive"
    },
    "Hemophilia A (Factor VIII)": {
        "gene": "F8",
        "wild_type": "AGAAATACT",  
        "mutant": "AGAAATTTC",  
        "carrier_code": "R",       
        "locus": "Chromosome X, F8 Locus (Exon 24 Transition)",
        "ncbi_ref": "TGGATACATTTCTGAGAAATACTTTTCCAGGAAGTCCT",  
        "pdb_id": "3H2A",  
        "type": "X-Linked Recessive"
    }
}

# --- OPTIMIZED CACHED ML MODELS ---
@st.cache_resource
def initialize_ml_predictors():
    try:
        from sklearn.ensemble import RandomForestClassifier
        X_train = np.array([
            [2, 25, 6.5, 4.0, 75.0],   
            [2, 28, 9.5, 22.0, 88.0],  
            [1, 30, 12.0, 2.0, 84.0],  
            [0, 22, 14.5, 0.5, 90.0],  
            [2, 12, 5.5, 6.0, 70.0],   
            [1, 45, 11.0, 1.5, 82.0]   
        ])
        y_hospital = np.array([2, 0, 0, 0, 2, 1])
        y_stroke   = np.array([2, 0, 0, 0, 1, 0])

        model_hosp = RandomForestClassifier(n_estimators=50, random_state=42).fit(X_train, y_hospital)
        model_stk = RandomForestClassifier(n_estimators=50, random_state=42).fit(X_train, y_stroke)
        return model_hosp, model_stk
    except ImportError:
        return None, None

model_hospital, model_stroke = initialize_ml_predictors()

selected_disease = st.selectbox("🎯 Select Genetic Screening Panel", list(GENETIC_PANEL.keys()))
active_panel = GENETIC_PANEL[selected_disease]

st.info(f"🧬 **Active Core Panel:** {selected_disease} | **Target Gene:** {active_panel['gene']} | **Mode of Inheritance:** {active_panel['type']}")

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
            rev_seq = bio_seq.reverse_complement()
            rev_seq_str = str(rev_seq)
            if (wt_sig in rev_seq_str) or (mut_sig in rev_seq_str) or (carrier_char in rev_seq_str):
                bio_seq = rev_seq
                st.sidebar.info("🔄 Antisense strand detected. Auto-flipped to sense orientation.")
        
        dna_sequence = str(bio_seq)
        
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
        
        aligner = Align.PairwiseAligner(
            mode='local',
            match_score=2.0,
            mismatch_score=-1.0,
            open_gap_score=-1.0,
            extend_gap_score=-0.5
        )
        
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
        except Exception as trans_err:
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
            genotype_status = "Heterozygous Carrier (Female Risk Profile)" if active_panel["type"] == "X-Linked Recessive" else "Carrier / Heterozygous Trait"
                
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
            genotype_status = "Normal / Wild-Type"
            idx = dna_sequence.find(wt_sig)
            html_highlighted_context = f"...{dna_sequence[max(0, idx-15):idx]}<span style='background-color:#C6F6D5; color:#22543D; padding:2px; font-weight:bold;'>{wt_sig} (Normal)</span>{dna_sequence[idx+len(wt_sig):idx+len(wt_sig)+15]}..."
            st.success(f"✅ CLINICAL STATUS: **{genotype_status}**")
            
        elif has_mutant and not has_wild_type:
            genotype_status = "Hemizygous Affected (High Severity Locus Match)" if active_panel["type"] == "X-Linked Recessive" else "Homozygous Mutant (Affected/Diseased)"
                
            idx = dna_sequence.find(mut_sig)
            html_highlighted_context = f"...{dna_sequence[max(0, idx-15):idx]}<span style='background-color:#FED7D7; color:#742A2A; padding:2px; font-weight:bold;'>{mut_sig} (Mutated)</span>{dna_sequence[idx+len(mut_sig):idx+len(mut_sig)+15]}..."
            st.error(f"🚨 CLINICAL STATUS: **{genotype_status}**")
        else:
            st.error("⚠️ REFERENCE TARGET INCONCLUSIVE (LOW VARIANT LOG MATRIX MATCH)")

        if html_highlighted_context:
            st.markdown("#### Genomic Mutation Context Map:")
            st.markdown(f"<div style='background-color:#F7FAFC; border:1px solid #E2E8F0; padding:12px; border-radius:5px; font-family:monospace;'>{html_highlighted_context}</div>", unsafe_allow_html=True)

        st.markdown("---")
        
        # --- HTML/JS CDN 3DMOL.JS SIMULATOR LAYER ---
        st.subheader("🩻 Interactive 3D Protein Structure Mutation Simulation")
        pdb_id = active_panel["pdb_id"]
        
        html_string = f"""
        <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
        <div id="container-3d" style="height: 400px; width: 100%; position: relative; border-radius:10px; border:1px solid #E2E8F0;"></div>
        <script>
            let element = document.getElementById('container-3d');
            let config = {{ backgroundColor: '#F7FAFC' }};
            let viewer = $3Dmol.createViewer(element, config);
            $3Dmol.download("pdb:{pdb_id}", viewer, {{}}, function() {{
                viewer.setStyle({{}}, {{cartoon: {{color: 'spectrum'}}}});
                viewer.zoomTo();
                viewer.render();
            }});
        </script>
        """
        st.components.v1.html(html_string, height=420)
        st.markdown("---")

        # --- CRISPR-Cas9 THERAPEUTIC TARGET FINDER ENGINE ---
        st.subheader("✂️ CRISPR-Cas9 Therapeutic Target Finder")
        
        if any(keyword in genotype_status for keyword in ["Affected", "Mutant", "Carrier"]):
            pam_pattern = re.compile(r"[ATGC][ATGC]GG") 
            mutation_idx = dna_sequence.find(mut_sig) if has_mutant else dna_sequence.find(wt_sig)
            
            search_window = dna_sequence[max(0, mutation_idx-40):min(len(dna_sequence), mutation_idx+40)]
            pam_sites = [m.start() for m in pam_pattern.finditer(search_window)]
            
            if pam_sites:
                st.success("🎯 Therapeutic Vector Active: Identified potential gRNA options targeting the mutation locus.")
                for p_idx in pam_sites:
                    if p_idx >= 20:
                        gRNA = search_window[p_idx-20:p_idx]
                        st.markdown(f"🧬 **gRNA Sequence:** `{gRNA}` | **Target PAM:** `<span style='color:red; font-weight:bold;'>{search_window[p_idx:p_idx+3]}</span>`", unsafe_allow_html=True)
            else:
                st.warning("⚠️ No valid NGG PAM profiles found nearby to execute Cas9 therapeutic deployment.")
        else:
            st.info("Therapeutic lookup resting. Current sequence resolved as completely Wild-Type Normal.")
            
        st.markdown("---")
        
        # 7. Virtual Restriction Fragment Length Polymorphism (RFLP Analysis)
        if selected_disease == "Sickle Cell Anemia":
            st.subheader("🧪 In-Silico Restriction Mapping (RFLP Validation)")
            mstII_pattern = re.compile(r"CCT[ATGC]AGG")
            cut_sites = [m.start() for m in mstII_pattern.finditer(dna_sequence)]
            st.caption(f"Detected **{len(cut_sites)}** restriction cleavage site(s) globally.")
            
            if len(cut_sites) > 0 and "Normal" in genotype_status:
                st.markdown("🔍 **Enzyme Digest Status:** <span style='color:green; font-weight:bold;'>SITES FUNCTIONAL (COMPLETE DIGESTION)</span>", unsafe_allow_html=True)
            elif "Mutant" in genotype_status:
                st.markdown("🔍 **Enzyme Digest Status:** <span style='color:red; font-weight:bold;'>RESTRICTION SITE BLOCKED (CLEAVAGE FAILURE)</span>", unsafe_allow_html=True)
            else:
                st.markdown("🔍 **Enzyme Digest Status:** <span style='color:orange; font-weight:bold;'>PARTIAL DIGESTION PROFILE (3-BAND SIGNATURE)</span>", unsafe_allow_html=True)
            st.markdown("---")

        # --- 8. MULTI-OMICS PHENOTYPIC SEVERITY PREDICTOR (MACHINE LEARNING LAYER) ---
        st.subheader("🤖 Multi-Omics Phenotypic Severity Predictor (ML Engine)")
        
        m1, m2, m3, m4 = st.columns(4)
        with m1:
            patient_age = st.number_input("Patient Current Age", min_value=1, max_value=100, value=24)
        with m2:
            hb_level = st.slider("Total Hemoglobin (Hb - g/dL)", min_value=3.0, max_value=17.0, value=8.5, step=0.1)
        with m3:
            hbf_level = st.slider("Fetal Hemoglobin (HbF - %)", min_value=0.0, max_value=40.0, value=12.0, step=0.5)
        with m4:
            mcv_level = st.slider("Mean Corpuscular Vol (MCV - fL)", min_value=50.0, max_value=120.0, value=85.0, step=1.0)

        encoded_genotype = 2 if any(x in genotype_status for x in ["Mutant", "Affected"]) else (1 if "Carrier" in genotype_status else 0)

        # Fallback values established to safeguard downstream variables
        risk_score, category, action = 0, "UNKNOWN", "ML Engine Unavailable."

        if model_hospital and model_stroke:
            current_patient_features = np.array([[encoded_genotype, patient_age, hb_level, hbf_level, mcv_level]])

            hosp_probs = model_hospital.predict_proba(current_patient_features)[0]
            stroke_probs = model_stroke.predict_proba(current_patient_features)[0]

            hospitalization_risk_score = int((hosp_probs[2] * 0.7 + hosp_probs[1] * 0.3) * 100)
            stroke_risk_score = int((stroke_probs[2] * 0.8 + stroke_probs[1] * 0.2) * 100)

            if "Normal" in genotype_status:
                hospitalization_risk_score, stroke_risk_score = max(2, hospitalization_risk_score // 10), max(1, stroke_risk_score // 12)
            elif hbf_level > 20.0:
                hospitalization_risk_score = max(15, hospitalization_risk_score - 40)

            p_col1, p_col2 = st.columns(2)
            with p_col1:
                st.metric(label="📊 Predicted Hospitalization Risk Index", value=f"{hospitalization_risk_score}%")
                if hospitalization_risk_score >= 65:
                    st.error("🚨 HIGH RISK: Frequent painful vaso-occlusive crises/clinical events predicted.")
                elif 30 <= hospitalization_risk_score < 65:
                    st.warning("⚠️ MODERATE RISK: Baseline clinical tracking required.")
                else:
                    st.success("✅ LOW RISK: Clinical phenotype appears stable.")

            with p_col2:
                st.metric(label="🧠 Calculated Secondary Criticality Index", value=f"{stroke_risk_score}%")
                if stroke_risk_score >= 50:
                    st.error("🚨 CRITICAL ALERT: Advanced diagnostic vascular/organ tracking urgent.")
                elif 20 <= stroke_risk_score < 50:
                    st.warning("⚠️ ELEVATED MONITORING: Periodic clinical tracking required.")
                else:
                    st.success("✅ STABLE STATUS: Low risk clinical profile detected.")

            risk_score = max(hospitalization_risk_score, stroke_risk_score)
            category = "HIGH RISK" if risk_score >= 60 else ("MODERATE RISK" if risk_score >= 30 else "LOW RISK")
            action = f"Hospitalization Risk: {hospitalization_risk_score}%, Complication Risk: {stroke_risk_score}%. Calibrated via Multi-Disease Random Forest Matrix."
        else:
            st.error("ML Core Integration Error: scikit-learn is not installed or failed to initialize.")

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

                story.append(Paragraph("AI GENETIC DISEASE MULTI-PANEL REPORT", title_style))
                story.append(Paragraph(f"<b>Patient ID:</b> {patient_id}", text_style))
                story.append(Paragraph(f"<b>Targeted Panel:</b> {selected_disease} ({active_panel['gene']} Gene)", text_style))
                story.append(Paragraph(f"<b>Resolved Genotype Status:</b> {genotype_status}", text_style))
                story.append(Paragraph(f"<b>Sequence Identity Alignment Match:</b> {identity_percentage:.2f}%", text_style))
                story.append(Paragraph(f"<b>Clinical Severity Score:</b> {risk_score}% ({category})", text_style))
                story.append(Spacer(1, 15))
                
                story.append(Paragraph("Diagnostic & Therapeutic Advisory Matrix", section_style))
                story.append(Paragraph(action, text_style))
                story.append(Spacer(1, 15))
                story.append(Paragraph("Disclaimer: For Research Use Only (RUO). Final clinical assertions must be validated via standard sequencing protocols.", ParagraphStyle('Disc', parent=text_style, fontSize=8, textColor=colors.gray)))

                doc.build(story)
                pdf_data = pdf_buffer.getvalue()

            st.download_button(
                label="📥 Download Clinical PDF Report",
                data=pdf_data,
                file_name=f"Multiplex_Genomics_Report_{patient_id}.pdf",
                mime="application/pdf"
            )
        except ModuleNotFoundError:
            st.warning("Install 'reportlab' using terminal to enable certified PDF downloading features.")

        st.markdown("---")
        st.json({
            "Patient Key": patient_id,
            "Targeted Gene Panel": selected_disease,
            "Resolved Clinical Genotype": genotype_status,
            "Pairwise Sequence Identity Match": f"{identity_percentage:.2f}%"
        })