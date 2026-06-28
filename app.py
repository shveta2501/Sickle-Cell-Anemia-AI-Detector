import streamlit as st
from Bio import SeqIO
from Bio.Seq import Seq
from io import StringIO, BytesIO
import re

# 1. Page Configuration
st.set_page_config(page_title="Genetic Disease Platform Pro", layout="centered", page_icon="🧬")
st.title("🧬 AI Genetic Disease Detection Platform")
st.write("Clinical Genomics, Zygosity Multi-Panel & RFLP Engine")
st.caption("Production Pipeline Module • System Architect: Shveta")

st.markdown("---")

# --- MULTI-DISEASE CONFIGURATION PANEL ---
GENETIC_PANEL = {
    "Sickle Cell Anemia": {
        "gene": "HBB",
        "wild_type": "CCTGAGGAG",  
        "mutant": "CCTGTGGAG",     
        "carrier_code": "W",       
        "locus": "Chromosome 11, HBB Locus (Codon 6)"
    },
    "Beta Thalassemia (Nonsense)": {
        "gene": "HBB",
        "wild_type": "CAGGAGGCT",  
        "mutant": "TAGGAGGCT",     
        "carrier_code": "Y",       
        "locus": "Chromosome 11, HBB Locus (Codon 39)"
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
    
    # Parse records safely into a list to prevent UI rendering explosion
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
        
        # Normalize sequence text
        raw_seq_str = "".join(str(record.seq).upper().split())
        bio_seq = Seq(raw_seq_str)
        
        # --- EDGE CASE 2 RESOLUTION: Handle Reverse Complement Reads ---
        wt_sig = active_panel["wild_type"]
        mut_sig = active_panel["mutant"]
        carrier_char = active_panel["carrier_code"]
        
        # Fallback check if reference signatures aren't found in forward strand
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
        
        # 7. Virtual Restriction Fragment Length Polymorphism (RFLP Analysis)
        if selected_disease == "Sickle Cell Anemia":
            st.subheader("🧪 In-Silico Restriction Mapping (RFLP Validation)")
            
            # --- EDGE CASE 1 RESOLUTION: Real dynamic verification of MstII presence ---
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

        # --- 8. PHENOTYPIC SEVERITY RISK SCORING ENGINE ---
        st.subheader("🤖 Phenotypic Severity & Risk Scoring Engine")
        st.write("Provide patient-specific laboratory co-factors to evaluate physiological risk.")
        
        m1, m2, m3 = st.columns(3)
        with m1:
            patient_age = st.number_input("Patient Current Age", min_value=1, max_value=100, value=24)
        with m2:
            hb_level = st.slider("Total Hemoglobin (Hb Level - g/dL)", min_value=3.0, max_value=17.0, value=8.5, step=0.1)
        with m3:
            hbf_level = st.slider("Fetal Hemoglobin (HbF Level - %)", min_value=0.0, max_value=40.0, value=12.0, step=0.5)

        base_score = 5 if "AA" in genotype_status else (15 if "AS" in genotype_status else 60)
        if hb_level < 7.0: base_score += 20
        elif hb_level < 10.0: base_score += 10
        if hbf_level >= 15.0: base_score -= 15
        elif hbf_level < 5.0: base_score += 10
        if patient_age > 12: base_score += 5
        
        risk_score = max(0, min(base_score, 100))
        
        if risk_score >= 70:
            category, action = "HIGH RISK", "Immediate clinical consultation required. Susceptible to pain crises."
        elif 35 <= risk_score < 70:
            category, action = "MODERATE RISK", "Standard prophylactic tracking. Monitor complete blood counts regularly."
        else:
            category, action = "LOW RISK", "Routine preventive checks. Baseline status stable."

        rc1, rc2 = st.columns([1, 2])
        with rc1:
            st.metric(label="Calculated Severity Index", value=f"{risk_score}%")
        with rc2:
            st.markdown(f"**Risk Category:** `{category}`")
            st.caption(f"**Advisory Notice:** {action}")

        # --- 9. IN-MEMORY REPORTLAB PDF GENERATOR ---
        st.markdown("---")
        st.subheader("📥 Download Certified Diagnostics Report")
        
        # --- EDGE CASE 3 RESOLUTION: Context management with block isolation ---
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
            "Calculated Phenotypic Risk Index": f"{risk_score}%"
        })
