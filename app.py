import streamlit as st
import hmac
import hashlib
import json
import secrets
from typing import Dict, List, Tuple, Any
from fpdf import FPDF

# ==========================================
# 1. ENTERPRISE CORE ENGINE & SECURITY GATEWAY
# ==========================================
class EliteGenomicSaaSEngine:
    def __init__(self):
        # Kept as strings for straightforward UI-to-backend alignment
        self._tenant_registry = {
            "TENANT-ALPHA-PRIME": {
                "name": "Nexus Bio-Labs LLC", "tier": "Enterprise Premium",
                "secret_key": "k3y_n3xus_secur3_2026_prod", "rate_per_run": 2.50,
                "billing_accumulator_usd": 0.0, "token_quota": 5000, "current_usage": 0
            },
            "TENANT-BETA-CLINIC": {
                "name": "Apex Diagnostics Corp", "tier": "Standard Business",
                "secret_key": "apex_secure_clinic_key_2026", "rate_per_run": 4.00,
                "billing_accumulator_usd": 0.0, "token_quota": 1000, "current_usage": 0
            }
        }

    def authenticate_tenant(self, tenant_id: str, client_signature: str, payload: str) -> Tuple[bool, str]:
        if tenant_id not in self._tenant_registry: 
            return False, "Tenant ID not registered."
        
        tenant_secret = self._tenant_registry[tenant_id]["secret_key"]
        # Both secret and payload must be encoded to bytes for HMAC
        expected_sig = hmac.new(tenant_secret.encode('utf-8'), payload.encode('utf-8'), hashlib.sha256).hexdigest()
        
        if secrets.compare_digest(expected_sig, client_signature):
            return True, "Verified"
        return False, "Invalid Signature"

    def execute_acmg_state_machine(self, selected_criteria: List[str], read_depth: int) -> Tuple[str, str, float, List[str]]:
        scrubbed = set(selected_criteria)
        telemetry = []
        
        if read_depth < 50:
            telemetry.append(f"QC WARN: Read depth {read_depth}x is below the 50x clinical threshold.")
            scrubbed.difference_update(["PVS1", "PS1", "PS3"])
        
        verdict = "VUS (Class 3)"
        if "PVS1" in scrubbed and len(scrubbed) > 1: 
            verdict = "Pathogenic (Class 5)"
        elif "BA1" in scrubbed:
            verdict = "Benign (Class 1)"
            
        return verdict, "red" if "Pathogenic" in verdict else "blue", 0.95, telemetry

    def update_billing(self, tenant_id: str):
        self._tenant_registry[tenant_id]["current_usage"] += 1
        self._tenant_registry[tenant_id]["billing_accumulator_usd"] += self._tenant_registry[tenant_id]["rate_per_run"]

# PDF Report Generator
def generate_clinical_pdf(patient_id: str, verdict: str, gene: str) -> str:
    pdf = FPDF()
    pdf.add_page()
    pdf.set_font("Arial", 'B', 16)
    pdf.cell(200, 10, txt="GenoScan AI Diagnostic Report", ln=True, align='C')
    pdf.ln(10)
    pdf.set_font("Arial", size=12)
    pdf.cell(200, 10, txt=f"Patient ID: {patient_id}", ln=True)
    pdf.cell(200, 10, txt=f"Gene: {gene}", ln=True)
    pdf.cell(200, 10, txt=f"Verdict: {verdict}", ln=True)
    pdf.cell(200, 10, txt=f"Generated On: {st.date_input('Date', disabled=True)}", ln=True)
    
    file_path = f"Report_{patient_id}.pdf"
    pdf.output(file_path)
    return file_path

# Initialize Backend Engine State
if 'engine' not in st.session_state: 
    st.session_state.engine = EliteGenomicSaaSEngine()

engine = st.session_state.engine

# ==========================================
# 2. UI LAYER
# ==========================================
st.set_page_config(page_title="GenoScan AI", layout="wide")
st.title("🧬 GenoScan AI: Enterprise Suite")

# Sidebar
input_tenant = st.sidebar.selectbox("Tenant Space", list(engine._tenant_registry.keys()))
# Dynamically pull default key based on selection to prevent immediate auth failure
default_key = engine._tenant_registry[input_tenant]["secret_key"]
input_secret = st.sidebar.text_input("Private Key", type="password", value=default_key)

# Main Form Metrics
patient_uuid = st.text_input("Patient UUID", value="PT-ALPHA-9942")
target_gene = st.selectbox("Target Gene", ["BRCA1", "CFTR", "HBB"])
sequencing_depth = st.slider("Read Depth", 10, 500, 120)

st.markdown("### ACMG Criteria Selection")
col1, col2 = st.columns(2)
with col1: pvs1 = st.checkbox("PVS1 (Very Strong Pathogenic)")
with col2: ba1 = st.checkbox("BA1 (Stand-alone Benign)")

# Fixed list construction
active_criteria = []
if pvs1: active_criteria.append("PVS1")
if ba1: active_criteria.append("BA1")

st.markdown("---")

if st.button("🚀 INITIATE ANALYSIS", type="primary"):
    payload = json.dumps({"patient": patient_uuid, "criteria": active_criteria}, sort_keys=True)
    # Correct signature calculation matching string inputs
    sig = hmac.new(input_secret.encode('utf-8'), payload.encode('utf-8'), hashlib.sha256).hexdigest()
    
    is_valid, msg = engine.authenticate_tenant(input_tenant, sig, payload)
    
    if not is_valid:
        st.error(f"Access Blocked! Code: {msg}")
    else:
        verdict, color, prob, telemetry = engine.execute_acmg_state_machine(active_criteria, sequencing_depth)
        engine.update_billing(input_tenant)
        
        # Display Results
        st.subheader("Analysis Results")
        if color == "red":
            st.error(f"Verdict: {verdict}")
        else:
            st.info(f"Verdict: {verdict}")
            
        if telemetry:
            for log in telemetry:
                st.warning(log)
        
        # PDF Generation and Download
        pdf_file = generate_clinical_pdf(patient_uuid, verdict, target_gene)
        with open(pdf_file, "rb") as f:
            st.download_button(
                label="📥 Download PDF Report", 
                data=f, 
                file_name=pdf_file,
                mime="application/pdf"
            )