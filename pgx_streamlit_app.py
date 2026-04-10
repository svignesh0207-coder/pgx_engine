#!/usr/bin/env python
# coding: utf-8

# In[1]:


# pgx_streamlit_app.py
import streamlit as st
import pandas as pd
from pathlib import Path
from pgx_engine import extract_pgx_variants, clean_pgx_variants, call_pgx_alleles, generate_pgx_report

st.set_page_config(page_title="PGx Report Engine", layout="wide", page_icon="🧬")

st.title("🧬 PGx Report Engine")
st.subheader("Advanced Pharmacogenomics Analyzer | 11 Genes | CPIC Guidelines")
st.markdown("**South Asian Context Included** | Illumina GSA + Standard VCF Support")

# Sidebar
st.sidebar.header("Upload Genetic Data")
uploaded_file = st.sidebar.file_uploader("Upload VCF file", type=["vcf", "vcf.gz"])

sample_index = st.sidebar.number_input("Sample Index", min_value=0, value=0, step=1)

if st.sidebar.button("🚀 Generate Advanced PGx Report", type="primary"):
    if uploaded_file is None:
        st.error("Please upload a VCF file")
    else:
        with st.spinner("Analyzing VCF and generating comprehensive PGx report..."):
            # Save uploaded file
            temp_path = "temp_upload.vcf"
            with open(temp_path, "wb") as f:
                f.write(uploaded_file.getbuffer())

            try:
                report = generate_pgx_report(temp_path, sample_index)

                st.success("✅ Comprehensive PGx Report Generated Successfully!")

                # Tabs for Advanced UI
                tab1, tab2, tab3 = st.tabs(["📊 Variants Detected", "🧬 Gene Analysis", "⚠️ Summary & Disclaimer"])

                with tab1:
                    st.subheader("Detected PGx Variants")
                    st.dataframe(report["variants"], use_container_width=True)

                with tab2:
                    st.subheader("Detailed Gene Analysis")
                    for gene, data in report["results"].items():
                        with st.expander(f"🔬 {gene}", expanded=True):
                            col1, col2 = st.columns([1, 1])
                            with col1:
                                for key in ["diplotype", "phenotype", "genotype", "cyp2c9_diplotype", "vkorc1_genotype"]:
                                    if key in data:
                                        st.write(f"**{key.replace('_', ' ').title()}**: {data[key]}")
                            with col2:
                                st.write(f"**Drugs**: {data.get('drugs', 'N/A')}")
                                st.write(f"**Conditions**: {data.get('conditions', 'N/A')}")
                            
                            if data.get("risk"):
                                st.warning(f"**Risk**: {data['risk']}")
                            if data.get("note"):
                                st.info(f"**Note**: {data['note']}")
                            if data.get("sa_note"):
                                st.success(f"**South Asian Context**: {data['sa_note']}")

                with tab3:
                    st.subheader("Summary & Important Disclaimer")
                    st.info("""
                    This report is generated based on variants detectable on Illumina GSA array.
                    Genes marked as '(assumed)' have limited coverage.
                    All drug recommendations follow CPIC clinical guidelines.
                    """)
                    st.warning("""
                    **Disclaimer**: This tool is for educational and research purposes only. 
                    It is not a substitute for professional medical advice. 
                    Always consult a qualified clinician or certified PGx specialist before changing any medication.
                    """)

                # Download Section
                st.subheader("Download Report")
                col1, col2 = st.columns(2)
                with col1:
                    csv = report["variants"].to_csv(index=False)
                    st.download_button("📥 Download Variants (CSV)", csv, "pgx_variants.csv", "text/csv")

            except Exception as e:
                st.error(f"Error: {str(e)}")

else:
    st.info("👈 Upload a VCF file from the sidebar and click 'Generate Advanced PGx Report'")
    st.markdown("### Features Included:")
    st.markdown("- 11 pharmacogenomic genes\n- CPIC guideline-based recommendations\n- South Asian allele frequency context\n- Advanced tabbed interface\n- Clean professional report")

st.caption("PGx Report Engine | Professional Modular Version | Built by Vignesh Sivakumar")


# In[ ]:




