# PGx Report Engine

**Advanced Pharmacogenomics Report Generator**  
Built for Illumina GSA array and standard VCF files with South Asian focus.

### Features
- Analyzes **11 clinically important pharmacogenomic genes**
- Follows **CPIC guidelines** and FDA pharmacogenetic labels
- Detects key variants for:
  - CYP2C19 (Clopidogrel, antidepressants, PPIs)
  - Warfarin sensitivity (CYP2C9 + VKORC1)
  - SLCO1B1 (Statin-induced myopathy)
  - CYP2D6 and others
- Includes South Asian allele frequency context (especially CYP2C19*2)
- Clean, interactive Streamlit dashboard with tabbed interface

### Project Structure
pgx-report-engine/
├── pgx_engine.py           # Backend: all analysis functions
├── pgx_streamlit_app.py    # Streamlit frontend
├── requirements.txt
└── README.md
text### How to Run Locally

```bash
# 1. Clone the repository
git clone <your-repo-url>
cd pgx-report-engine

# 2. Install dependencies
pip install -r requirements.txt

# 3. Run the app
streamlit run pgx_streamlit_app.py
Deployment
This app is ready to deploy on Streamlit Community Cloud.

Push the code to GitHub
Go to Streamlit Cloud
Connect your GitHub repository
Select pgx_streamlit_app.py as the main file
Deploy

Important Disclaimer
This tool is developed for educational and research purposes only.
It has limited variant coverage due to the nature of Illumina GSA genotyping arrays.
Genes marked as "(assumed)" could not be fully assessed due to missing variants on the array.
Always consult a qualified clinician or certified pharmacogenomics specialist before making any medication decisions based on this report.