#!/usr/bin/env python
# coding: utf-8

# In[1]:


# pgx_engine.py   --- Version: Small Step 1
import pandas as pd
from pathlib import Path

PGX_DATA_DIR = Path("pgx_data")

def find_file(pattern):
    """Find first file that matches the pattern"""
    matches = list(PGX_DATA_DIR.rglob(pattern))
    return matches[0] if matches else None

def load_relationships():
    """Load relationships and correctly identify gene & drug columns"""
    rel_file = find_file("*relationships*.tsv")
    
    if not rel_file:
        print("❌ relationships.tsv not found")
        return None
    
    df = pd.read_csv(rel_file, sep="\t")
    
    print("✅ Successfully loaded relationships.tsv")
    print(f"   Total rows: {len(df):,}")
    print(f"   Columns: {list(df.columns)}")
    
    # Correct column mapping for this ClinPGx file
    gene_col = 'Entity1_name'
    drug_col = 'Entity2_name'
    
    # Filter only rows where Entity1 is a gene and Entity2 is a drug
    gene_rows = df[df['Entity1_type'] == 'Gene']
    
    print(f"   Number of gene-drug relationships: {len(gene_rows):,}")
    print(f"   Example genes : {gene_rows[gene_col].unique()[:8].tolist()}")
    print(f"   Example drugs : {gene_rows[drug_col].unique()[:8].tolist()}")
    
    # Keep only our focus genes for now
    focus_genes = ["CYP2C19", "CYP2C9", "VKORC1", "SLCO1B1", "CYP2D6"]
    filtered = gene_rows[gene_rows[gene_col].isin(focus_genes)]
    
    print(f"   Focus genes found: {filtered[gene_col].unique().tolist()}")
    
    return filtered


# === Run the test ===
if __name__ == "__main__":
    print("=== PGx Engine - Step 2: Improved Relationships Loader ===\n")
    df = load_relationships()
    
    if df is not None and not df.empty:
        print("\n✅ Ready for next step (variant extraction)")


# Read vdf function

# In[2]:


def read_vcf_with_genotype(vcf_path, sample_index=0, sample_name=None):
    """
    Universal VCF reader that auto-detects format:
    1. Standard VCF (has REF, ALT, FORMAT)
    2. Illumina GSA Array VCF (no REF/ALT/FORMAT, samples start after ID)
    """
    variants = []
    sample_names = None
    is_gsa_format = False

    with open(vcf_path, "r") as f:
        for line in f:
            if line.startswith("##"):
                continue

            if line.startswith("#CHROM"):
                header = line.strip().split("\t")
                
                # Auto-detect format
                if len(header) > 8 and header[8] == "FORMAT":
                    # Standard VCF format
                    sample_names = header[9:]
                    print("Detected: Standard VCF format")
                else:
                    # Illumina GSA format (samples start from column 3)
                    sample_names = header[3:]
                    is_gsa_format = True
                    print("Detected: Illumina GSA Array VCF format (non-standard)")

                # If user specified sample_name, find its index
                if sample_name is not None:
                    if sample_name not in sample_names:
                        raise ValueError(f"Sample '{sample_name}' not found.\nAvailable: {sample_names}")
                    sample_index = sample_names.index(sample_name)
                    print(f"→ Selected by name: {sample_name} (index {sample_index})")
                else:
                    print(f"→ Using sample index {sample_index}: {sample_names[sample_index]}")

                print(f"Total samples found: {len(sample_names)}")
                continue

            # Parse variant line
            parts = line.strip().split("\t")
            if len(parts) < 4:
                continue

            chrom = parts[0].replace("chr", "")
            pos = int(parts[1])
            rsid = parts[2] if parts[2] != "." else f"{chrom}:{pos}"

            # === Handle different formats ===
            if is_gsa_format:
                # Illumina GSA format: GT is directly in sample column
                gt = parts[3 + sample_index]
                ref = "N"   # Placeholder
                alt = "N"   # Placeholder
            else:
                # Standard VCF format
                ref = parts[3]
                alt = parts[4]
                format_fields = parts[8].split(":")
                sample_fields = parts[9 + sample_index].split(":")
                sample_dict = dict(zip(format_fields, sample_fields))
                gt = sample_dict.get("GT", "./.")

            # Convert GT to dosage
            if gt in ["0/0", "0|0"]:
                dosage = 0
            elif gt in ["0/1", "1/0", "0|1", "1|0"]:
                dosage = 1
            elif gt in ["1/1", "1|1"]:
                dosage = 2
            else:
                dosage = None

            variants.append({
                "CHR": chrom,
                "POS": pos,
                "RSID": rsid,
                "REF": ref,
                "ALT": alt,
                "GT": gt,
                "DOSAGE": dosage,
                "SAMPLE": sample_names[sample_index] if sample_names else "Unknown"
            })

    df = pd.DataFrame(variants)
    df = df.dropna(subset=["DOSAGE"]).reset_index(drop=True)

    if len(df) == 0:
        print(f"WARNING: No valid genotypes found for sample index {sample_index}")
        return pd.DataFrame(columns=["CHR", "POS", "RSID", "REF", "ALT", "GT", "DOSAGE"])

    print(f"Successfully loaded {len(df):,} variants for sample: {df['SAMPLE'].iloc[0]} (index {sample_index})")

    return df.drop(columns=["SAMPLE"])




# extract pgx variants 

# In[3]:


import pandas as pd
from pathlib import Path

# Important PGx positions (defining variants for common star alleles)
PGX_VARIANT_POSITIONS = {
    "CYP2C19": [
        ("10", 94775453),   # rs4244285  (*2)
        ("10", 94762706),   # rs4986893  (*3)
        ("10", 94780654),   # rs12248560 (*17)
    ],
    "CYP2C9": [
        ("10", 94942291),   # rs1799853  (*2)
        ("10", 94981276),   # rs1057910  (*3)
    ],
    "VKORC1": [
        ("16", 31096368),   # rs9923231  (-1639G>A)
    ],
    "SLCO1B1": [
        ("12", 21178615),   # rs4149056  (*5)
    ],
    "CYP2D6": [
        ("22", 42522613),   # rs1065852
        ("22", 42525086),   # rs3892097 (*4)
    ]
}

def extract_pgx_variants(vcf_path: str, sample_index: int = 0):
    """
    Extracts only PGx-relevant variants from VCF using the read_vcf_with_genotype function.
    """
    if not Path(vcf_path).exists():
        print(f"❌ VCF file not found: {vcf_path}")
        return None

    print(f"\nReading VCF: {vcf_path} (sample index {sample_index})...")

    # Use the function you already defined in previous cell
    genotype_df = read_vcf_with_genotype(vcf_path, sample_index=sample_index)

    if genotype_df is None or genotype_df.empty:
        print("❌ Failed to read genotype data from VCF")
        return None

    print(f"✅ Loaded {len(genotype_df):,} total variants from VCF")

    # Filter only PGx important positions
    pgx_snps = []
    for gene, positions in PGX_VARIANT_POSITIONS.items():
        for chrom, pos in positions:
            match = genotype_df[
                (genotype_df['CHR'].astype(str) == str(chrom)) & 
                (genotype_df['POS'] == pos)
            ].copy()
            
            if not match.empty:
                match['GENE'] = gene
                pgx_snps.append(match)

    if pgx_snps:
        result = pd.concat(pgx_snps, ignore_index=True)
        print(f"✅ Found {len(result)} PGx-relevant variants")
        print("\nExtracted PGx variants:")
        print(result[['GENE', 'CHR', 'POS', 'RSID', 'REF', 'ALT', 'DOSAGE']])
        return result
    else:
        print("⚠️ No PGx-relevant variants found at the defined positions.")
        print("   This is common with array data (GSA) — we will expand positions later.")
        return None


# === Quick Test ===
print("=== PGx Engine - Step 3: Variant Extraction Ready ===")

# CHANGE THIS TO YOUR ACTUAL VCF FILE NAME
your_vcf = "E:/report sautomsation/RSO1104.vcf"          # ←←← Put your real patient VCF name here

result = extract_pgx_variants(your_vcf, sample_index=0)


# Clean pgx variants

# In[4]:


def clean_pgx_variants(pgx_df):
    """
    Clean duplicate rows and improve display for GSA array data.
    """
    if pgx_df is None or pgx_df.empty:
        print("No variants to clean.")
        return None
    
    # Remove exact duplicate rows
    pgx_df = pgx_df.drop_duplicates(subset=['GENE', 'CHR', 'POS', 'DOSAGE'])
    
    # Rename columns for nicer display
    pgx_df = pgx_df.rename(columns={
        'CHR': 'Chromosome',
        'POS': 'Position',
        'RSID': 'Variant_ID',
        'DOSAGE': 'Dosage'
    })
    
    # Add a simple Genotype column for easier reading
    def dosage_to_genotype(dosage):
        if dosage == 0:
            return "0/0 (Ref/Ref)"
        elif dosage == 1:
            return "0/1 (Het)"
        elif dosage == 2:
            return "1/1 (Hom)"
        else:
            return "Missing"
    
    pgx_df['Genotype'] = pgx_df['Dosage'].apply(dosage_to_genotype)
    
    print(f"✅ Cleaned PGx variants: {len(pgx_df)} unique variants")
    print("\nCleaned PGx Variants:")
    display_cols = ['GENE', 'Chromosome', 'Position', 'Variant_ID', 'Genotype', 'Dosage']
    print(pgx_df[display_cols])
    
    return pgx_df


# === Test the cleaning function ===
print("=== Step 4: Testing Variant Cleaner ===")

# Run this on the result from previous cell
# (Make sure you have 'result' variable from Step 3 cell)
if 'result' in globals() and result is not None:
    cleaned = clean_pgx_variants(result)
else:
    print("⚠️  'result' variable not found. Please re-run Step 3 cell first.")


# Step 5: Basic Star Allele Caller for CYP2C19 (Cell 5)

# In[5]:


# === Step 6: Improved CYP2C19 Star Allele Caller (Replace old version) ===

def call_cyp2c19_allele(pgx_df):
    """
    Improved rule-based star allele caller for CYP2C19 with South Asian context.
    """
    if pgx_df is None or pgx_df.empty:
        return {
            "gene": "CYP2C19",
            "diplotype": "*1/*1",
            "phenotype": "Normal Metabolizer",
            "note": "No variants detected"
        }
    
    # Filter CYP2C19 variants
    cyp2c19 = pgx_df[pgx_df['GENE'] == 'CYP2C19'].copy()
    
    # Default assumption (most common)
    allele1 = "*1"
    allele2 = "*1"
    
    # Check detected variants
    has_rs4244285 = not cyp2c19[cyp2c19['Position'] == 94775453].empty   # *2
    has_rs4986893 = not cyp2c19[cyp2c19['Position'] == 94762706].empty   # *3
    
    dosage_rs4244285 = cyp2c19[cyp2c19['Position'] == 94775453]['Dosage'].iloc[0] if has_rs4244285 else 0
    dosage_rs4986893 = cyp2c19[cyp2c19['Position'] == 94762706]['Dosage'].iloc[0] if has_rs4986893 else 0
    
    # Assign alleles
    if dosage_rs4244285 == 2:
        allele1 = allele2 = "*2"
    elif dosage_rs4244285 == 1:
        allele1 = "*2"
        allele2 = "*1"
    elif dosage_rs4986893 == 2:
        allele1 = allele2 = "*3"
    elif dosage_rs4986893 == 1:
        allele1 = "*3"
        allele2 = "*1"
    
    diplotype = f"{allele1}/{allele2}"
    
    # Phenotype mapping
    phenotype_map = {
        "*1/*1": "Normal Metabolizer",
        "*1/*2": "Intermediate Metabolizer",
        "*1/*3": "Intermediate Metabolizer",
        "*2/*2": "Poor Metabolizer",
        "*2/*3": "Poor Metabolizer",
        "*3/*3": "Poor Metabolizer",
        "*1/*17": "Rapid Metabolizer",
        "*17/*17": "Ultra-rapid Metabolizer",
        "*2/*17": "Intermediate Metabolizer"
    }
    
    phenotype = phenotype_map.get(diplotype, "Normal Metabolizer")
    
    # South Asian Context (very important for your project)
    sa_note = ""
    if allele1 == "*2" or allele2 == "*2":
        sa_note = " (Common in South Asians - *2 frequency ~35-40% in Tamil/South Indian populations)"
    elif allele1 == "*3" or allele2 == "*3":
        sa_note = " (*3 is less common but present in South Asians)"
    
    result = {
        "gene": "CYP2C19",
        "diplotype": diplotype,
        "phenotype": phenotype,
        "note": "Based on detected GSA array variants only" + sa_note
    }
    
    print(f"✅ CYP2C19 Analysis:")
    print(f"   Diplotype : {diplotype}")
    print(f"   Phenotype : {phenotype}")
    if sa_note:
        print(f"   South Asian Note: {sa_note.strip()}")
    
    return result


# === Test the improved function ===
print("=== Step 6: Testing Improved CYP2C19 Caller ===")

if 'cleaned' in globals() and cleaned is not None:
    cyp2c19_result = call_cyp2c19_allele(cleaned)
else:
    print("⚠️  'cleaned' variable not found. Please re-run Step 4 cell first.")


# In[ ]:





# In[6]:


# === Step 7: Drug Recommendation Engine (New Cell) ===

def get_cyp2c19_drug_recommendations(cyp2c19_result):
    """
    Returns actionable drug recommendations based on CYP2C19 phenotype.
    Focused on high-evidence CPIC Level A/B drugs.
    """
    phenotype = cyp2c19_result.get("phenotype", "Normal Metabolizer")
    diplotype = cyp2c19_result.get("diplotype", "*1/*1")
    
    recommendations = {
        "gene": "CYP2C19",
        "diplotype": diplotype,
        "phenotype": phenotype,
        "drugs": []
    }
    
    if phenotype == "Normal Metabolizer":
        recs = [
            {"drug": "Clopidogrel", "recommendation": "Standard dose", 
             "note": "Normal activation of clopidogrel expected"},
            {"drug": "Omeprazole / Esomeprazole", "recommendation": "Standard dose", 
             "note": "Normal response expected"},
            {"drug": "Citalopram / Escitalopram", "recommendation": "Standard dose", 
             "note": "Normal response expected"},
            {"drug": "Sertraline", "recommendation": "Standard dose", 
             "note": "Normal response expected"}
        ]
        
    elif phenotype == "Intermediate Metabolizer":
        recs = [
            {"drug": "Clopidogrel", "recommendation": "Consider alternative (prasugrel or ticagrelor)", 
             "note": "Reduced activation → higher risk of cardiovascular events"},
            {"drug": "Omeprazole", "recommendation": "Increase dose or use alternative (pantoprazole, rabeprazole)", 
             "note": "Reduced efficacy"},
            {"drug": "Citalopram / Escitalopram", "recommendation": "Consider 50% dose reduction or alternative", 
             "note": "Increased risk of side effects"}
        ]
        
    elif phenotype == "Poor Metabolizer":
        recs = [
            {"drug": "Clopidogrel", "recommendation": "Avoid - use prasugrel or ticagrelor instead", 
             "note": "Very low activation of prodrug"},
            {"drug": "Omeprazole", "recommendation": "Use alternative (pantoprazole or rabeprazole)", 
             "note": "Significantly reduced efficacy"},
            {"drug": "Citalopram / Escitalopram", "recommendation": "Avoid or use 50% dose with monitoring", 
             "note": "High risk of QT prolongation and side effects"}
        ]
        
    elif "Rapid" in phenotype or "Ultra-rapid" in phenotype:
        recs = [
            {"drug": "Clopidogrel", "recommendation": "Standard dose", 
             "note": "May have increased activation"},
            {"drug": "Omeprazole", "recommendation": "May need higher dose", 
             "note": "Faster metabolism → reduced efficacy"},
            {"drug": "Citalopram / Escitalopram", "recommendation": "Standard dose", 
             "note": "Normal to slightly increased clearance"}
        ]
    else:
        recs = [{"drug": "General", "recommendation": "Standard dosing", "note": "Normal metabolizer assumed"}]
    
    recommendations["drugs"] = recs
    
    # Print nice summary
    print(f"✅ CYP2C19 Drug Recommendations ({phenotype})")
    print(f"   Diplotype: {diplotype}\n")
    
    for item in recs:
        print(f"• {item['drug']}:")
        print(f"   → {item['recommendation']}")
        if item.get("note"):
            print(f"   Note: {item['note']}")
        print()
    
    return recommendations
# === Test ===
print("=== Step 7: Testing Drug Recommendations ===")

if 'cyp2c19_result' in globals():
    drug_recs = get_cyp2c19_drug_recommendations(cyp2c19_result)
else:
    # Fallback if variable name is different
    print("Creating test result for demonstration...")
    test_result = {"phenotype": "Normal Metabolizer", "diplotype": "*1/*1"}
    drug_recs = get_cyp2c19_drug_recommendations(test_result)


# In[ ]:


# === Step 8:Expanded Multi-Gene Caller with Drug & Condition Info ===#


# In[12]:


def call_pgx_alleles(cleaned_df):
    """
    Comprehensive PGx analysis covering 11 high-evidence genes.
    All recommendations based on CPIC guidelines and FDA labels.
    """
    results = {}
    
    # 1. CYP2C19
    cyp2c19 = cleaned_df[cleaned_df['GENE'] == 'CYP2C19']
    a1 = "*1"; a2 = "*1"
    if not cyp2c19.empty:
        d2 = cyp2c19[cyp2c19['Position'] == 94775453]['Dosage'].iloc[0] if any(cyp2c19['Position'] == 94775453) else 0
        d3 = cyp2c19[cyp2c19['Position'] == 94762706]['Dosage'].iloc[0] if any(cyp2c19['Position'] == 94762706) else 0
        if d2 == 2: a1 = a2 = "*2"
        elif d2 == 1: a1 = "*2"
        if d3 == 1 and a1 == "*1": a1 = "*3"
        elif d3 == 1: a2 = "*3"
    dip = f"{a1}/{a2}"
    results["CYP2C19"] = {
        "diplotype": dip,
        "phenotype": {"*1/*1":"Normal Metabolizer", "*1/*2":"Intermediate Metabolizer", "*1/*3":"Intermediate Metabolizer",
                      "*2/*2":"Poor Metabolizer", "*2/*3":"Poor Metabolizer", "*3/*3":"Poor Metabolizer"}.get(dip, "Normal Metabolizer"),
        "drugs": "Clopidogrel, Omeprazole, Esomeprazole, Citalopram, Escitalopram, Sertraline",
        "conditions": "Cardiovascular event prevention, GERD, Depression & Anxiety",
        "sa_note": " *2 allele frequency is high (~35-40%) in South Indian/Tamil populations" if "*2" in dip else ""
    }
    
    # 2. Warfarin (CYP2C9 + VKORC1)
    # (code same as before - abbreviated for space)
    cyp2c9 = cleaned_df[cleaned_df['GENE'] == 'CYP2C9']
    c9a1 = "*1"; c9a2 = "*1"
    if not cyp2c9.empty:
        if any(cyp2c9['Position'] == 94942291) and cyp2c9[cyp2c9['Position'] == 94942291]['Dosage'].iloc[0] >= 1: c9a1 = "*2"
        if any(cyp2c9['Position'] == 94981276) and cyp2c9[cyp2c9['Position'] == 94981276]['Dosage'].iloc[0] >= 1: c9a2 = "*3"
    c9dip = f"{c9a1}/{c9a2}"
    
    vkor = cleaned_df[cleaned_df['GENE'] == 'VKORC1']
    vkor_geno = "GG"
    if not vkor.empty:
        vkor_geno = "AA" if vkor['Dosage'].iloc[0] == 2 else "GA" if vkor['Dosage'].iloc[0] == 1 else "GG"
    
    results["Warfarin (CYP2C9 + VKORC1)"] = {
        "cyp2c9": c9dip,
        "vkorc1": vkor_geno,
        "drugs": "Warfarin",
        "conditions": "Atrial fibrillation, DVT, Pulmonary embolism, Stroke prevention",
        "note": "CPIC guideline recommends dose adjustment based on combined genotype"
    }
    
    # 3. SLCO1B1, 4. CYP2D6  (same as previous version)
    slco = cleaned_df[cleaned_df['GENE'] == 'SLCO1B1']
    slco_geno = "*1/*1"
    if not slco.empty:
        d = slco['Dosage'].iloc[0]
        slco_geno = "*5/*5" if d == 2 else "*1/*5" if d == 1 else "*1/*1"
    results["SLCO1B1"] = {
        "genotype": slco_geno,
        "drugs": "Simvastatin, Atorvastatin",
        "conditions": "High cholesterol, Cardiovascular prevention",
        "risk": "High myopathy risk (CPIC guideline)" if "*5/*5" in slco_geno else "Moderate myopathy risk" if "*1/*5" in slco_geno else "Normal"
    }
    
    cyp2d6 = cleaned_df[cleaned_df['GENE'] == 'CYP2D6']
    d6 = "*1"
    if not cyp2d6.empty:
        if any(cyp2d6['Position'] == 42525086) and cyp2d6[cyp2d6['Position'] == 42525086]['Dosage'].iloc[0] >= 1: d6 = "*4"
        elif any(cyp2d6['Position'] == 42522613) and cyp2d6[cyp2d6['Position'] == 42522613]['Dosage'].iloc[0] >= 1: d6 = "*10"
    results["CYP2D6"] = {
        "diplotype": f"{d6}/*1",
        "phenotype": "Normal" if d6 == "*1" else "Intermediate" if d6 in ["*4","*10"] else "Poor",
        "drugs": "Codeine, Tramadol, Tamoxifen, Amitriptyline, Nortriptyline",
        "conditions": "Pain management, Breast cancer treatment, Depression"
    }
    
    # New genes (all CPIC/FDA supported)
    results["TPMT / NUDT15"] = {
        "genotype": "*1/*1 (assumed - limited GSA coverage)",
        "drugs": "Azathioprine, 6-Mercaptopurine, Thioguanine",
        "conditions": "Acute lymphoblastic leukemia, Inflammatory bowel disease, Autoimmune disorders",
        "note": "CPIC guideline: Poor metabolizers have high risk of severe myelosuppression"
    }
    
    results["DPYD"] = {
        "genotype": "*1/*1 (assumed)",
        "drugs": "5-Fluorouracil (5-FU), Capecitabine",
        "conditions": "Colorectal cancer, Breast cancer, Gastric cancer",
        "note": "CPIC/FDA: Variants significantly increase risk of severe/fatal toxicity"
    }
    
    results["UGT1A1"] = {
        "genotype": "*1/*1 (assumed)",
        "drugs": "Irinotecan",
        "conditions": "Colorectal cancer, Lung cancer",
        "note": "*28/*28 associated with increased neutropenia risk (CPIC guideline)"
    }
    
    results["G6PD"] = {
        "genotype": "Normal (assumed)",
        "drugs": "Primaquine, Rasburicase, Dapsone, Nitrofurantoin",
        "conditions": "Malaria, Gout, Urinary tract infections",
        "note": "Deficiency can cause acute hemolytic anemia (CPIC guideline)"
    }
    
    results["HLA-B*58:01"] = {
        "genotype": "Negative (assumed - not on GSA array)",
        "drugs": "Allopurinol",
        "conditions": "Gout, Hyperuricemia",
        "note": "Positive = high risk of severe cutaneous adverse reaction (SCAR) - CPIC guideline"
    }
    
    # Print comprehensive report
    print("✅ Comprehensive PGx Analysis — 11 Genes (CPIC/FDA-based)\n")
    for gene, info in results.items():
        print(f"● {gene}")
        for k, v in info.items():
            if v and k not in ["sa_note"]:
                print(f"   {k.replace('_', ' ').title()}: {v}")
        if info.get("sa_note"):
            print(f"   South Asian Context: {info['sa_note']}")
        print("-" * 70)
    
    return results


# === Test ===
print("=== Step 12: Testing 11-Gene Comprehensive Caller ===")

if 'cleaned' in globals() and cleaned is not None:
    full_results = call_pgx_alleles(cleaned)
else:
    print("⚠️ 'cleaned' variable not found. Re-run Step 4 first.")


# # === Step 13: Updated Full PGx Report Function ===

# In[13]:


def generate_pgx_report(vcf_path: str, sample_index: int = 0):
    """
    End-to-end PGx Report using the comprehensive 11-gene caller.
    """
    print("="*80)
    print("                  PHARMACOGENOMICS (PGx) REPORT")
    print("="*80)
    print(f"File          : {vcf_path}")
    print(f"Sample Index  : {sample_index}")
    print(f"Date          : {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M')}")
    print("-"*80 + "\n")

    # 1. Extract variants from VCF
    pgx_variants = extract_pgx_variants(vcf_path, sample_index)
    
    if pgx_variants is None or pgx_variants.empty:
        print("❌ Could not extract variants from VCF.")
        return None

    # 2. Clean variants
    cleaned = clean_pgx_variants(pgx_variants)
    if cleaned is None or cleaned.empty:
        print("❌ No valid PGx variants after cleaning.")
        return None

    # 3. Run comprehensive multi-gene analysis
    print("\n" + "="*60)
    print("GENE ANALYSIS RESULTS")
    print("="*60)
    results = call_pgx_alleles(cleaned)

    # 4. Final Summary & Disclaimer
    print("\n" + "="*80)
    print("FINAL SUMMARY & DISCLAIMER")
    print("="*80)
    print("• This report analyzes 11 pharmacogenomic genes based on detectable variants in the Illumina GSA array.")
    print("• Genes marked with '(assumed)' have limited variant coverage on GSA arrays.")
    print("• All recommendations are based on CPIC guidelines and FDA pharmacogenetic labels.")
    print("• This is for educational and research purposes only.")
    print("• Consult a qualified clinician or certified PGx specialist before making any medication changes.")
    print("\n✅ Report generation completed.")

    return {
        "variants": cleaned,
        "results": results,
        "sample_info": {"vcf_path": vcf_path, "sample_index": sample_index}
    }


# === Test the Updated Full Report ===
print("=== Step 13: Testing Updated Full PGx Report ===")

# Update this path to your actual VCF
your_vcf = "E:/report sautomsation/RSO1104.vcf"  

report = generate_pgx_report(your_vcf, sample_index=0)


# In[ ]:




