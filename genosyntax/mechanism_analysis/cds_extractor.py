# genosyntax/mechanism_analysis/cds_extractor.py

import os
import json
import pandas as pd
from tqdm import tqdm

# =============================================================================
# --- 1. FUNCTIONAL CLASSIFICATION LOGIC ---
# =============================================================================

FUNCTIONAL_MAP = {
    "Structural & Infection": [
        'tail', 'baseplate', 'capsid', 'head', 'portal', 'terminase', 
        'lysis', 'holin', 'endolysin', 'sheath', 'fiber', 'tape measure',
        'structural', 'virion'
    ],
    "Replication & Regulation": [
        'dna polymerase', 'helicase', 'primase', 'ligase', 'replication',
        'rna polymerase', 'transcription', 'sigma factor', 'regulator',
        'repressor', 'dna binding'
    ],
    "Genome Maintenance & Metabolism": [
        'recombination', 'exonuclease', 'endonuclease', 'repair', 
        'nucleotide', 'synthesis', 'kinase', 'synthase', 'reductase',
        'metabolism', 'folate', 'integrase', 'recombinase', 'transposase'
    ],
    "Hypothetical & Unknown": [
        'hypothetical', 'unknown', 'uncharacterized', 'predicted protein'
    ]
}

def _get_function_category(product_description: str) -> str:
    """
    Assigns a high-level functional category to a gene based on its product description.
    """
    if not isinstance(product_description, str):
        return "Hypothetical & Unknown"
        
    desc_lower = product_description.lower()
    
    for category, keywords in FUNCTIONAL_MAP.items():
        if any(keyword in desc_lower for keyword in keywords):
            return category
            
    # If no keywords match, default to a general category
    return "Other"

# =============================================================================
# --- 2. MAIN SCRIPT LOGIC  ---
# =============================================================================

def run_cds_extraction(output_dir, arch_def, gene_coord_db, gene_annotation_map):
    """
    Main function to load GenoSyntax outputs and generate a comprehensive CDS table.
    """
    reports_dir = os.path.join(output_dir, "reports")
    os.makedirs(reports_dir, exist_ok=True)
    
    # --- Step 1: Extract the list of high-quality genomes ---
    passed_genomes = arch_def.get('metadata', {}).get('qc_summary', {}).get('passed_genomes_list', [])
    if not passed_genomes:
        print("  [INFO] No high-quality genomes found in the architecture definition file. Skipping CDS layout extraction.")
        return
    
    print(f"  -> Found {len(passed_genomes)} high-quality genomes to process.")

    # --- Step 2: Iterate through genomes and genes to build the table ---
    all_cds_data = []
    
    for genome_id in tqdm(passed_genomes, desc="  Extracting CDS", leave=False):
        gene_list = gene_coord_db.get(genome_id)
        if not gene_list:
            continue

        for gene_info in gene_list:
            group_id = gene_info.get('group_id')
            if not group_id:
                continue
                
            product = gene_annotation_map.get(group_id, "Unknown product")
            category = _get_function_category(product)
            
            all_cds_data.append({
                'genome_id': genome_id,
                'start': gene_info.get('start'),
                'end': gene_info.get('end'),
                'strand': gene_info.get('strand'),
                'product': product,
                'function_category': category,
                'gene_family_id': group_id # Add gene family ID for extra utility
            })
    
    # --- Step 3: Convert to DataFrame and save to CSV ---
    if not all_cds_data:
        print("  [INFO] No CDS data was extracted. The output file will not be created.")
        return

    df = pd.DataFrame(all_cds_data)
    
    # Ensure a useful column order
    df = df[['genome_id', 'start', 'end', 'strand', 'gene_family_id', 'product', 'function_category']]
    
    output_path = os.path.join(reports_dir, "all_genomes_cds_layout.csv")
    df.to_csv(output_path, index=False)
    
    print(f"  -> Comprehensive CDS table for {df['genome_id'].nunique()} genomes saved to: {output_path}")
