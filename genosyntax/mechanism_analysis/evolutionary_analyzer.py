import os
import json
import subprocess
from io import StringIO
import pandas as pd
from Bio import SeqIO, AlignIO
from tqdm import tqdm
from itertools import combinations
from collections import Counter
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

# --- Helper Function for Protein Distance ---
def _calculate_pairwise_identity(alignment):
    """Calculates the average pairwise identity from a Bio.AlignIO alignment object."""
    num_seqs = len(alignment)
    if num_seqs < 2: return None
    
    total_identity, num_pairs = 0, 0
    for i, j in combinations(range(num_seqs), 2):
        seq1, seq2 = str(alignment[i].seq), str(alignment[j].seq)
        matches, length = 0, 0
        for char1, char2 in zip(seq1, seq2):
            if char1 != '-' or char2 != '-':
                length += 1
                if char1 == char2:
                    matches += 1
        if length > 0:
            total_identity += (matches / length)
            num_pairs += 1
            
    return total_identity / num_pairs if num_pairs > 0 else None

# --- Core Analysis Functions ---
def _calculate_protein_evolutionary_distance(accessory_modules, gene_coord_db, protein_seqs, temp_dir):
    """
    Calculates Protein Evolutionary Distance (PED) for each accessory module.
    PED = 1 - Average Pairwise Identity.
    """
    print("  -> Sub-task: Calculating Protein Evolutionary Distance (PED)...")
    results = []
    
    for module in tqdm(accessory_modules, desc="  Calculating PED"):
        module_id, member_genes = module['module_id'], module['members']
        
        all_family_records = []
        for gene_family in member_genes:
            family_records_for_gene = []
            for genome_id, gene_list in gene_coord_db.items():
                gene_info = next((g for g in gene_list if g.get('group_id') == gene_family), None)
                if gene_info and gene_info.get('locus_tag'):
                    locus_tag = gene_info['locus_tag']
                    if genome_id in protein_seqs and locus_tag in protein_seqs[genome_id]:
                        family_records_for_gene.append(protein_seqs[genome_id][locus_tag])
            all_family_records.append(family_records_for_gene)

        module_distances = []
        for i, records in enumerate(all_family_records):
            if len(records) < 2: continue
            
            gene_family = member_genes[i]
            unaligned_path = os.path.join(temp_dir, f"{gene_family}.faa")
            SeqIO.write(records, unaligned_path, "fasta")
            
            try:
                mafft_process = subprocess.run(
                    ['mafft', '--auto', '--quiet', unaligned_path], 
                    capture_output=True, text=True, check=True, encoding='utf-8'
                )
                if not mafft_process.stdout: continue
                alignment = AlignIO.read(StringIO(mafft_process.stdout), "fasta")
                avg_identity = _calculate_pairwise_identity(alignment)
                if avg_identity is not None:
                    module_distances.append(1 - avg_identity)
            except (subprocess.CalledProcessError, ValueError, FileNotFoundError) as e:
                if isinstance(e, FileNotFoundError):
                    print("\n[FATAL ERROR] MAFFT command not found. Is it installed and in your PATH?")
                    exit(1)
                continue
        
        if module_distances:
            avg_module_ped = sum(module_distances) / len(module_distances)
            results.append({'module_id': module_id, 'protein_evolutionary_distance': avg_module_ped})

    return pd.DataFrame(results, columns=['module_id', 'protein_evolutionary_distance'])

def _analyze_high_spi_slot_enrichment(df_spi, spi_threshold):
    """
    Performs Fisher's Exact Test to find if high-SPI modules are enriched in certain slot types.
    """
    print(f"  -> Sub-task: Analyzing slot enrichment for high-SPI (>{spi_threshold}) modules...")
    df = df_spi.copy()
    df = df[df['observed_slots'] != 'None'].copy()
    df['spi_group'] = df['spi'].apply(lambda x: 'high' if x > spi_threshold else 'low')
    
    high_spi_slots, low_spi_slots = [], []
    for _, row in df.iterrows():
        slots = [s.strip() for s in row['observed_slots'].split(';')]
        if row['spi_group'] == 'high':
            high_spi_slots.extend(slots)
        else:
            low_spi_slots.extend(slots)

    high_counts, low_counts = Counter(high_spi_slots), Counter(low_spi_slots)
    all_slots = sorted(list(set(high_counts.keys()) | set(low_counts.keys())))
    
    enrichment_results = []
    total_high_obs, total_low_obs = sum(high_counts.values()), sum(low_counts.values())
    
    if total_high_obs == 0:
        print("    [INFO] No high-SPI modules found. Skipping enrichment analysis.")
        return pd.DataFrame()

    for slot in all_slots:
        a = high_counts.get(slot, 0)
        b = total_high_obs - a
        c = low_counts.get(slot, 0)
        d = total_low_obs - c
        
        odds_ratio, p_value = fisher_exact([[a, b], [c, d]], alternative='greater')
        enrichment_results.append({
            'slot_type': slot, 'high_spi_count_in_slot': a, 'low_spi_count_in_slot': c,
            'odds_ratio': odds_ratio, 'p_value': p_value
        })
        
    df_enrichment = pd.DataFrame(enrichment_results)
    if not df_enrichment.empty:
        _, q_values, _, _ = multipletests(df_enrichment['p_value'], alpha=0.05, method='fdr_bh')
        df_enrichment['q_value'] = q_values
        df_enrichment = df_enrichment.sort_values('p_value').reset_index(drop=True)

    return df_enrichment

# --- Main Entrypoint Function ---
def run_evolutionary_analysis(output_dir, df_spi, accessory_modules, gene_coord_db, pharokka_dir, spi_threshold):
    """
    Executes all evolutionary dynamics analyses.
    """
    # Setup directories
    reports_dir = os.path.join(output_dir, "reports")
    temp_dir = os.path.join(output_dir, "intermediate", "mafft_temp")
    os.makedirs(reports_dir, exist_ok=True)
    os.makedirs(temp_dir, exist_ok=True)

    # Pre-load all protein sequences for efficiency
    print("  -> Pre-loading protein sequences from Pharokka .faa files...")
    protein_seqs = {}
    genome_dirs = [d for d in os.listdir(pharokka_dir) if os.path.isdir(os.path.join(pharokka_dir, d))]
    for genome_id in tqdm(genome_dirs, desc="  Loading proteins"):
        faa_path = os.path.join(pharokka_dir, genome_id, "phanotate.faa")
        if os.path.exists(faa_path):
            try:
                protein_seqs[genome_id] = SeqIO.to_dict(SeqIO.parse(faa_path, "fasta"))
            except Exception:
                continue

    # Run analyses
    df_prot_dist = _calculate_protein_evolutionary_distance(accessory_modules, gene_coord_db, protein_seqs, temp_dir)
    df_slot_enrichment = _analyze_high_spi_slot_enrichment(df_spi, spi_threshold)

    # Merge and save final reports
    df_master = pd.merge(df_spi, df_prot_dist, on='module_id', how='left')
    
    master_path = os.path.join(reports_dir, "spi_vs_ped_analysis.csv")
    enrichment_path = os.path.join(reports_dir, "high_spi_slot_enrichment.csv")
    
    df_master.to_csv(master_path, index=False)
    if not df_slot_enrichment.empty:
        df_slot_enrichment.to_csv(enrichment_path, index=False)

    print(f"  -> Evolutionary dynamics report saved to: {master_path}")
    if not df_slot_enrichment.empty:
        print(f"  -> High-SPI slot enrichment report saved to: {enrichment_path}")
