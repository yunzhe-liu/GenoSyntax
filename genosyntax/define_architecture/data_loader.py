import pandas as pd
import os
import json
from collections import defaultdict
from tqdm import tqdm
from Bio import SeqIO

def load_panaroo_data(panaroo_file, strict_core_threshold=1.0, soft_core_threshold=0.98):
    """
    Parses the Panaroo gene_presence_absence_roary.csv file.
    """
    print("--- Task 1.1: Loading Panaroo Data ---")
    try:
        df_gpa = pd.read_csv(panaroo_file, low_memory=False)
    except FileNotFoundError:
        print(f"[ERROR] Panaroo file not found at: {panaroo_file}"); exit(1)

    genome_ids = df_gpa.columns[14:].tolist()
    num_genomes = len(genome_ids)
    print(f"Found {num_genomes} genomes in the Panaroo matrix.")

    locus_to_group_map = {}
    print("Building locus_tag -> group_id map...")
    for _, row in tqdm(df_gpa.iterrows(), total=len(df_gpa), desc="  Parsing rows"):
        group_id = row['Gene']
        for genome_id in genome_ids:
            if pd.notna(row[genome_id]):
                for tag in str(row[genome_id]).split(';'):
                    if tag: locus_to_group_map[tag.strip()] = group_id
    
    print("Stratifying genes into core/accessory categories...")
    strict_core, soft_core, accessory = [], [], []
    df_gpa['prevalence'] = df_gpa['No. isolates'] / num_genomes
    for _, row in df_gpa.iterrows():
        if row['prevalence'] >= strict_core_threshold: strict_core.append(row['Gene'])
        elif row['prevalence'] >= soft_core_threshold: soft_core.append(row['Gene'])
        else: accessory.append(row['Gene'])
            
    stratified_genes = {'strict_core': strict_core, 'soft_core': soft_core, 'accessory': accessory}
    print("Panaroo data loading complete.")
    return genome_ids, stratified_genes, locus_to_group_map


def load_gff_data(gff_dir, locus_to_group_map, genome_ids):
    """
    Parses all GFF files from a directory to build gene order and coordinate databases.
    """
    print("\n--- Task 1.2: Loading GFF Data ---")
    gene_order_db, gene_coord_db = {}, defaultdict(list)
    available_subdirs = {d for d in os.listdir(gff_dir) if os.path.isdir(os.path.join(gff_dir, d))}
    
    for genome_id in tqdm(genome_ids, desc="Parsing GFFs"):
        if genome_id not in available_subdirs:
            print(f"  [Warning] No subdirectory named '{genome_id}' found in {gff_dir}. Skipping."); gene_order_db[genome_id], gene_coord_db[genome_id] = [], []; continue

        gff_path = os.path.join(gff_dir, genome_id, f"{genome_id}.gff")
        if not os.path.exists(gff_path):
            print(f"  [Warning] GFF file not found at {gff_path}. Skipping."); gene_order_db[genome_id], gene_coord_db[genome_id] = [], []; continue
            
        genes_with_pos = []
        with open(gff_path, 'r') as f:
            for line in f:
                if line.startswith('#') or '\tCDS\t' not in line: continue
                parts = line.strip().split('\t')
                attributes = {k: v for k, v in [attr.split('=') for attr in parts[8].split(';') if '=' in attr]}
                locus_tag = attributes.get('ID')
                if locus_tag and locus_tag in locus_to_group_map:
                    group_id = locus_to_group_map[locus_tag]
                    gene_info = { 'group_id': group_id, 'contig': parts[0], 'start': int(parts[3]), 'end': int(parts[4]), 'strand': parts[6], 'locus_tag': locus_tag }
                    genes_with_pos.append(gene_info)
                    
        genes_with_pos.sort(key=lambda x: x['start'])
        gene_coord_db[genome_id].extend(genes_with_pos)
        gene_order_db[genome_id] = [g['group_id'] for g in genes_with_pos]
        
    print("GFF data loading complete.")
    return gene_order_db, gene_coord_db


def load_protein_sequences(protein_faa_dir, locus_to_group_map):
    """
    Loads all protein sequences from phanotate.faa files for ASI calculation.
    """
    print("\n--- Task: Loading Protein Sequences for ASI ---")
    protein_db = defaultdict(dict)
    genome_dirs = [d for d in os.listdir(protein_faa_dir) if os.path.isdir(os.path.join(protein_faa_dir, d))]
    
    for genome_id in tqdm(genome_dirs, desc="Loading FAAs"):
        faa_path = os.path.join(protein_faa_dir, genome_id, "phanotate.faa")
        if not os.path.exists(faa_path): continue
            
        for record in SeqIO.parse(faa_path, "fasta"):
            locus_tag = record.id
            if locus_tag in locus_to_group_map:
                group_id = locus_to_group_map[locus_tag]
                protein_db[group_id][genome_id] = str(record.seq)
                
    print(f"Loaded protein sequences for {len(protein_db)} gene families.")
    return protein_db


def load_genome_sequences(fasta_dir, genome_ids):
    """
    Loads all genome sequences from a directory for boundary signal analysis.
    """
    print("\n--- Task: Loading Genome Sequences for Boundary Analysis ---")
    genome_seqs = {}
    
    for genome_id in tqdm(genome_ids, desc="Loading FASTAs"):
        fasta_path = None
        for ext in ['.fasta', '.fa', '.fna']:
            potential_path = os.path.join(fasta_dir, f"{genome_id}{ext}")
            if os.path.exists(potential_path):
                fasta_path = potential_path
                break
        
        if fasta_path is None: continue
            
        try:
            genome_seqs[genome_id] = SeqIO.read(fasta_path, "fasta")
        except Exception as e:
            print(f"  [Error] Could not parse FASTA for {genome_id}: {e}")

    print(f"Loaded {len(genome_seqs)} genome sequences.")
    return genome_seqs
