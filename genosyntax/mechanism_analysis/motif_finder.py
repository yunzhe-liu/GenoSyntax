import os
import json
from tqdm import tqdm
from Bio import SeqIO
import pandas as pd
from collections import defaultdict

def _sanitize_for_header(text: str) -> str:
    """
    Sanitizes a string to be a single, safe word for a FASTA header.
    """
    text = text.replace('(', '').replace(')', '')
    text = text.replace(' -> ', '-to-')
    text = text.replace(' ', '-')
    text = text.replace(':', '-')
    return text


def _sanitize_fasta_for_meme(file_path: str, min_length: int):
    """
    Reads a FASTA file, filters out sequences shorter than min_length,
    and overwrites the file with linearized, clean FASTA records.
    Write records in a format accepted by MEME-ChIP.
    """
    # Check if the file was created and is not empty
    if not os.path.exists(file_path) or os.path.getsize(file_path) == 0:
        return # Do nothing if the file doesn't exist or is empty

    try:
        # Read all records, keeping only those that meet the length requirement
        kept_records = []
        for record in SeqIO.parse(file_path, "fasta"):
            if len(record.seq) >= min_length:
                kept_records.append(record)

        # Overwrite the original file with the filtered and sanitized records
        # SeqIO.write handles the linearization automatically.
        with open(file_path, 'w') as f:
            SeqIO.write(kept_records, f, "fasta")

        if kept_records:
             print(f"  - Sanitized {os.path.basename(file_path)}: Kept {len(kept_records)} sequences >= {min_length}bp.")

    except Exception as e:
        print(f"  [ERROR] Could not sanitize file {file_path}: {e}")


def extract_boundary_sequences(
    slot_database,
    gene_order_db,
    gene_coord_db,
    pharokka_output_dir,
    output_dir,
    window_size=200,
    p_value_threshold=0.05,
    min_meme_seq_length=20, 
    statistical_report_df=None,
    target_slot_instances=None
):
    """Extracts boundary sequences from genomic slots."""

    
    # --- Determine Operational Mode ---
    if target_slot_instances is not None:
        mode = "Targeted"
        database_to_process = target_slot_instances
        print(f"\n--- Task: Extracting Boundary Sequences (Targeted Mode, Window: {window_size}bp) ---")
    elif statistical_report_df is not None and not statistical_report_df.empty:
        mode = "Hotspot/Coldspot"
        database_to_process = slot_database
        print(f"\n--- Task: Extracting Boundary Sequences (Hotspot/Coldspot Mode, Window: {window_size}bp) ---")
        print("  Step 1: Classifying slot types...")
        hot_slots = set(statistical_report_df[statistical_report_df['p_value_enrichment'] < p_value_threshold]['slot_type'])
        cold_slots = set(statistical_report_df[statistical_report_df['p_value_depletion'] < p_value_threshold]['slot_type'])
    else:
        print("\n[ERROR] No valid input provided for motif extraction.")
        return

    print("  Step 2: Loading genome sequences...")
    all_possible_genomes = list(slot_database.keys())
    nucleotide_seqs = {}
    for genome_id in tqdm(all_possible_genomes, desc="  Loading NT sequences"):
        fna_path = next((os.path.join(pharokka_output_dir, genome_id, f) for f in os.listdir(os.path.join(pharokka_output_dir, genome_id)) if f.endswith(('.fna', '.fasta', 'fa')) and 'reoriented' in f), None)
        if not fna_path:
            fna_path = next((os.path.join(pharokka_output_dir, genome_id, f) for f in os.listdir(os.path.join(pharokka_output_dir, genome_id)) if f.endswith(('.fna', '.fasta', 'fa'))), None)
        if fna_path and os.path.exists(fna_path):
            try:
                records = list(SeqIO.parse(fna_path, "fasta"))
                if len(records) == 1: nucleotide_seqs[genome_id] = records[0]
            except Exception: # Ignore parsing errors for corrupted files
                pass

    print("\n  Step 3: Extracting boundary sequences...")
    output_seqs = defaultdict(list)
    
    for genome_id, slots in tqdm(database_to_process.items(), desc="  Processing genomes"):
        if genome_id not in nucleotide_seqs: continue
        full_gene_order = gene_order_db.get(genome_id, [])
        coords_by_group = {g['group_id']: g for g in reversed(gene_coord_db.get(genome_id, []))}
        genome_sequence = str(nucleotide_seqs[genome_id].seq)
        
        for slot_id, slot_info in slots.items():
            if mode == "Hotspot/Coldspot":
                if slot_info['slot_id_low_res'] in hot_slots: slot_class = "hot"
                elif slot_info['slot_id_low_res'] in cold_slots: slot_class = "cold"
                else: slot_class = "neutral"
            else:
                slot_class = "high_spi_freq"

            if not slot_info['content']: continue
            try:
                order_map = {gene: i for i, gene in enumerate(full_gene_order)}
                first_gene_global_idx = order_map[slot_info['content'][0]]
                last_gene_global_idx = first_gene_global_idx + len(slot_info['content']) - 1
            except (KeyError, ValueError):
                continue

            backend_gene_id = full_gene_order[first_gene_global_idx - 1] if first_gene_global_idx > 0 else None
            frontend_gene_id = full_gene_order[last_gene_global_idx + 1] if last_gene_global_idx < len(full_gene_order) - 1 else None
            sanitized_slot_id = _sanitize_for_header(slot_info['slot_id_low_res'])
            
            if backend_gene_id and coords_by_group.get(backend_gene_id):
                backend_coord = coords_by_group[backend_gene_id]
                first_gene_coord = coords_by_group.get(slot_info['content'][0])
                if first_gene_coord:
                    extract_start = backend_coord['end']
                    limit = first_gene_coord['start'] - 1
                    extract_end = min(extract_start + window_size, limit)
                    if extract_end > extract_start:
                        seq = genome_sequence[extract_start:extract_end]
                        header = f">{genome_id}__{sanitized_slot_id}__backend__abs_start={extract_start}"
                        output_seqs[f"{slot_class}_backend"].append(f"{header}\n{seq}")

            if frontend_gene_id and coords_by_group.get(frontend_gene_id):
                frontend_coord = coords_by_group[frontend_gene_id]
                last_gene_coord = coords_by_group.get(slot_info['content'][-1])
                if last_gene_coord:
                    extract_end = frontend_coord['start'] - 1
                    limit = last_gene_coord['end']
                    extract_start = max(extract_end - window_size, limit)
                    if extract_end > extract_start:
                        seq = genome_sequence[extract_start:extract_end]
                        header = f">{genome_id}__{sanitized_slot_id}__frontend__abs_start={extract_start}"
                        output_seqs[f"{slot_class}_frontend"].append(f"{header}\n{seq}")

    print("\n  Step 4: Writing initial output FASTA files...")
    created_files = []
    for category, sequences in output_seqs.items():
        if sequences:
            file_path = os.path.join(output_dir, f"{category}_boundaries.fasta")
            with open(file_path, 'w') as f:
                f.write("\n".join(sequences))
            print(f"  - Wrote {len(sequences)} sequences to {os.path.basename(file_path)}")
            created_files.append(file_path)

    # Final sanitization pass to ensure MEME-ChIP compatibility
    print(f"\n  Step 5: Sanitizing FASTA files for MEME-ChIP (min length: {min_meme_seq_length}bp)...")
    if not created_files:
        print("  - No files were created, skipping sanitization.")
    else:
        for f_path in created_files:
            _sanitize_fasta_for_meme(f_path, min_length=min_meme_seq_length)

    print("\n--- Motif Finder data preparation finished ---")
