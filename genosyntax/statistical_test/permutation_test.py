import numpy as np
import pandas as pd
from collections import Counter
from tqdm import tqdm
import random

def run_permutation_test(slot_database, accessory_modules, num_permutations=10000):
    """Performs a permutation test to identify significant hotspots and coldspots."""
    print("\n--- Task: Running Permutation Test for Slot Significance ---")

    # ==============================================================================
    # --- ANALYSIS A: ALL ACCESSORY MODULES (>= 1 GENE) ---
    # This is the original, unchanged analysis.
    # ==============================================================================
    print("\n--- Analysis A: Testing ALL accessory modules (single- and multi-gene) ---")
    
    # --- 1A. Tally Observations and Prepare Data Pools ---
    print("  Step 1A: Tallying observations for all accessory genes...")
    
    all_module_member_sets = {g for m in accessory_modules for g in m['members']}
    
    observed_counts_per_slot_type = Counter()
    all_genes_in_slots = []
    slot_type_of_each_gene = []

    for genome_id, slots in slot_database.items():
        for slot_id, slot_info in slots.items():
            slot_type = slot_info['slot_id_low_res']
            for gene in slot_info['content']:
                all_genes_in_slots.append(gene)
                slot_type_of_each_gene.append(slot_type)
                if gene in all_module_member_sets:
                    observed_counts_per_slot_type[slot_type] += 1

    total_genes_in_slots = len(all_genes_in_slots)
    total_module_genes_in_slots = sum(1 for g in all_genes_in_slots if g in all_module_member_sets)
    unique_slot_types = sorted(list(set(slot_type_of_each_gene)))
    
    print(f"  Found {total_module_genes_in_slots} total accessory gene instances among {total_genes_in_slots} available positions.")
    
    report_df_all = pd.DataFrame() # Initialize in case of early exit
    if total_module_genes_in_slots > 0 and total_genes_in_slots > 0:
        # --- 2A. Run Permutations ---
        print(f"\n  Step 2A: Running {num_permutations} permutations for all accessory modules...")
        
        simulation_matrix = np.zeros((num_permutations, len(unique_slot_types)), dtype=int)
        slot_type_map = {st: i for i, st in enumerate(unique_slot_types)}
        all_indices = list(range(total_genes_in_slots))
        
        for i in tqdm(range(num_permutations), desc="  Simulating (All Modules)"):
            random_indices = random.sample(all_indices, total_module_genes_in_slots)
            simulated_slot_types = [slot_type_of_each_gene[j] for j in random_indices]
            sim_counts = Counter(simulated_slot_types)
            for slot_type, count in sim_counts.items():
                if slot_type in slot_type_map:
                    simulation_matrix[i, slot_type_map[slot_type]] = count
        
        # --- 3A. Calculate p-values ---
        print("\n  Step 3A: Calculating significance for all accessory modules...")
        results_all = []
        total_genes_per_slot_type = Counter(slot_type_of_each_gene)
        for i, slot_type in enumerate(unique_slot_types):
            observed_count = observed_counts_per_slot_type.get(slot_type, 0)
            simulation_counts = simulation_matrix[:, i]
            p_enrichment = (np.sum(simulation_counts >= observed_count) + 1) / (num_permutations + 1)
            p_depletion = (np.sum(simulation_counts <= observed_count) + 1) / (num_permutations + 1)
            mean_simulated = np.mean(simulation_counts)
            fold_enrichment = observed_count / mean_simulated if mean_simulated > 0 else float('inf')
            results_all.append({
                'slot_type': slot_type, 'observed_module_gene_count': observed_count,
                'mean_simulated_count': round(mean_simulated, 2), 'total_genes_in_slot': total_genes_per_slot_type[slot_type],
                'p_value_enrichment': round(p_enrichment, 5), 'p_value_depletion': round(p_depletion, 5),
                'fold_enrichment': round(fold_enrichment, 2) if fold_enrichment != float('inf') else 'inf'
            })
        report_df_all = pd.DataFrame(results_all).sort_values(by='p_value_enrichment', ascending=True)


    print("\n--- Analysis B: Testing MULTI-GENE modules ONLY ---")

    # --- 1B. Filter for multi-gene modules and re-tally observations ---
    print("  Step 1B: Tallying observations for multi-gene module members...")
    
    multi_gene_modules = [m for m in accessory_modules if m.get('num_genes', 0) > 1]
    multi_gene_member_sets = {g for m in multi_gene_modules for g in m['members']}

    observed_counts_per_slot_type_multi = Counter()
    for gene, slot_type in zip(all_genes_in_slots, slot_type_of_each_gene):
        if gene in multi_gene_member_sets:
            observed_counts_per_slot_type_multi[slot_type] += 1
    
    total_multi_gene_module_events = sum(observed_counts_per_slot_type_multi.values())
    
    print(f"  Found {total_multi_gene_module_events} multi-gene module instances among {total_genes_in_slots} available positions.")

    report_df_multi = pd.DataFrame() # Initialize in case of early exit
    if total_multi_gene_module_events > 0 and total_genes_in_slots > 0:
        # --- 2B. Run Permutations for multi-gene modules ---
        print(f"\n  Step 2B: Running {num_permutations} permutations for multi-gene modules...")
        
        # Reuse the same simulation setup, just change the number of events to sample
        simulation_matrix_multi = np.zeros((num_permutations, len(unique_slot_types)), dtype=int)
        
        for i in tqdm(range(num_permutations), desc="  Simulating (Multi-Gene)"):
            random_indices = random.sample(all_indices, total_multi_gene_module_events)
            simulated_slot_types = [slot_type_of_each_gene[j] for j in random_indices]
            sim_counts = Counter(simulated_slot_types)
            for slot_type, count in sim_counts.items():
                if slot_type in slot_type_map:
                    simulation_matrix_multi[i, slot_type_map[slot_type]] = count

        # --- 3B. Calculate p-values for multi-gene modules ---
        print("\n  Step 3B: Calculating significance for multi-gene modules...")
        results_multi = []
        for i, slot_type in enumerate(unique_slot_types):
            observed_count = observed_counts_per_slot_type_multi.get(slot_type, 0)
            simulation_counts = simulation_matrix_multi[:, i]
            p_enrichment = (np.sum(simulation_counts >= observed_count) + 1) / (num_permutations + 1)
            p_depletion = (np.sum(simulation_counts <= observed_count) + 1) / (num_permutations + 1)
            mean_simulated = np.mean(simulation_counts)
            fold_enrichment = observed_count / mean_simulated if mean_simulated > 0 else float('inf')
            results_multi.append({
                'slot_type': slot_type, 'observed_module_gene_count': observed_count,
                'mean_simulated_count': round(mean_simulated, 2), 'total_genes_in_slot': total_genes_per_slot_type[slot_type],
                'p_value_enrichment': round(p_enrichment, 5), 'p_value_depletion': round(p_depletion, 5),
                'fold_enrichment': round(fold_enrichment, 2) if fold_enrichment != float('inf') else 'inf'
            })
        report_df_multi = pd.DataFrame(results_multi).sort_values(by='p_value_enrichment', ascending=True)

    print("\nPermutation tests finished.")
    return report_df_all, report_df_multi
