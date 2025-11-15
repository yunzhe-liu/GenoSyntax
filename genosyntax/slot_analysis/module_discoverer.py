import pandas as pd
from itertools import combinations
from scipy.stats import fisher_exact
import networkx as nx
import community as community_louvain
from tqdm import tqdm
import numpy as np
from collections import defaultdict

def _calculate_phi_coefficient(contingency_table):
    """Calculates Phi coefficient from a 2x2 contingency table."""
    n11, n10, n01, n00 = contingency_table.ravel()
    denominator = (n11 + n10) * (n01 + n00) * (n11 + n01) * (n10 + n00)
    return (n11 * n00 - n10 * n01) / np.sqrt(denominator) if denominator > 0 else 0

def discover_multi_gene_modules(
    pangenome_strat, 
    gene_order_db, 
    gene_coord_db, 
    passed_genomes, 
    phi_threshold=0.8, 
    p_value_threshold=0.01, 
    adjacency_distance=250, 
    support_threshold=0.8
):
    """
    Discovers multi-gene accessory modules (AMs) using a dual-validation workflow.
    1. Statistical co-occurrence network clustering.
    2. Physical adjacency validation.
    """
    print("\n--- Task 3: Discovering Multi-Gene Accessory Modules ---")
    
    accessory_genes = set(pangenome_strat['accessory'])
    num_genomes = len(passed_genomes)
    
    # --- 1a. Build gene presence matrix for accessory genes ---
    print("  Step 1: Calculating co-occurrence...")
    presence_matrix = defaultdict(list)
    for gene in accessory_genes:
        for genome in passed_genomes:
            presence_matrix[gene].append(1 if gene in gene_order_db.get(genome, []) else 0)
    
    df_presence = pd.DataFrame(presence_matrix, index=passed_genomes).T

    # --- 1b. Calculate significant co-occurring pairs ---
    significant_pairs = []
    for gene1, gene2 in tqdm(combinations(accessory_genes, 2), total=len(accessory_genes)*(len(accessory_genes)-1)//2, desc="  Finding pairs"):
        v1, v2 = df_presence.loc[gene1], df_presence.loc[gene2]
        n11, n10, n01, n00 = ((v1==1)&(v2==1)).sum(), ((v1==1)&(v2==0)).sum(), ((v1==0)&(v2==1)).sum(), ((v1==0)&(v2==0)).sum()
        
        contingency_table = [[n11, n10], [n01, n00]]
        phi = _calculate_phi_coefficient(np.array(contingency_table))
        
        if phi >= phi_threshold:
            _, p_val = fisher_exact(contingency_table, alternative='greater')
            if p_val < p_value_threshold:
                significant_pairs.append({'gene1': gene1, 'gene2': gene2, 'phi': phi, 'p_value': p_val})

    if not significant_pairs:
        print("  No significant co-occurring pairs found. No multi-gene modules discovered.")
        return []

    # --- 2. Community Detection ---
    print(f"\n  Step 2: Found {len(significant_pairs)} significant pairs. Detecting communities...")
    edges_df = pd.DataFrame(significant_pairs)
    G = nx.from_pandas_edgelist(edges_df, 'gene1', 'gene2', edge_attr='phi')
    partition = community_louvain.best_partition(G, weight='phi')
    
    candidate_modules = defaultdict(list)
    for gene, mod_id in partition.items():
        candidate_modules[mod_id].append(gene)
    
    # Filter out single-gene communities, as they will be handled later
    multi_gene_candidates = [members for members in candidate_modules.values() if len(members) > 1]
    print(f"  Found {len(multi_gene_candidates)} candidate multi-gene modules.")

    # --- 3. Physical Adjacency Validation ---
    print("\n  Step 3: Performing physical adjacency validation...")
    validated_modules = []
    for module_genes in tqdm(multi_gene_candidates, desc="  Validating modules"):
        module_genes_set = set(module_genes)
        genomes_to_check = [
            gid for gid in passed_genomes 
            if module_genes_set.issubset(set(gene_order_db.get(gid, [])))
        ]
        
        if not genomes_to_check:
            continue
            
        support_count = 0
        for genome_id in genomes_to_check:
            coords = gene_coord_db[genome_id]
            module_coords = [g for g in coords if g['group_id'] in module_genes_set]
            
            if len(module_coords) != len(module_genes): continue # Handle paralogs simply by size check

            module_coords.sort(key=lambda x: x['start'])
            
            # Check if all are on the same contig
            if len(set(g['contig'] for g in module_coords)) != 1:
                continue

            # Check adjacency
            is_adjacent = True
            for i in range(len(module_coords) - 1):
                gene_curr, gene_next = module_coords[i], module_coords[i+1]
                distance = gene_next['start'] - gene_curr['end']
                if distance < 0 or distance > adjacency_distance: # Genes overlap or are too far
                    is_adjacent = False
                    break
                
                # Check for intervening genes
                intervening_genes = [
                    g for g in coords 
                    if g['contig'] == gene_curr['contig'] and 
                       g['start'] > gene_curr['end'] and 
                       g['end'] < gene_next['start']
                ]
                if intervening_genes:
                    is_adjacent = False
                    break
            
            if is_adjacent:
                support_count += 1
        
        support_fraction = support_count / len(genomes_to_check) if genomes_to_check else 0
        if support_fraction >= support_threshold:
            validated_modules.append(sorted(module_genes))

    print(f"  Validation complete. Discovered {len(validated_modules)} high-confidence multi-gene modules.")
    return validated_modules


def define_accessory_modules(
    pangenome_strat, 
    multi_gene_modules,
    gene_annotation_map
):
    """
    Finalizes the Accessory Module (AM) list by combining multi-gene modules
    and defining all other accessory genes as single-gene modules.
    """
    print("\n--- Task 4: Finalizing Accessory Module Definitions ---")
    
    accessory_genes = set(pangenome_strat['accessory'])
    genes_in_multi_gene_modules = {gene for mod in multi_gene_modules for gene in mod}
    
    accessory_modules = []
    module_counter = 1
    
    # --- 1. Add validated multi-gene modules ---
    for members in multi_gene_modules:
        # Find a consensus function
        functions = [gene_annotation_map.get(g, 'hypothetical protein') for g in members]
        consensus_func = next((f for f in functions if 'hypothetical' not in f.lower()), 'hypothetical protein')
        
        accessory_modules.append({
            'module_id': f"AM_{module_counter:04d}",
            'module_type': 'multi-gene',
            'num_genes': len(members),
            'members': members,
            'consensus_function': consensus_func
        })
        module_counter += 1

    # --- 2. Add remaining accessory genes as single-gene modules ---
    single_genes = accessory_genes - genes_in_multi_gene_modules
    for gene in sorted(list(single_genes)):
        accessory_modules.append({
            'module_id': f"AM_{module_counter:04d}",
            'module_type': 'single-gene',
            'num_genes': 1,
            'members': [gene],
            'consensus_function': gene_annotation_map.get(gene, 'hypothetical protein')
        })
        module_counter += 1
        
    print(f"Defined a total of {len(accessory_modules)} accessory modules ({len(multi_gene_modules)} multi-gene, {len(single_genes)} single-gene).")
    return accessory_modules
