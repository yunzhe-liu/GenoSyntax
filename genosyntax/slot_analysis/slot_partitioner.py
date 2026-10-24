import json
from collections import defaultdict, Counter
from tqdm import tqdm
import pandas as pd

def _find_anchor_location(skeleton, genome_order, all_core_genes, tolerance):
    """
    Find the start and end indices of an anchor
    skeleton within a genome's gene order.
    """
    skel_len = len(skeleton)
    order_len = len(genome_order)
    try:
        first_gene_indices = [i for i, g in enumerate(genome_order) if g == skeleton[0]]
    except ValueError:
        return None, None
    for start_pos in first_gene_indices:
        skel_ptr, genome_ptr, acc_count = 1, start_pos + 1, 0
        found_indices = [start_pos]
        while skel_ptr < skel_len and genome_ptr < order_len:
            if genome_order[genome_ptr] == skeleton[skel_ptr]:
                found_indices.append(genome_ptr)
                skel_ptr += 1
                acc_count = 0
            elif genome_order[genome_ptr] not in all_core_genes:
                acc_count += 1
                if acc_count > tolerance: break
            else: break
            genome_ptr += 1
        if skel_ptr == skel_len:
            return min(found_indices), max(found_indices)
    return None, None


def partition_slot_space(architecture_def, gene_order_db, gene_coord_db):
    """
    Partitions the slot space based on the provided architecture definition.
    """
    print("\n--- Task 1: Partitioning Slot Space ---")
    
    slot_database = defaultdict(dict)
    anchor_defs = {a['anchor_id']: a for a in architecture_def['anchor_definitions']}
    
    try:
        pangenome_strat = architecture_def['metadata']['pangenome_stratification']
        all_core_genes = set(pangenome_strat.get('strict_core', []))
        params = architecture_def.get('metadata', {}).get('analysis_parameters', {})
        tolerance = params.get('tolerance', 2)
        qc_summary = architecture_def.get('metadata', {}).get('qc_summary', {})
        passed_genomes = qc_summary.get('passed_genomes_list', [])
    except KeyError as e:
        print(f"[FATAL] A required key was not found in the architecture_definition.json metadata: {e}")
        return {}
        
    if not passed_genomes:
        print("[ERROR] No high-quality genomes found in the architecture definition. Aborting.")
        return {}

    for genome_id in tqdm(passed_genomes, desc="Partitioning genomes"):
        gene_order = gene_order_db.get(genome_id)
        coords = gene_coord_db.get(genome_id)
        if not gene_order or not coords: continue

        located_anchors = []
        for anchor_id, anchor_info in anchor_defs.items():
            start_idx, end_idx = _find_anchor_location(
                tuple(anchor_info['skeleton']), gene_order, all_core_genes, tolerance
            )
            if start_idx is not None:
                located_anchors.append({'id': anchor_id, 'start_idx': start_idx, 'end_idx': end_idx})
        
        if not located_anchors: continue
        located_anchors.sort(key=lambda x: x['start_idx'])
        
        def get_slot_size_bp(start_gene_idx, end_gene_idx):
            # This function requires the 'coords' list, hence the check above.
            if start_gene_idx >= len(coords) or end_gene_idx > len(coords): return 0
            if start_gene_idx >= end_gene_idx: return 0
            end_bp = coords[start_gene_idx-1]['end'] if start_gene_idx > 0 else 0
            start_bp = coords[end_gene_idx]['start'] if end_gene_idx < len(coords) else (coords[-1]['end'] if coords else 0)
            return abs(start_bp - end_bp)

        is_consistent = qc_summary.get('orientation_analysis', {}).get('is_consistently_oriented', False)
        
        if is_consistent:
            first_anchor = located_anchors[0]
            if first_anchor['start_idx'] > 0:
                slot_id_high = "Slot(Upstream)"
                slot_content = gene_order[0:first_anchor['start_idx']]
                slot_database[genome_id][slot_id_high] = {
                    'slot_id_high_res': slot_id_high, 'slot_id_low_res': 'Slot(Upstream)',
                    'slot_type': 'Upstream', 'content': slot_content,
                    'slot_size_genes': len(slot_content),
                    'slot_size_bp': get_slot_size_bp(0, first_anchor['start_idx'])
                }
            last_anchor = located_anchors[-1]
            if last_anchor['end_idx'] < len(gene_order) - 1:
                slot_id_high = "Slot(Downstream)"
                slot_content = gene_order[last_anchor['end_idx']+1:]
                slot_database[genome_id][slot_id_high] = {
                    'slot_id_high_res': slot_id_high, 'slot_id_low_res': 'Slot(Downstream)',
                    'slot_type': 'Downstream', 'content': slot_content,
                    'slot_size_genes': len(slot_content),
                    'slot_size_bp': get_slot_size_bp(last_anchor['end_idx']+1, len(gene_order))
                }
                
        for i in range(len(located_anchors) - 1):
            anchor1, anchor2 = located_anchors[i], located_anchors[i+1]
            start_idx, end_idx = anchor1['end_idx'] + 1, anchor2['start_idx']
            if start_idx >= end_idx: continue
            
            slot_id_high = f"Slot({anchor1['id']} -> {anchor2['id']})"
            slot_content = gene_order[start_idx:end_idx]
            slot_database[genome_id][slot_id_high] = {
                'slot_id_high_res': slot_id_high, 'slot_id_low_res': slot_id_high,
                'slot_type': 'Inter-Anchor', 'content': slot_content,
                'slot_size_genes': len(slot_content),
                'slot_size_bp': get_slot_size_bp(start_idx, end_idx)
            }
            
        for anchor in located_anchors:
            anchor_id = anchor['id']
            skeleton = set(anchor_defs[anchor_id]['skeleton'])
            
            skeleton_indices_in_span = [i for i in range(anchor['start_idx'], anchor['end_idx'] + 1) if gene_order[i] in skeleton]
            skeleton_indices_sorted = sorted(skeleton_indices_in_span)
            
            for i in range(len(skeleton_indices_sorted) - 1):
                start_idx, end_idx = skeleton_indices_sorted[i] + 1, skeleton_indices_sorted[i+1]
                if start_idx >= end_idx: continue
                
                gene1 = gene_order[skeleton_indices_sorted[i]]
                gene2 = gene_order[skeleton_indices_sorted[i+1]]
                
                slot_id_high = f"Slot({gene1} -> {gene2}, in: {anchor_id})"
                slot_id_low = f"Slot(in: {anchor_id})"
                slot_content = gene_order[start_idx:end_idx]
                
                slot_database[genome_id][slot_id_high] = {
                    'slot_id_high_res': slot_id_high, 'slot_id_low_res': slot_id_low,
                    'slot_type': 'Intra-Anchor', 'content': slot_content,
                    'slot_size_genes': len(slot_content),
                    'slot_size_bp': get_slot_size_bp(start_idx, end_idx)
                }
    
    print("Slot space partitioning complete.")
    return slot_database


def analyze_grey_zone_genes(slot_database, pangenome_strat, anchor_defs):
    """Analyzes the distribution of 'Soft Core' and 'Orphan Core' genes within the slot space."""
    print("\n--- Task 2: Analyzing 'Grey Zone' Gene Distribution ---")
    
    soft_core_genes = set(pangenome_strat.get('soft_core', []))
    all_anchor_core_genes = {gene for anchor in anchor_defs for gene in anchor.get('skeleton', [])}
    orphan_core_genes = set(pangenome_strat.get('strict_core', [])) - all_anchor_core_genes
    
    soft_core_dist = defaultdict(Counter)
    orphan_core_dist = defaultdict(Counter)
    
    for genome_id, slots in slot_database.items():
        for slot_id, slot_info in slots.items():
            for gene in slot_info.get('content', []):
                if gene in soft_core_genes:
                    soft_core_dist[gene][slot_info['slot_id_low_res']] += 1
                elif gene in orphan_core_genes:
                    orphan_core_dist[gene][slot_info['slot_id_low_res']] += 1
                    
    def format_distribution_report(dist_dict, gene_type):
        rows = []
        for gene, counter in dist_dict.items():
            for slot_id, count in counter.items():
                rows.append({
                    'gene_id': gene,
                    'gene_type': gene_type,
                    'found_in_slot_id': slot_id,
                    'observation_count': count
                })
        
        # If no rows were created, return an empty DataFrame WITH the correct columns defined.
        if not rows:
            return pd.DataFrame(columns=['gene_id', 'gene_type', 'found_in_slot_id', 'observation_count'])
        return pd.DataFrame(rows)

    df_soft_core = format_distribution_report(soft_core_dist, 'Soft_Core')
    df_orphan_core = format_distribution_report(orphan_core_dist, 'Orphan_Core')
    
    print(f"Found {df_orphan_core['gene_id'].nunique()} Orphan Core and {df_soft_core['gene_id'].nunique()} Soft Core genes in slots.")
    return df_orphan_core, df_soft_core
