import os
from tqdm import tqdm
import numpy as np
from scipy.stats import spearmanr
from itertools import combinations
from genosyntax.utils.file_ops import calculate_asi_for_gene, calculate_gc_content

def discover_anchors(gene_order_db, strict_core_genes, passed_genomes, 
                     tolerance, coverage_threshold, min_length):
    print("\n--- Task 3: Discovering Anchor Modules ---")
    if not passed_genomes: return []
    reference_genome_id = passed_genomes[0]
    reference_order = gene_order_db[reference_genome_id]
    candidates, current_block, accessory_count = [], [], 0
    for gene in reference_order:
        if gene in strict_core_genes:
            current_block.append(gene); accessory_count = 0
        else:
            accessory_count += 1
            if accessory_count > tolerance:
                if len(current_block) >= min_length: candidates.append(tuple(current_block))
                current_block = []
    if len(current_block) >= min_length: candidates.append(tuple(current_block))
    unique_candidates = sorted(list(set(candidates)), key=len, reverse=True)
    print(f"Identified {len(unique_candidates)} unique candidate skeletons from reference '{reference_genome_id}'.")
    print("Verifying candidates across all genomes...")
    skeleton_to_ref_pos = {sk: reference_order.index(sk[0]) for sk in unique_candidates}
    valid_anchors, num_genomes, min_genomes_required = [], len(passed_genomes), len(passed_genomes) * coverage_threshold
    for skeleton in tqdm(unique_candidates, desc="  Verifying skeletons"):
        count = sum(1 for genome_id in passed_genomes if _verify_skeleton_in_genome(skeleton, gene_order_db[genome_id], strict_core_genes, tolerance))
        if count >= min_genomes_required:
            valid_anchors.append({"skeleton": list(skeleton), "core_gene_count": len(skeleton), "coverage_count": count, "coverage_fraction": round(count / num_genomes, 4),"reference_position": skeleton_to_ref_pos[skeleton]})
    print(f"Found {len(valid_anchors)} anchors meeting the {coverage_threshold:.0%} coverage threshold.")
    return valid_anchors

def _verify_skeleton_in_genome(skeleton, genome_order, core_genes, tolerance):
    skel_ptr, skel_len, order_ptr, order_len = 0, len(skeleton), 0, len(genome_order)
    while order_ptr < order_len:
        try: start_pos = genome_order.index(skeleton[0], order_ptr)
        except ValueError: return False
        skel_ptr, genome_ptr, acc_count = 1, start_pos + 1, 0
        while skel_ptr < skel_len and genome_ptr < order_len:
            if genome_order[genome_ptr] == skeleton[skel_ptr]: skel_ptr += 1; acc_count = 0
            elif genome_order[genome_ptr] not in core_genes:
                acc_count += 1
                if acc_count > tolerance: break
            else: break;
            genome_ptr += 1
        if skel_ptr == skel_len: return True
        order_ptr = start_pos + 1
    return False


def characterize_anchors(anchors, protein_db, gene_order_db, gene_coord_db, 
                       genome_seqs, 
                       passed_genomes, temp_dir, boundary_window=200, tolerance=2):

    print("\n--- Task 4: Characterizing Anchor Modules ---")
    print("  Pre-calculating global background GC content...")
    global_gc = [calculate_gc_content(str(genome_seqs[gid].seq)) for gid in passed_genomes if gid in genome_seqs and str(genome_seqs[gid].seq)]
    background_gc_mean, background_gc_std = np.mean(global_gc) if global_gc else 0.0, np.std(global_gc) if global_gc else 1.0
    anchor_presence_map = {}
    anchors.sort(key=lambda x: x['reference_position'])
    all_core_genes = {gene for sk in anchors for gene in sk['skeleton']}
    for i, anchor in enumerate(anchors):
        anchor_id = f"Anchor_{i+1}"; anchor['anchor_id'] = anchor_id
        anchor_presence_map[anchor_id] = [gid for gid in passed_genomes if _verify_skeleton_in_genome(tuple(anchor['skeleton']), gene_order_db[gid], all_core_genes, tolerance)]
    for anchor in tqdm(anchors, desc="  Characterizing anchors"):
        skeleton, present_genomes = anchor['skeleton'], anchor_presence_map[anchor['anchor_id']]
        all_gene_asis = [asi for core_gene in skeleton if (asi := calculate_asi_for_gene([s for gid in present_genomes if (s := protein_db.get(core_gene, {}).get(gid))], temp_dir)) is not None]
        anchor['asi'] = round(np.mean(all_gene_asis), 4) if all_gene_asis else 0.0
        rank_data = [ranks for gid in present_genomes if len(ranks := [gene_order_db[gid].index(g) for g in skeleton if g in gene_order_db[gid]]) == len(skeleton)]
        if len(rank_data) < 2: anchor['aoi'] = 1.0
        else:
            rank_matrix = np.array(rank_data)
            correlations = [spearmanr(rank_matrix[:, i], rank_matrix[:, j])[0] for i, j in combinations(range(len(skeleton)), 2)]
            anchor['aoi'] = round(np.mean([c for c in correlations if not np.isnan(c)]), 4) if correlations else 1.0
        boundary_gc_zscores = []
        first_gene, last_gene = skeleton[0], skeleton[-1]
        for genome_id in present_genomes:
            if genome_id not in genome_seqs or genome_id not in gene_coord_db: continue
            genome_sequence, coords = str(genome_seqs[genome_id].seq), gene_coord_db[genome_id]
            first_coord, last_coord = next((g for g in coords if g['group_id'] == first_gene), None), next((g for g in coords if g['group_id'] == last_gene), None)
            if not first_coord or not last_coord: continue
            up_end, up_start = first_coord['start'] - 1, max(0, first_coord['start'] - 1 - boundary_window)
            if up_end > up_start: boundary_gc_zscores.append((calculate_gc_content(genome_sequence[up_start:up_end]) - background_gc_mean) / background_gc_std if background_gc_std > 0 else 0)
            down_start, down_end = last_coord['end'], min(len(genome_sequence), last_coord['end'] + boundary_window)
            if down_end > down_start: boundary_gc_zscores.append((calculate_gc_content(genome_sequence[down_start:down_end]) - background_gc_mean) / background_gc_std if background_gc_std > 0 else 0)
        anchor['boundary_signals'] = {"window_size_bp": boundary_window, "avg_gc_enrichment_zscore": round(np.mean(boundary_gc_zscores), 4) if boundary_gc_zscores else 0.0}
    return anchors
