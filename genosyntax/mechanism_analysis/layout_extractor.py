import json
import pandas as pd
from tqdm import tqdm
from collections import defaultdict
import os
import re

def run_layout_extraction(output_dir, architecture_def, gene_order_db, gene_coord_db):
    """Main function to extract the true structural layouts of all high-quality genomes."""
    print("  -> Extracting genome layouts from anchor and slot records...")
    

    
    # Get the parent directory of `output_dir` 
    base_analysis_dir = os.path.dirname(os.path.abspath(output_dir.rstrip('/')))
    
    # Construct the correct path to slot_database.json
    slot_db_path = os.path.join(base_analysis_dir, "2_slot_analysis", "intermediate", "slot_database.json")


    try:
        with open(slot_db_path, 'r') as f:
            slot_db = json.load(f)
    except FileNotFoundError:
        print(f"  [ERROR] Cannot determine true layouts. Essential file 'slot_database.json' not found at the expected path: {slot_db_path}")
        return

    # --- STEP 1: Determine the set of PRESENT anchors in each genome ---
    genome_to_present_anchors = defaultdict(set)
    anchor_pattern = re.compile(r'(Anchor_\d+)')

    for genome_id, slots in slot_db.items():
        for slot_id in slots.keys():
            found_anchors = anchor_pattern.findall(slot_id)
            for anchor in found_anchors:
                genome_to_present_anchors[genome_id].add(anchor)

    # --- STEP 2: Get the theoretical anchor ORDER for each genome ---
    genome_to_order_template = {}
    macro_arch = architecture_def['macro_architectures']
    
    major_arch_group = macro_arch.get('major_architecture', {})
    if major_arch_group:
        for genome_id in major_arch_group.get('supporting_genomes', []):
            genome_to_order_template[genome_id] = major_arch_group.get('permutation', [])
            
    for minor_arch_group in macro_arch.get('minor_architectures', []):
        for genome_id in minor_arch_group.get('supporting_genomes', []):
            genome_to_order_template[genome_id] = minor_arch_group.get('permutation', [])

    # --- STEP 3: Combine PRESENCE and ORDER to define the TRUE LAYOUT ---
    true_layouts = {}
    for genome_id in genome_to_order_template.keys():
        order_template = genome_to_order_template.get(genome_id, [])
        present_anchors = genome_to_present_anchors.get(genome_id, set())
        true_layouts[genome_id] = [anchor for anchor in order_template if anchor in present_anchors]

    # --- STEP 4: Convert true layouts to physical coordinates ---
    layout_data = []
    print(f"  -> Processing {len(true_layouts)} genomes to extract physical layouts...")

    anchors = {a['anchor_id']: a for a in architecture_def['anchor_definitions']}
    
    genomes_to_process = sorted(true_layouts.keys())

    for genome_id in tqdm(genomes_to_process, desc="  Extracting Layouts", leave=False):
        current_layout = true_layouts.get(genome_id, [])
        if genome_id not in gene_coord_db or not gene_coord_db[genome_id]: continue
            
        genome_coords = gene_coord_db[genome_id]
        genome_length = max((c['end'] for c in genome_coords), default=0)
        
        located_anchors_with_coords = []
        for anchor_id in current_layout:
            skeleton_genes = set(anchors.get(anchor_id, {}).get('skeleton', []))
            coords_for_this_anchor = [c for c in genome_coords if c['group_id'] in skeleton_genes]
            
            if coords_for_this_anchor:
                start_bp = min(c['start'] for c in coords_for_this_anchor)
                end_bp = max(c['end'] for c in coords_for_this_anchor)
                located_anchors_with_coords.append({'id': anchor_id, 'start_bp': start_bp, 'end_bp': end_bp})
        
        located_anchors_with_coords.sort(key=lambda x: x['start_bp'])
        
        last_pos = 1
        # Upstream Slot
        if located_anchors_with_coords:
            first_anchor_start = located_anchors_with_coords[0]['start_bp']
            if first_anchor_start > last_pos:
                layout_data.append({'genome_id': genome_id, 'element_id': 'Slot_Upstream', 'element_type': 'Slot', 'start': last_pos, 'end': first_anchor_start - 1})
        
        # Anchors and Inter-Anchor Slots
        for i, anchor in enumerate(located_anchors_with_coords):
            layout_data.append({'genome_id': genome_id, 'element_id': anchor['id'], 'element_type': 'Anchor', 'start': anchor['start_bp'], 'end': anchor['end_bp']})
            
            if i < len(located_anchors_with_coords) - 1:
                current_anchor_end = anchor['end_bp']
                next_anchor_start = located_anchors_with_coords[i+1]['start_bp']
                
                if next_anchor_start > current_anchor_end + 1:
                    next_anchor_id = located_anchors_with_coords[i+1]['id']
                    slot_id = f"Slot({anchor['id']}->{next_anchor_id})"
                    layout_data.append({'genome_id': genome_id, 'element_id': slot_id, 'element_type': 'Slot', 'start': current_anchor_end + 1, 'end': next_anchor_start - 1})

        # Downstream Slot
        last_anchor_end = located_anchors_with_coords[-1]['end_bp'] if located_anchors_with_coords else 0
        if last_anchor_end < genome_length:
            layout_data.append({'genome_id': genome_id, 'element_id': 'Slot_Downstream', 'element_type': 'Slot', 'start': last_anchor_end + 1, 'end': genome_length})
            
    # --- STEP 5: Save the output CSV ---
    df_layout = pd.DataFrame(layout_data)
    if not df_layout.empty:
        df_layout['length'] = df_layout['end'] - df_layout['start'] + 1
    
    reports_dir = os.path.join(output_dir, "reports")
    os.makedirs(reports_dir, exist_ok=True)
    output_path = os.path.join(reports_dir, "all_genome_layouts.csv")
    df_layout.to_csv(output_path, index=False)
    
    print(f"  -> Genome layout data for {df_layout['genome_id'].nunique()} genomes saved to: {output_path}")
