# /cluster/home/yunzliu/Project/GenoSynax_test/slot_analysis/src/dynamics_analyzer.py

import pandas as pd
from collections import defaultdict
from tqdm import tqdm

def analyze_module_dynamics(
    accessory_modules, 
    slot_database, 
    pangenome_strat,
    passed_genomes,
    spi_resolution='low'
):
    """
    Analyzes the structural dynamics of accessory modules (AMs).
    VERSION 2 : Implements a strict definition for module presence in a slot.
    """
    print("\n--- Task 5: Analyzing Modular Dynamics ---")
    print(f"  Using SPI resolution: '{spi_resolution}'")

    module_location_map = defaultdict(set) # {module_id: set_of_observed_slot_keys}
    module_occurrence_map = defaultdict(set) # {module_id: set_of_genome_ids}
    
    # Create a reverse map from gene to module for quick lookup
    gene_to_module_map = {}
    for module in accessory_modules:
        for gene in module['members']:
            gene_to_module_map[gene] = module['module_id']

    # Pre-map modules to their member sets for efficient checking
    module_member_sets = {m['module_id']: set(m['members']) for m in accessory_modules}

    print("  Mapping module locations with strict presence definition...")
    for genome_id, slots in tqdm(slot_database.items(), desc="  Mapping module locations"):
        for slot_id, slot_info in slots.items():
            slot_content_set = set(slot_info['content'])
            slot_key = slot_info['slot_id_low_res'] if spi_resolution == 'low' else slot_info['slot_id_high_res']
            
            
            # Find all unique modules that have at least one member in the slot
            candidate_module_ids = {gene_to_module_map[g] for g in slot_content_set if g in gene_to_module_map}
            
            for module_id in candidate_module_ids:
                # A module is considered present in a slot ONLY IF
                # all of its members are found within that slot's content.
                if module_member_sets[module_id].issubset(slot_content_set):
                    module_location_map[module_id].add(slot_key)
                    module_occurrence_map[module_id].add(genome_id)

    print("\n  Calculating SPI and summarizing module dynamics...")
    results = []
    num_total_genomes = len(passed_genomes)

    for module in tqdm(accessory_modules, desc="  Analyzing modules"):
        module_id = module['module_id']
        
        # Get data from our new, correct maps
        observed_slots = module_location_map.get(module_id, set())
        num_genomes_with_module = len(module_occurrence_map.get(module_id, set()))
        
        frequency = num_genomes_with_module / num_total_genomes if num_total_genomes > 0 else 0
        spi = len(observed_slots)
        
        behavioral_class = "N/A"
        if spi > 1: behavioral_class = "Multi-Position"
        elif spi == 1:
            if frequency > 0.5: behavioral_class = "Structurally Conserved"
            else: behavioral_class = "Niche-Specific"
        elif spi == 0:
             behavioral_class = "Not Observed in Slots"

        results.append({
            'module_id': module['module_id'],
            'module_type': module['module_type'],
            'num_genes': module['num_genes'],
            'members': ", ".join(module['members']),
            'consensus_function': module['consensus_function'],
            'frequency': round(frequency, 4),
            'occurrence_count': num_genomes_with_module,
            'spi': spi,
            'behavioral_class': behavioral_class,
            'observed_slots': "; ".join(sorted(list(observed_slots))) if observed_slots else "None"
        })

    report_df = pd.DataFrame(results).sort_values(by=['spi', 'frequency'], ascending=[False, False])
    
    print("Modular dynamics analysis complete.")
    return report_df
