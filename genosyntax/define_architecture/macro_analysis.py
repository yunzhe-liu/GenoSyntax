# /GenoSynax_test/genosyntax/define_architecture/macro_analysis.py

from collections import Counter, defaultdict
from scipy.stats import entropy

# --- Robust Import ---
# It is critical to use the same anchor verification logic as the rest of the pipeline.
try:
    # This works when running as part of the genosyntax package.
    from .anchor_discovery import _verify_skeleton_in_genome
except ImportError:
    # This is a fallback for standalone testing of the script.
    from anchor_discovery import _verify_skeleton_in_genome

# --- Helper Functions ---

def _find_anchor_start_pos(skeleton, genome_order):
    """
    Helper function to find the starting position of an anchor's skeleton.
    Uses the position of the first gene in the skeleton as a proxy for sorting.
    Returns -1 if the first gene is not found, ensuring stability.
    """
    try:
        return genome_order.index(skeleton[0])
    except ValueError:
        return -1

def is_subsequence(child_perm, parent_perm):
    """
    Checks if one tuple/list ('child_perm') is an ordered subsequence of another ('parent_perm').
    This is the core of the flexible grouping logic.
    e.g., is_subsequence(('A', 'C'), ('A', 'B', 'C')) -> True
    e.g., is_subsequence(('C', 'A'), ('A', 'B', 'C')) -> False
    """
    it = iter(parent_perm)
    return all(element in it for element in child_perm)

# --- Main Analysis Function ---

def analyze_macro_architecture(anchors, gene_order_db, passed_genomes, tolerance=2):
    """
    Analyzes the permutation order of anchor modules across all high-quality genomes.

    VERSION 3:
    This version implements a robust, two-stage process to correctly classify genomic architectures:
    1.  It first uses a strict verification function (`_verify_skeleton_in_genome`) to accurately
        identify the presence of structurally complete anchors in each genome.
    2.  It then groups genomes using an "Ordered Subset" (subsequence) matching algorithm.
        This ensures that genomes with significant anchor loss (e.g., containing only 2 of 10 anchors)
        are correctly classified as part of a major architecture, provided their remaining anchors
        maintain the correct relative order. This prevents incorrect fragmentation into multiple
        "chaotic" architectures.

    Args:
        anchors (list): List of anchor definition dictionaries from previous steps.
        gene_order_db (dict): Maps genome_id to its ordered list of gene group IDs.
        passed_genomes (list): List of genome IDs that passed quality control.
        tolerance (int): The number of non-core genes tolerated within an anchor skeleton.

    Returns:
        dict: A report summarizing the macro-architectural landscape, including major/minor
              architectures and global permutation entropy.
    """
    print("\n--- Task 5: Analyzing Macro-Architecture ---")

    if not anchors:
        print("  [Warning] No anchor definitions provided. Skipping macro-analysis.")
        return {}

    # --- STEP 1: Robustly locate all valid anchors and their order in each genome ---
    print("  Step 1: Accurately locating structurally complete anchors in each genome...")
    genome_permutations = {}
    anchor_map = {a['anchor_id']: a for a in anchors}
    all_core_genes = {gene for anchor in anchors for gene in anchor['skeleton']}

    for genome_id in passed_genomes:
        gene_order = gene_order_db.get(genome_id)
        if not gene_order: continue

        located_anchors = []
        for anchor_id, anchor_info in anchor_map.items():
            skeleton = tuple(anchor_info['skeleton'])
            
            # Use the robust verification function to ensure the entire anchor is present.
            if _verify_skeleton_in_genome(skeleton, gene_order, all_core_genes, tolerance):
                start_pos = _find_anchor_start_pos(skeleton, gene_order)
                if start_pos != -1:
                    located_anchors.append({'pos': start_pos, 'id': anchor_id})
        
        if located_anchors:
            # Sort the found anchors by their physical position in the genome.
            located_anchors.sort(key=lambda x: x['pos'])
            genome_permutations[genome_id] = tuple(a['id'] for a in located_anchors)

    if not genome_permutations:
        print("  [Warning] No valid anchor permutations could be determined for any genome.")
        return {}

    # --- STEP 2: Group architectures based on the "Ordered Subset" principle ---
    print("  Step 2: Grouping architectures via ordered subset matching...")
    
    unique_perms = sorted(Counter(genome_permutations.values()).keys(), key=len, reverse=True)
    
    # This map will store the "ultimate parent" for every observed permutation.
    parent_map = {p: p for p in unique_perms} # Initially, each permutation is its own parent.

    # Two-pass approach to resolve the full hierarchy of parent-child relationships.
    for p_parent in unique_perms:
        for p_child in unique_perms:
            if len(p_child) >= len(p_parent): continue # A child cannot be longer than its parent.
            if p_child == p_parent: continue

            # If the child is a subsequence of the parent, map the child to the parent.
            if is_subsequence(p_child, parent_map[p_parent]):
                 parent_map[p_child] = parent_map[p_parent]

    # Second pass for path compression: ensure all mappings point to the ultimate root parent.
    for perm in unique_perms:
        parent = parent_map[perm]
        while parent != parent_map[parent]:
            parent = parent_map[parent]
        parent_map[perm] = parent
        
    # --- STEP 3: Assemble the final report from the grouped architectural families ---
    print(f"  Step 3: Assembling final report for {len(set(parent_map.values()))} distinct architectural families...")
    
    architecture_families = defaultdict(list)
    for genome_id, perm in genome_permutations.items():
        final_parent_perm = parent_map[perm]
        architecture_families[final_parent_perm].append(genome_id)

    final_groups = []
    for parent_perm, genomes in architecture_families.items():
        final_groups.append({
            "permutation": list(parent_perm),
            "supporting_genomes": genomes,
        })

    # Sort final architecture groups by the number of supporting genomes.
    sorted_groups = sorted(final_groups, key=lambda x: len(x['supporting_genomes']), reverse=True)
    
    # --- STEP 4: Calculate global entropy and finalize the report ---
    major_architecture = sorted_groups[0] if sorted_groups else {}
    minor_architectures = sorted_groups[1:] if len(sorted_groups) > 1 else []

    total_genomes_with_arch = len(genome_permutations)
    frequencies = [len(g['supporting_genomes']) / total_genomes_with_arch for g in sorted_groups]
    global_entropy = entropy(frequencies, base=2) if frequencies else 0.0

    # Add frequency information to the report.
    if major_architecture:
        major_architecture['frequency'] = round(len(major_architecture['supporting_genomes']) / total_genomes_with_arch, 4)
    for group in minor_architectures:
        group['frequency'] = round(len(group['supporting_genomes']) / total_genomes_with_arch, 4)

    macro_report = {
        "global_permutation_entropy": round(global_entropy, 4),
        "total_architectures_found": len(sorted_groups),
        "major_architecture": major_architecture,
        "minor_architectures": minor_architectures
    }
    
    print("Macro-architecture analysis complete.")
    return macro_report
