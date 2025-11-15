from collections import Counter

# Define a dictionary of keywords for functional super-classes.
# This is the "knowledge base" of our annotation engine.
# It should be easily configurable.
FUNCTIONAL_MAP = {
    "Virion Morphogenesis & Infection": [
        'tail', 'baseplate', 'capsid', 'head', 'portal', 'terminase', 
        'lysis', 'holin', 'endolysin', 'sheath', 'fiber'
    ],
    "Host Takeover & Genome Replication": [
        'dna polymerase', 'helicase', 'primase', 'ligase', 'replication',
        'rna polymerase', 'transcription', 'sigma factor'
    ],
    "Genome Maintenance & Host Interaction": [
        'recombination', 'exonuclease', 'endonuclease', 'repair', 
        'nucleotide', 'synthesis', 'kinase', 'synthase', 'reductase',
        'metabolism', 'folate'
    ],
    "Unknown & Other": [
        'hypothetical', 'unknown', 'uncharacterized'
    ]
}

def annotate_anchor_functions(anchors, gene_annotation_map):
    """
    Assigns a high-level functional category to each anchor module based on the
    functions of its member genes.

    Args:
        anchors (list): The list of characterized anchor definition dictionaries.
        gene_annotation_map (dict): A map from group_id to functional description.

    Returns:
        list: The same list of anchors, with 'functional_class' and 
              'function_profile' fields added to each anchor dictionary.
    """
    print("\n--- Task 6: Annotating Anchor Functions ---")
    
    for anchor in anchors:
        skeleton = anchor['skeleton']
        function_counts = Counter()

        # Get annotations for all genes in the anchor's skeleton
        gene_functions = [gene_annotation_map.get(gene_id, 'unknown') for gene_id in skeleton]
        
        # Classify each gene's function
        for func_desc in gene_functions:
            assigned_class = "Unknown & Other" # Default
            for category, keywords in FUNCTIONAL_MAP.items():
                # Use lower case for robust matching
                if any(keyword in func_desc.lower() for keyword in keywords):
                    assigned_class = category
                    break
            function_counts[assigned_class] += 1
            
        # Determine the dominant functional class for the anchor
        # We prioritize non-unknown functions
        non_unknown_counts = {k: v for k, v in function_counts.items() if k != "Unknown & Other"}
        
        if non_unknown_counts:
            # The class is the one with the most member genes
            dominant_class = max(non_unknown_counts, key=non_unknown_counts.get)
        else:
            dominant_class = "Unknown & Other"
            
        anchor['functional_class'] = dominant_class
        anchor['function_profile'] = dict(function_counts)

    print("Functional annotation complete.")
    return anchors
