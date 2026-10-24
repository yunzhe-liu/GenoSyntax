import os
import pandas as pd
from collections import defaultdict

def run_adjacency_analysis(output_dir, accessory_modules, gene_annotation_map):
    """Export annotations and membership tables for multi-gene modules."""
    analysis_dir = os.path.join(output_dir, "adjacency_module_analysis")
    os.makedirs(analysis_dir, exist_ok=True)
    
    adjacency_modules = [m for m in accessory_modules if m.get('module_type') == 'multi-gene']
    print(f"  -> Found {len(adjacency_modules)} adjacency (multi-gene) modules to analyze.")
    
    if not adjacency_modules:
        print("  [INFO] No adjacency modules were found in the input. Analysis complete.")
        # Create empty files to show the task ran
        open(os.path.join(analysis_dir, "adjacency_modules_narrative_report.txt"), 'w').write("No adjacency modules found.\n")
        pd.DataFrame().to_csv(os.path.join(analysis_dir, "adjacency_modules_summary.csv"), index=False)
        return

    report_sections = defaultdict(list)
    node_data, edge_data, heatmap_data, summary_data = [], [], [], []
    non_informative_keywords = ["hypothetical protein", "unknown function"]

    for module in sorted(adjacency_modules, key=lambda x: x.get('module_id', '')):
        module_id = module.get('module_id', 'Unknown_ID')
        members = module.get('members', [])
        num_genes = len(members)
        
        informative_annot_count = 0
        gene_details = [] 

        for i, gene_id in enumerate(members):
            annotation = gene_annotation_map.get(gene_id, "unknown function")
            is_informative = not any(keyword in annotation.lower() for keyword in non_informative_keywords)
            if is_informative:
                informative_annot_count += 1
            
            gene_details.append({'id': gene_id, 'annotation': annotation})
            heatmap_data.append({'module_id': module_id, 'member_index': i, 'gene_id': gene_id, 'function_description': annotation})
            node_data.append({'NodeID': gene_id, 'ModuleID': module_id, 'FunctionDescription': annotation})

        if informative_annot_count == num_genes and num_genes > 0:
            category, confidence = "Category A: All Members Annotated", "Complete"
        elif informative_annot_count > 0:
            category, confidence = "Category B: Some Members Annotated", "Partial"
        else:
            category, confidence = "Category C: No Informative Annotations", "None"

        narrative = f"### Module: {module_id} ({num_genes} genes)\n"
        narrative += f"Annotation Coverage: {confidence} ({informative_annot_count}/{num_genes} informative)\n\n"
        narrative += "| Gene ID    | Function Description (from Panaroo)\n"
        narrative += "|------------|--------------------------------------------------\n"
        for gene in gene_details:
            narrative += f"| {gene['id']:<10} | {gene['annotation']}\n"
        
        report_sections[category].append(narrative)

        summary_data.append({'module_id': module_id, 'num_genes': num_genes, 'annotation_coverage_category': category, 'informative_gene_count': informative_annot_count, 'member_ids': ", ".join(members)})

        for gene_id in members:
            edge_data.append({'Source': module_id, 'Target': gene_id,
                              'InteractionType': 'ModuleMembership', 'ModuleID': module_id})
        node_data.append({'NodeID': module_id, 'ModuleID': module_id,
                          'FunctionDescription': 'Accessory module'})

    final_report = "====== Adjacency Module Narrative Functional Report ======\n\n"
    final_report += "This report categorizes multi-gene adjacency modules based on the coverage of their functional annotations, sourced from Panaroo.\n\n"
    category_order = ["Category A: All Members Annotated", "Category B: Some Members Annotated", "Category C: No Informative Annotations"]

    for cat_name in category_order:
        if cat_name in report_sections:
            final_report += f"===========================================================\n"
            final_report += f"## {cat_name}\n"
            final_report += f"===========================================================\n\n"
            separator = "\n---\n\n"
            final_report += separator.join(report_sections[cat_name])
            final_report += "\n\n"

    with open(os.path.join(analysis_dir, "adjacency_modules_narrative_report.txt"), 'w') as f: f.write(final_report)
    print(f"  -> Detailed narrative report saved.")
    
    pd.DataFrame(summary_data).to_csv(os.path.join(analysis_dir, "adjacency_modules_summary.csv"), index=False)
    print(f"  -> CSV summary report saved.")
    
    pd.DataFrame(heatmap_data).to_csv(os.path.join(analysis_dir, "adjacency_module_heatmap_data.csv"), index=False)
    print(f"  -> Heatmap data saved.")
    
    pd.DataFrame(node_data).to_csv(os.path.join(analysis_dir, "adjacency_nodes.csv"), index=False)
    print(f"  -> Cytoscape node table saved.")
    
    pd.DataFrame(edge_data).to_csv(os.path.join(analysis_dir, "adjacency_edges.csv"), index=False)
    print(f"  -> Cytoscape edge table saved.")
    print(f"  -> All outputs are in: {analysis_dir}")
