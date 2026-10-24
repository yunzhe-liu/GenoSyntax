import os
import json
import pandas as pd

from genosyntax.define_architecture.data_loader import load_panaroo_data, load_gff_data, load_protein_sequences, load_genome_sequences
from genosyntax.define_architecture.quality_control import run_qc
from genosyntax.define_architecture.anchor_discovery import discover_anchors, characterize_anchors
from genosyntax.define_architecture.functional_annotator import annotate_anchor_functions
from genosyntax.define_architecture.macro_analysis import analyze_macro_architecture
from genosyntax.utils.file_ops import save_json

def run(args):
    """Entrypoint function for the 'define_architecture' module."""
    # --- Setup Directories ---
    output_dir = args.output_dir
    intermediate_dir = os.path.join(output_dir, "intermediate")
    reports_dir = os.path.join(output_dir, "reports")
    os.makedirs(intermediate_dir, exist_ok=True)
    os.makedirs(reports_dir, exist_ok=True)

    if args.dry_run:
        return True

    print("====== GenoSyntax: Running Module 'define_architecture' ======")
    
    # --- Step 1 & 2 : Data Loading & QC  ---
    genome_ids, stratified_genes, locus_to_group_map = load_panaroo_data(
        args.panaroo_file,
        strict_core_threshold=args.strict_core_threshold,
        soft_core_threshold=args.soft_core_threshold
    )

    gene_order_db, gene_coord_db = load_gff_data(args.pharokka_dir, locus_to_group_map, genome_ids)
    df_gpa = pd.read_csv(args.panaroo_file, low_memory=False)
    gene_annotation_map = pd.Series(df_gpa['Annotation'].values, index=df_gpa['Gene']).to_dict()
    protein_db = load_protein_sequences(args.faa_dir, locus_to_group_map)
    genome_seqs = load_genome_sequences(args.fasta_dir, genome_ids)
    passed_genomes, qc_report = run_qc(
        genome_ids, args.fasta_dir, gene_order_db, gene_coord_db, 
        max_contigs=getattr(args, 'max_contigs', 1),
        orientation_threshold=getattr(args, 'orientation_consistency_threshold', 0.9)
    )

    if not passed_genomes:
        save_json(qc_report, os.path.join(reports_dir, 'qc_report.json'))
        print('[ERROR] No genomes passed QC; see reports/qc_report.json.')
        return False

    # Step 3: Anchor Discovery
    params_da = {'tolerance': args.tolerance, 'coverage_threshold': args.coverage_threshold, 'min_length': args.min_length}
    anchors_basic = discover_anchors(gene_order_db, stratified_genes['strict_core'], passed_genomes, **params_da)
    
    # Check if any anchors were found. If not, we stop further analysis but create a valid output file.
    if not anchors_basic:
        print('[INFO] No anchors satisfy the current criteria; saving an empty architecture report.')

        # Create an empty but valid output structure
        final_output = {
            "metadata": {
                "analysis_parameters": vars(args), 
                "qc_summary": qc_report,
                "pangenome_stratification": stratified_genes,
                "status_message": "No anchors found. Analysis halted."
            },
            "anchor_definitions": [],
            "macro_architectures": {}
        }
    else:
        # If anchors were found, proceed with the normal workflow
        anchors_characterized = characterize_anchors(anchors_basic, protein_db, gene_order_db, gene_coord_db, genome_seqs, passed_genomes, intermediate_dir, tolerance=args.tolerance)
        anchors_functional = annotate_anchor_functions(anchors_characterized, gene_annotation_map)
        macro_report = analyze_macro_architecture(anchors_functional, gene_order_db, passed_genomes, tolerance=args.tolerance)
        
        final_output = {
            "metadata": {
                "analysis_parameters": vars(args), 
                "qc_summary": qc_report,
                "pangenome_stratification": stratified_genes
            },
            "anchor_definitions": anchors_functional,
            "macro_architectures": macro_report
        }

    # Save intermediate files regardless, as they can be useful for diagnostics
    print("\n[*] Saving intermediate databases for downstream modules or diagnostics...")
    save_json(stratified_genes, os.path.join(intermediate_dir, "pangenome_stratification.json"))
    save_json(locus_to_group_map, os.path.join(intermediate_dir, "locus_to_group_map.json"))
    save_json(gene_order_db, os.path.join(intermediate_dir, "gene_order_database.json"))
    save_json(gene_coord_db, os.path.join(intermediate_dir, "gene_coordinate_database.json"))
    save_json(gene_annotation_map, os.path.join(intermediate_dir, "gene_annotation_map.json"))
    print("[*] Intermediate files saved successfully.")

    # Save the final (possibly empty) architecture definition
    final_output_path = os.path.join(output_dir, "architecture_definition.json")
    save_json(final_output, final_output_path)
    
    print(f"\n====== 'define_architecture' Module Finished Successfully! ======")
    print(f"[*] Final output saved to: {final_output_path}")
    return True
