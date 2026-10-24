import os
import json

from genosyntax.utils.file_ops import load_json, save_json
from genosyntax.slot_analysis.slot_partitioner import partition_slot_space, analyze_grey_zone_genes
from genosyntax.slot_analysis.module_discoverer import discover_multi_gene_modules, define_accessory_modules
from genosyntax.slot_analysis.dynamics_analyzer import analyze_module_dynamics

def run(args):
    """Entrypoint function for the 'analyze_slots' module."""
    if args.dry_run:
        return True

    output_dir = args.output_dir
    intermediate_dir = os.path.join(output_dir, "intermediate")
    reports_dir = os.path.join(output_dir, "reports")
    os.makedirs(intermediate_dir, exist_ok=True)
    os.makedirs(reports_dir, exist_ok=True)

    print("====== GenoSyntax: Running Module 'analyze_slots' ======")

    # Step 0: Load data
    architecture_def = load_json(args.architecture_file)
    if not architecture_def:
        print(f"[FATAL] Could not load architecture definition file from: {args.architecture_file}")
        return False # Indicate failure for pipeline

    # Check if the previous step actually found any anchors.
    if not architecture_def.get("anchor_definitions"):
        print("\n[INFO] The architecture definition file contains no anchor definitions.")
        print("This indicates that the 'define_architecture' step did not find a stable backbone.")
        print("[ACTION] Skipping 'analyze_slots' module as it cannot run without anchors.")
        print(f"\n====== 'analyze_slots' Module Skipped. ======")
        return True # Return True because this is an expected outcome, not a crash.
        
    gene_order_db = load_json(os.path.join(args.intermediate_dir_prev, 'gene_order_database.json'))
    gene_coord_db = load_json(os.path.join(args.intermediate_dir_prev, 'gene_coordinate_database.json'))
    pangenome_strat = load_json(os.path.join(args.intermediate_dir_prev, 'pangenome_stratification.json'))
    gene_annotation_map = load_json(os.path.join(args.intermediate_dir_prev, 'gene_annotation_map.json'))

    if not all([gene_order_db, gene_coord_db, pangenome_strat, gene_annotation_map]):
        print(f"[FATAL] Could not load one or more database files from: {args.intermediate_dir_prev}")
        return False
        
    # Steps 1-5: The actual analysis logic
    slot_database = partition_slot_space(architecture_def, gene_order_db, gene_coord_db)
    df_orphan, df_soft = analyze_grey_zone_genes(slot_database, pangenome_strat, architecture_def['anchor_definitions'])
    passed_genomes = architecture_def['metadata']['qc_summary']['passed_genomes_list']
    multi_gene_modules = discover_multi_gene_modules(
        pangenome_strat, gene_order_db, gene_coord_db, passed_genomes,
        phi_threshold=args.phi_threshold,
        p_value_threshold=args.p_value_threshold,
        adjacency_distance=args.adjacency_distance,
        support_threshold=args.support_threshold,
        seed=getattr(args, 'seed', 0)
    )
    
    accessory_modules = define_accessory_modules(pangenome_strat, multi_gene_modules, gene_annotation_map)
    module_dynamics_df = analyze_module_dynamics(accessory_modules, slot_database, pangenome_strat, passed_genomes, spi_resolution=args.spi_resolution)

    # Save outputs
    save_json(slot_database, os.path.join(intermediate_dir, "slot_database.json"))
    save_json(accessory_modules, os.path.join(intermediate_dir, "accessory_module_definitions.json"))
    if not df_orphan.empty: df_orphan.to_csv(os.path.join(reports_dir, "orphan_core_gene_distribution.csv"), index=False)
    if not df_soft.empty: df_soft.to_csv(os.path.join(reports_dir, "soft_core_gene_distribution.csv"), index=False)
    module_dynamics_df.to_csv(os.path.join(reports_dir, "module_dynamics_report.csv"), index=False)

    print(f"\n====== 'analyze_slots' Module Finished Successfully! ======")
    print(f"[*] Final reports saved to: {reports_dir}")
    return True
