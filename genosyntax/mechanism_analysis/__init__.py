import os
import pandas as pd
from collections import defaultdict

from genosyntax.utils.file_ops import load_json
from .motif_finder import extract_boundary_sequences
from .evolutionary_analyzer import run_evolutionary_analysis
from .adjacency_analyzer import run_adjacency_analysis
from .layout_extractor import run_layout_extraction
from .cds_extractor import run_cds_extraction

def run(args):
    """Entrypoint for 'explore_mechanisms' module."""
    if args.dry_run:
        print(f"--- Dry Run for 'explore_mechanisms' ---")
        return True
        
    output_dir = args.output_dir
    os.makedirs(output_dir, exist_ok=True)

    print("\n====== GenoSyntax: Running Module 'explore_mechanisms' ======")
    
    print("--- Loading prerequisite data from previous modules ---")
    def_arch_dir = args.def_arch_dir
    slot_analysis_dir = args.slot_analysis_dir
    
    architecture_def = load_json(os.path.join(def_arch_dir, 'architecture_definition.json'))
    gene_order_db = load_json(os.path.join(def_arch_dir, 'intermediate/gene_order_database.json'))
    gene_coord_db = load_json(os.path.join(def_arch_dir, 'intermediate/gene_coordinate_database.json'))
    
    if not all([architecture_def, gene_order_db, gene_coord_db]):
        print("[FATAL] Could not load essential data from 'define_architecture' step. Aborting.")
        return False

    # --- Sub-module Dispatcher ---
    
    if 'extract_layouts' in args.analyses:
        print("\n--- Running Task: Extracting All Genome Layouts ---")
        run_layout_extraction(output_dir, architecture_def, gene_order_db, gene_coord_db)
    
    if 'extract_cds_layout' in args.analyses:
        print("\n--- Running Task: Extracting All CDS Layouts ---")
        gene_annotation_map = load_json(os.path.join(def_arch_dir, 'intermediate/gene_annotation_map.json'))
        if not gene_annotation_map:
            print("[FATAL] Could not load gene_annotation_map.json, required for CDS layout extraction. Aborting task.")
        else:
            run_cds_extraction(output_dir, architecture_def, gene_coord_db, gene_annotation_map)

    if 'evo_dynamics' in args.analyses:
        print("\n--- Running Task: Evolutionary Dynamics Analysis ---")
        df_spi = pd.read_csv(os.path.join(slot_analysis_dir, 'reports/module_dynamics_report.csv'))
        accessory_modules = load_json(os.path.join(slot_analysis_dir, 'intermediate/accessory_module_definitions.json'))
        run_evolutionary_analysis(
            output_dir, df_spi, accessory_modules, gene_coord_db, 
            args.pharokka_dir, getattr(args, 'spi_enrichment_threshold', 1.0)
        )

    if 'adjacency_analysis' in args.analyses:
        print("\n--- Running Task: Adjacency Module Functional Analysis ---")
        accessory_modules = load_json(os.path.join(slot_analysis_dir, 'intermediate/accessory_module_definitions.json'))
        gene_annotation_map = load_json(os.path.join(def_arch_dir, 'intermediate/gene_annotation_map.json'))
        run_adjacency_analysis(output_dir, accessory_modules, gene_annotation_map)

    if 'motif' in args.analyses:
        print("\n--- Running Task: Motif Discovery Preparation (All Modes) ---")
        slot_database = load_json(os.path.join(slot_analysis_dir, 'intermediate/slot_database.json'))
        
        # Mode 1: ALL accessory modules
        print("  -> Mode 1/3: Checking for 'all module' hotspots...")
        stats_report_all_path = os.path.join(args.stats_dir, 'reports/permutation_test_slot_significance.csv')
        if os.path.exists(stats_report_all_path):
            stats_report_all = pd.read_csv(stats_report_all_path)
            motif_dir_all = os.path.join(output_dir, "motif_analysis_all_modules")
            os.makedirs(motif_dir_all, exist_ok=True)
            extract_boundary_sequences(
                slot_database, gene_order_db, gene_coord_db, args.pharokka_dir, motif_dir_all,
                statistical_report_df=stats_report_all
            )
        else:
            print("     [INFO] Statistics report for 'all modules' not found. Skipping.")
        
        print("  -> Mode 2/3: Checking for 'multi-gene module' hotspots...")
        stats_report_multi_path = os.path.join(args.stats_dir, 'reports/permutation_test_multigene_module_significance.csv')
        if os.path.exists(stats_report_multi_path):
            stats_report_multi = pd.read_csv(stats_report_multi_path)
            motif_dir_multi = os.path.join(output_dir, "motif_analysis_multigene_modules")
            os.makedirs(motif_dir_multi, exist_ok=True)
            extract_boundary_sequences(
                slot_database, gene_order_db, gene_coord_db, args.pharokka_dir, motif_dir_multi,
                statistical_report_df=stats_report_multi
            )
        else:
            print("     [INFO] Statistics report for 'multi-gene modules' not found. Skipping.")
        
        print("  -> Mode 3/3: Checking for slots enriched in 'high-SPI' modules...")
        high_spi_enrichment_path = os.path.join(output_dir, 'reports/high_spi_slot_enrichment.csv')
        if os.path.exists(high_spi_enrichment_path):
            high_spi_enrichment_df = pd.read_csv(high_spi_enrichment_path)
            significant_slots = high_spi_enrichment_df[high_spi_enrichment_df['q_value'] < 0.05]['slot_type'].tolist()
            if significant_slots:
                targeted_slot_instances = defaultdict(dict)
                for genome_id, slots in slot_database.items():
                    for slot_id, slot_info in slots.items():
                        if slot_info['slot_id_low_res'] in significant_slots:
                            targeted_slot_instances[genome_id][slot_id] = slot_info
                motif_dir_high_spi = os.path.join(output_dir, "motif_analysis_high_spi_enriched")
                os.makedirs(motif_dir_high_spi, exist_ok=True)
                extract_boundary_sequences(
                    slot_database, gene_order_db, gene_coord_db, args.pharokka_dir, motif_dir_high_spi,
                    target_slot_instances=targeted_slot_instances
                )
            else:
                print("     [INFO] No statistically significant slots for high-SPI modules found. Skipping.")
        else:
            print("     [INFO] High-SPI enrichment report not found. Skipping. (Hint: Run 'evo_dynamics' analysis).")
        
    print("\n====== 'explore_mechanisms' Module Finished Successfully! ======")
    return True
