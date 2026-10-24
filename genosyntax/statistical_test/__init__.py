import os
import json

from genosyntax.utils.file_ops import load_json
from genosyntax.statistical_test.permutation_test import run_permutation_test
import pandas as pd

def run(args):
    """Entrypoint function for the 'run_stats' module."""
    if args.dry_run:
        print(f"--- Dry Run for 'run_stats' ---")
        print("Arguments received:", vars(args))
        if not os.path.exists(args.slot_analysis_dir):
            print(f"  [ERROR] Input path not found: {args.slot_analysis_dir}")
            return False
        print("  All input paths are valid.")
        print("Dry run complete.")
        return True
        
    output_dir = args.output_dir
    reports_dir = os.path.join(output_dir, "reports")
    os.makedirs(reports_dir, exist_ok=True)

    print("====== GenoSyntax: Running Module 'statistical_test' ======")
    
    # Load data from the 'slot_analysis' step
    slot_db_path = os.path.join(args.slot_analysis_dir, 'intermediate/slot_database.json')
    am_def_path = os.path.join(args.slot_analysis_dir, 'intermediate/accessory_module_definitions.json')
    
    print(f"[*] Loading slot database from: {slot_db_path}")
    slot_database = load_json(slot_db_path)
    
    print(f"[*] Loading accessory module definitions from: {am_def_path}")
    accessory_modules = load_json(am_def_path)
    
    if slot_database is None or accessory_modules is None:
        print("[FATAL] Could not load necessary input files. Please ensure the 'analyze_slots' module ran successfully.")
        return False

    results_df_all, results_df_multi = run_permutation_test(
        slot_database, accessory_modules, num_permutations=args.permutations
    )

    
    # Save the original, primary report for all accessory modules
    if not results_df_all.empty:
        report_path_all = os.path.join(reports_dir, "permutation_test_slot_significance.csv")
        results_df_all.to_csv(report_path_all, index=False)
        print(f"\n[SUCCESS] Report for ALL accessory modules saved to: {report_path_all}")
    else:
        print("\n[INFO] No results were generated for the 'all modules' analysis.")

    if not results_df_multi.empty:
        report_path_multi = os.path.join(reports_dir, "permutation_test_multigene_module_significance.csv")
        results_df_multi.to_csv(report_path_multi, index=False)
        print(f"[SUCCESS] Additional report for MULTI-GENE modules saved to: {report_path_multi}")
    else:
        print("[INFO] No results were generated for the 'multi-gene modules' analysis.")

    print(f"\n====== 'run_stats' Module Finished Successfully! ======")
    print(f"[*] Final reports saved to: {reports_dir}")

    return True
