# /GenoSynax_test/main.py

import argparse
import os
import sys
from argparse import Namespace

# --- Package Setup ---
project_root = os.path.dirname(os.path.abspath(__file__))
if project_root not in sys.path:
    sys.path.insert(0, project_root)

# --- Module Entrypoint Imports ---
try:
    from genosyntax.define_architecture import run as run_da
    from genosyntax.slot_analysis import run as run_sa
    from genosyntax.statistical_test import run as run_st
    from genosyntax.mechanism_analysis import run as run_ma
    from genosyntax import summarizer as run_summary
except ImportError as e:
    print(f"\\nFATAL: Could not import a necessary analysis module: {e}")
    sys.exit(1)

def run_pipeline(args):
    """Orchestrates the sequential execution of the GenoSyntax pipeline."""
    print("====== GenoSyntax v1.0: Initiating Integrated Pipeline Run ======\\n")

    main_output_dir = args.output_dir
    dir_step1 = os.path.join(main_output_dir, "1_define_architecture")
    dir_step2 = os.path.join(main_output_dir, "2_slot_analysis")
    dir_step3 = os.path.join(main_output_dir, "3_statistical_test")
    dir_step4 = os.path.join(main_output_dir, "4_mechanism_analysis")

    os.makedirs(dir_step1, exist_ok=True); os.makedirs(dir_step2, exist_ok=True)
    os.makedirs(dir_step3, exist_ok=True); os.makedirs(dir_step4, exist_ok=True)
    
    print("\n>>> EXECUTING STEP 1: define_architecture...")
    args_step1 = Namespace(**vars(args)); args_step1.output_dir = dir_step1
    if not run_da(args_step1):
        print("[ERROR] Pipeline halted due to failure in Step 1."); return
    if args.until_step == 1:
        print("\\n====== Pipeline finished after Step 1 as requested. ======"); return

    print("\n>>> EXECUTING STEP 2: analyze_slots...")
    args_step2 = Namespace(**vars(args))
    args_step2.architecture_file = os.path.join(dir_step1, "architecture_definition.json")
    args_step2.intermediate_dir_prev = os.path.join(dir_step1, "intermediate")
    args_step2.output_dir = dir_step2
    if not run_sa(args_step2):
        print("[ERROR] Pipeline halted due to failure in Step 2."); return
    if args.until_step == 2:
        print("\\n====== Pipeline finished after Step 2 as requested. ======"); return

    print("\n>>> EXECUTING STEP 3: run_stats...")
    args_step3 = Namespace(**vars(args))
    args_step3.slot_analysis_dir = dir_step2
    args_step3.output_dir = dir_step3
    if not run_st(args_step3):
        print("[ERROR] Pipeline halted due to failure in Step 3."); return

    print("\n>>> EXECUTING STEP 3.5: Generating Comprehensive Summary...")
    if not args.dry_run:
        if not run_summary.run(main_output_dir):
            print("[WARNING] Comprehensive summary generation failed, but pipeline will continue.")

    should_run_step_4 = (args.until_step == 4) or (args.mechanisms is not None)
    if not should_run_step_4:
        print("\\n====== Pipeline finished after Step 3 as requested. ======"); return

    print("\n>>> EXECUTING STEP 4: explore_mechanisms...")
    mechanisms_to_run = args.mechanisms if args.mechanisms else ['motif', 'evo_dynamics', 'adjacency_analysis', 'extract_layouts']
    
    args_step4 = Namespace(**vars(args))
    args_step4.analyses = mechanisms_to_run
    args_step4.def_arch_dir = dir_step1
    args_step4.slot_analysis_dir = dir_step2
    args_step4.stats_dir = dir_step3
    args_step4.output_dir = dir_step4
    
    if not run_ma(args_step4):
        print("[ERROR] Pipeline halted due to failure in Step 4."); return
        
    print("\\n====== Pipeline finished after Step 4. ======")


def main():
    """Main function to configure and parse command-line arguments."""
    parser = argparse.ArgumentParser(
        prog="genosyntax",
        description="GenoSyntax v1.0: A Framework for Genomic Architecture and Plasticity Analysis.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('-v', '--version', action='version', version='%(prog)s 1.0')

    subparsers = parser.add_subparsers(dest='command', required=True, help="Available commands")

    # --- Parent parser for shared STEP 1 arguments ---
    da_args = argparse.ArgumentParser(add_help=False)
    da_args.add_argument('--panaroo_file', required=True, help='Path to Panaroo gene_presence_absence_roary.csv file.')
    da_args.add_argument('--pharokka_dir', required=True, help='Path to Pharokka output directory (contains .gff, .faa, etc.).')
    da_args.add_argument('--fasta_dir', required=True, help='Path to directory with raw genome FASTA files (for QC).')
    da_args.add_argument('--faa_dir', help='Path to directory with protein FASTA (.faa) files.')
    da_args.add_argument('--output_dir', required=True, help='Path to the output directory for this module/pipeline.')
    da_args.add_argument('--max_contigs', type=int, default=1, help='(Step 1) Max contigs for a genome to pass QC (default: 1).')
    da_args.add_argument('--orientation_consistency_threshold', type=float, default=0.9, help='(Step 1) Threshold for consistent genome orientation (default: 0.9).')
    da_args.add_argument('--tolerance', type=int, default=2, help='(Step 1) Max accessory genes allowed within an anchor (default: 2).')
    da_args.add_argument('--coverage_threshold', type=float, default=0.9, help='(Step 1) Min fraction of genomes an anchor must be in (default: 0.9).')
    da_args.add_argument('--min_length', type=int, default=2, help='(Step 1) Min number of core genes in an anchor (default: 2).')

    # ============================ START OF S1 MODIFICATION ============================
    # Add new arguments for core gene stratification thresholds to the shared parent parser.
    da_args.add_argument('--strict_core_threshold', type=float, default=1.0, help="(Step 1) Prevalence threshold for a gene to be 'strict core' (default: 1.0).")
    da_args.add_argument('--soft_core_threshold', type=float, default=0.98, help="(Step 1) Prevalence threshold for a gene to be 'soft core', also defines accessory boundary (default: 0.98).")
    # ============================= END OF S1 MODIFICATION =============================

    # --- Sub-parser for the 'run' command ---
    p_run = subparsers.add_parser('run', help='(Recommended) Run the entire GenoSyntax pipeline sequentially.', parents=[da_args])
    p_run.add_argument('--spi_resolution', choices=['low', 'high'], default='low', help="(Step 2) Resolution for SPI calculation (default: 'low').")
    p_run.add_argument('--permutations', type=int, default=10000, help="(Step 3) Number of permutations for statistical tests (default: 10000).")
    p_run.add_argument('--spi_enrichment_threshold', type=float, default=1.0, help="(Step 4, evo_dynamics) SPI value above which a module is 'high-SPI' (default: 1.0).")
    p_run.add_argument('--until-step', type=int, default=3, choices=[1, 2, 3, 4], help='Stop the pipeline after this step. Default is 3.')
    p_run.add_argument('--mechanisms', nargs='+', choices=['motif', 'evo_dynamics', 'adjacency_analysis', 'extract_layouts', 'extract_cds_layout'], help='(Optional) Select one or more mechanism analyses to run as Step 4.')
    p_run.add_argument('--dry_run', action='store_true', help='Check inputs and logic for all steps without running heavy computations.')

    # ============================ START OF S1 MODIFICATION ============================
    # Add new arguments for Accessory Module discovery to the 'run' command parser.
    p_run.add_argument('--phi_threshold', type=float, default=0.8, help="(Step 2) Phi coefficient threshold for gene co-occurrence (default: 0.8).")
    p_run.add_argument('--p_value_threshold', type=float, default=0.01, help="(Step 2) P-value threshold for gene co-occurrence (default: 0.01).")
    p_run.add_argument('--adjacency_distance', type=int, default=250, help="(Step 2) Max distance (bp) between adjacent genes in a module (default: 250).")
    p_run.add_argument('--support_threshold', type=float, default=0.8, help="(Step 2) Min fraction of genomes where a module must be physically adjacent (default: 0.8).")
    # ============================= END OF S1 MODIFICATION =============================

    # ==============================================================================
    # --------------- DEFINING SUBPARSERS FOR INDIVIDUAL MODULES ------------------
    # ==============================================================================
    p_da = subparsers.add_parser('define_architecture', help='(Step 1) Identify and characterize the core genomic architecture.', parents=[da_args])
    p_da.add_argument('--dry_run', action='store_true', help='Check inputs and logic without running heavy computations.')

    p_sa = subparsers.add_parser('analyze_slots', help='(Step 2) Partition slot space and analyze accessory module dynamics.')
    p_sa.add_argument('--architecture_file', required=True, help='Path to the architecture_definition.json from the define_architecture step.')
    p_sa.add_argument('--intermediate_dir_prev', required=True, help="Path to the 'intermediate' output directory from the define_architecture step.")
    p_sa.add_argument('--output_dir', required=True, help='Path to the output directory for this module.')
    p_sa.add_argument('--spi_resolution', choices=['low', 'high'], default='low', help="Resolution for SPI calculation (default: 'low').")
    p_sa.add_argument('--dry_run', action='store_true', help='Check inputs and logic without running heavy computations.')

    # ============================ START OF S1 MODIFICATION ============================
    # Add the same Accessory Module arguments to the standalone 'analyze_slots' command.
    p_sa.add_argument('--phi_threshold', type=float, default=0.8, help="(Step 2) Phi coefficient threshold for gene co-occurrence (default: 0.8).")
    p_sa.add_argument('--p_value_threshold', type=float, default=0.01, help="(Step 2) P-value threshold for gene co-occurrence (default: 0.01).")
    p_sa.add_argument('--adjacency_distance', type=int, default=250, help="(Step 2) Max distance (bp) between adjacent genes in a module (default: 250).")
    p_sa.add_argument('--support_threshold', type=float, default=0.8, help="(Step 2) Min fraction of genomes where a module must be physically adjacent (default: 0.8).")
    # ============================= END OF S1 MODIFICATION =============================

    p_st = subparsers.add_parser('run_stats', help='(Step 3) Run statistical tests to identify hotspots/coldspots.')
    p_st.add_argument('--slot_analysis_dir', required=True, help="Path to the output directory from the analyze_slots step.")
    p_st.add_argument('--output_dir', required=True, help='Path to the output directory for this module.')
    p_st.add_argument('--permutations', type=int, default=10000, help="Number of permutations to run (default: 10000).")
    p_st.add_argument('--dry_run', action='store_true', help='Check inputs and logic without running heavy computations.')

    p_ma = subparsers.add_parser('explore_mechanisms', help='(Step 4) Prepare data and run downstream mechanism exploration.')
    p_ma.add_argument('--analyses', required=True, nargs='+', choices=['motif', 'evo_dynamics', 'adjacency_analysis', 'extract_layouts', 'extract_cds_layout'], help='Select one or more analyses to run.')
    p_ma.add_argument('--def_arch_dir', required=True, help="Path to the output directory from the define_architecture step.")
    p_ma.add_argument('--slot_analysis_dir', required=True, help="Path to the output directory from the analyze_slots step.")
    p_ma.add_argument('--stats_dir', required=True, help="Path to the output directory from the run_stats step.")
    p_ma.add_argument('--pharokka_dir', required=True, help='(Required for motif/evo_dynamics) Path to the Pharokka output directory.')
    p_ma.add_argument('--output_dir', required=True, help='Path to the output directory for this module.')
    p_ma.add_argument('--spi_enrichment_threshold', type=float, default=1.0, help="(For evo_dynamics) SPI value above which a module is 'high-SPI' (default: 1.0).")
    p_ma.add_argument('--dry_run', action='store_true', help='Check inputs and logic without running heavy computations.')

    # --- Argument Parsing and Dispatch ---\
    args = parser.parse_args()

    # Set faa_dir from pharokka_dir for compatibility
    if hasattr(args, 'pharokka_dir'):
        args.faa_dir = args.pharokka_dir
        
    if args.command == 'run':
        run_pipeline(args)
    elif args.command == 'define_architecture':
        run_da(args)
    elif args.command == 'analyze_slots':
        run_sa(args)
    elif args.command == 'run_stats':
        run_st(args)
    elif args.command == 'explore_mechanisms':
        run_ma(args)

if __name__ == '__main__':
    main()
