# /GenoSynax_test/genosyntax/summarizer.py

import json
import os
import pandas as pd
import numpy as np
from scipy.stats import entropy

def _load_data(path, file_type='json'):
    """Helper function to safely load JSON or CSV files."""
    if not os.path.exists(path):
        print(f"  [Info] Optional input file not found, will proceed without it: {os.path.basename(path)}")
        return None
    try:
        if file_type == 'json':
            with open(path, 'r') as f:
                return json.load(f)
        elif file_type == 'csv':
            return pd.read_csv(path)
    except Exception as e:
        print(f"  [ERROR] Could not load or parse file {path}: {e}")
        return None

def calculate_gri(arch_def, weights={'c': 0.6, 's': 0.3, 'm': 0.1}):
    """
    Calculates the Genomic Rigidity Index (GRI) based on the architecture definition.
    """
    if not arch_def or not arch_def.get("anchor_definitions"):
        return {
            "GRI_score": 0.0,
            "details": {
                "backbone_coverage": 0.0,
                "backbone_stability": 0.0,
                "macro_conservation": 0.0,
            },
            "status": "No anchor modules were found, rigidity is considered minimal."
        }

    # 1. Backbone Coverage (C_backbone)
    anchor_core_genes = {gene for anchor in arch_def['anchor_definitions'] for gene in anchor['skeleton']}
    total_core_genes = set(arch_def.get('metadata', {}).get('pangenome_stratification', {}).get('strict_core', []))
    c_backbone = len(anchor_core_genes) / len(total_core_genes) if total_core_genes else 0.0

    # 2. Backbone Stability (S_backbone)
    anchors = arch_def['anchor_definitions']
    avg_coverage = np.mean([a.get('coverage_fraction', 0) for a in anchors])
    avg_asi = np.mean([a.get('asi', 0) for a in anchors])
    avg_aoi = np.mean([a.get('aoi', 0) for a in anchors])
    s_backbone = (avg_coverage + avg_asi + avg_aoi) / 3 if anchors else 0.0

    # 3. Macro-structural Conservation (M_conservation)
    macro_arch = arch_def.get('macro_architectures', {})
    E = macro_arch.get('global_permutation_entropy', 0)
    N = macro_arch.get('total_architectures_found', 1)
    if N <= 1:
        m_conservation = 1.0
    else:
        E_max = np.log2(N)
        m_conservation = 1 - (E / E_max) if E_max > 0 else 1.0
        
    gri_score = (weights['c'] * c_backbone +
                 weights['s'] * s_backbone +
                 weights['m'] * m_conservation)

    return {
        "GRI_score": round(gri_score, 4),
        "details": {
            "backbone_coverage": round(c_backbone, 4),
            "backbone_stability": round(s_backbone, 4),
            "macro_conservation": round(m_conservation, 4),
        },
        "status": "Calculated successfully."
    }

def _calculate_hotspot_activity(stats_df):
    """Helper to calculate total activity score from a permutation test report."""
    if stats_df is None or stats_df.empty:
        return 0.0
    
    hotspots = stats_df[stats_df['p_value_enrichment'] < 0.05].copy()
    if hotspots.empty:
        return 0.0
        
    hotspots.loc[:, 'p_value_enrichment'] = hotspots['p_value_enrichment'].replace(0, 1e-300)
    hotspots.loc[:, 'fold_enrichment'] = pd.to_numeric(hotspots['fold_enrichment'], errors='coerce').fillna(0)
    
    hotspots['score'] = hotspots['fold_enrichment'] * (-np.log10(hotspots['p_value_enrichment']))
    return hotspots['score'].sum()

def calculate_gps(module_dynamics_report, stats_report_all, stats_report_multi, weights={'m': 0.4, 'h': 0.3, 'c': 0.3}):
    """
    Calculates the Genomic Plasticity Score (GPS).
    """
    # 1. Genetic Mobility (M_genetic)
    if module_dynamics_report is not None and not module_dynamics_report.empty:
        total_freq = module_dynamics_report['frequency'].sum()
        m_genetic = (module_dynamics_report['frequency'] * module_dynamics_report['spi']).sum() / total_freq if total_freq > 0 else 1.0
    else:
        m_genetic = 1.0

    # 2. Hotspot Activity (H_activity)
    h_activity = _calculate_hotspot_activity(stats_report_all)

    # 3. Cassette-Driven Plasticity (C_cassette)
    c_cassette = _calculate_hotspot_activity(stats_report_multi)

    gps_score = (weights['m'] * m_genetic +
                 weights['h'] * h_activity +
                 weights['c'] * c_cassette)
                 
    return {
        "GPS_score": round(gps_score, 4),
        "details": {
            "genetic_mobility": round(m_genetic, 4),
            "hotspot_activity_all_genes": round(h_activity, 4),
            "hotspot_activity_cassettes": round(c_cassette, 4),
        }
    }

def extract_key_findings(arch_def, module_dynamics_report, stats_report_all):
    """Extracts human-readable highlights from the analysis."""
    findings = {}
    
    # Architecture findings
    if arch_def and arch_def.get("anchor_definitions"):
        major_arch = arch_def.get('macro_architectures', {}).get('major_architecture', {})
        findings['dominant_architecture'] = {
            "permutation": " -> ".join(major_arch.get('permutation', [])),
            "frequency": major_arch.get('frequency', 0)
        }
        findings['anchor_module_count'] = len(arch_def["anchor_definitions"])
    else:
        findings['dominant_architecture'] = "No stable architecture found."
        findings['anchor_module_count'] = 0

    # Plasticity findings
    if module_dynamics_report is not None and not module_dynamics_report.empty:
        top_plastic = module_dynamics_report.sort_values(by='spi', ascending=False).head(5)
        findings['top_5_most_plastic_modules'] = top_plastic[['module_id', 'spi', 'frequency', 'consensus_function']].to_dict('records')
    
    if stats_report_all is not None and not stats_report_all.empty:
        top_hotspots = stats_report_all[stats_report_all['p_value_enrichment'] < 0.05].sort_values(by='p_value_enrichment').head(3)
        findings['top_3_insertion_hotspots'] = top_hotspots[['slot_type', 'p_value_enrichment', 'fold_enrichment']].to_dict('records')
        
    return findings


def run(main_output_dir):
    """
    Entrypoint function to generate the comprehensive summary report.
    """
    print("\n--- Running Step 3.5: Generating Comprehensive Summary Report ---")
    
    # Define paths
    arch_def_path = os.path.join(main_output_dir, "1_define_architecture", "architecture_definition.json")
    module_report_path = os.path.join(main_output_dir, "2_slot_analysis", "reports", "module_dynamics_report.csv")
    stats_all_path = os.path.join(main_output_dir, "3_statistical_test", "reports", "permutation_test_slot_significance.csv")
    stats_multi_path = os.path.join(main_output_dir, "3_statistical_test", "reports", "permutation_test_multigene_module_significance.csv")
    
    # Load data
    arch_def = _load_data(arch_def_path, 'json')
    module_dynamics_report = _load_data(module_report_path, 'csv')
    stats_report_all = _load_data(stats_all_path, 'csv')
    stats_report_multi = _load_data(stats_multi_path, 'csv')

    if not arch_def:
        print("  [FATAL] `architecture_definition.json` is essential and was not found. Cannot generate summary.")
        return False

    # Perform calculations
    gri_results = calculate_gri(arch_def)
    gps_results = calculate_gps(module_dynamics_report, stats_report_all, stats_report_multi)
    key_findings = extract_key_findings(arch_def, module_dynamics_report, stats_report_all)

    # Assemble final report
    summary_report = {
        "metadata": {
            "report_generator": "GenoSyntax v1.0 Summarizer",
            "qc_summary": arch_def.get("metadata", {}).get("qc_summary", {}),
            "analysis_parameters": arch_def.get("metadata", {}).get("analysis_parameters", {}),
        },
        "genomic_rigidity_index": gri_results,
        "genomic_plasticity_profile": {
            "genetic_mobility_score": gps_results["details"]["genetic_mobility"],
            "general_insertion_activity_score": gps_results["details"]["hotspot_activity_all_genes"],
            "cassette_driven_plasticity_score": gps_results["details"]["hotspot_activity_cassettes"],
            "legacy_gps_combined_score": gps_results["GPS_score"]
        },
        "key_findings_summary": key_findings
    }

    # Save the report
    report_path = os.path.join(main_output_dir, "comprehensive_summary_report.json")
    try:
        with open(report_path, 'w') as f:
            json.dump(summary_report, f, indent=4)
        print(f"[*] Comprehensive summary report saved successfully to: {report_path}")
    except IOError as e:
        print(f"  [FATAL] Could not write summary report to {report_path}: {e}")
        return False
        
    return True
