# /GenoSynax_test/genosyntax/define_architecture/quality_control.py

import os
from collections import Counter
from Bio import SeqIO

def analyze_fasta_for_qc(fasta_path):
    """
    A lightweight and robust function to analyze a FASTA file for QC purposes.
    VERSION 2 (ROBUST): Instead of relying on Bio.SeqIO.parse(), this version
    reads the file content directly and counts '>' occurrences. This is more
    resilient to non-standard line endings and other formatting issues that
    can trick the parser into seeing multiple contigs where there is only one.
    """
    if not os.path.exists(fasta_path):
        return {'num_contigs': -1, 'error': 'FASTA file not found'}
    
    try:
        contig_count = 0
        with open(fasta_path, 'r', errors='ignore') as f:
            # Read the file line by line to count headers
            for line in f:
                if line.startswith('>'):
                    contig_count += 1
        
        # If the file is empty or contains no headers, it's an issue.
        if contig_count == 0:
            return {'num_contigs': 0, 'error': 'No FASTA headers found in file'}
            
        return {'num_contigs': contig_count}

    except Exception as e:
        # This will catch file reading errors, encoding errors, etc.
        return {'num_contigs': -1, 'error': f'File reading failed: {e}'}


def run_qc(genome_ids, fasta_dir, gene_order_db, gene_coord_db, max_contigs=1, orientation_threshold=0.9):
    """
    Performs QC checks.
    This function now uses the robust analyze_fasta_for_qc v2.
    """
    print("\n--- Task 2: Running Quality Control ---")
    
    total_genomes_input = genome_ids
    excluded_genomes = {}
    
    # --- QC Check 1: Contig Count from SOURCE FASTA files ---
    print(f"QC 1: Checking contig counts from FASTA files in '{fasta_dir}' (threshold: max_contigs <= {max_contigs})...")
    genomes_after_contig_qc = []
    fasta_contig_counts = {}
    
    for genome_id in total_genomes_input:
        fasta_path = None
        # Robustly find the FASTA file with common extensions
        for ext in ['.fasta', '.fa', '.fna']:
            potential_path = os.path.join(fasta_dir, f"{genome_id}{ext}")
            if os.path.exists(potential_path):
                fasta_path = potential_path
                break
        
        # This now calls our new, robust function
        qc_result = analyze_fasta_for_qc(fasta_path)
        num_contigs = qc_result.get('num_contigs')
        fasta_contig_counts[genome_id] = num_contigs
        
        if num_contigs in [-1, 0]: # Handle file not found or empty/no-header files
            reason = qc_result.get('error', 'Unknown file error')
            excluded_genomes[genome_id] = f"QC Error: {reason}"
            continue
            
        if not gene_coord_db.get(genome_id):
            excluded_genomes[genome_id] = "Data Loading Failed (GFF not found or empty)"
            continue
            
        if num_contigs > max_contigs:
            excluded_genomes[genome_id] = f"Fragmented Assembly ({num_contigs} contigs)"
        else:
            genomes_after_contig_qc.append(genome_id)

    print(f"  {len(genomes_after_contig_qc)}/{len(total_genomes_input)} genomes passed all filtering criteria.")

    # --- QC Check 2: Orientation Consistency ---
    print(f"QC 2: Checking orientation consistency (threshold: >{orientation_threshold*100}%)...")
    if not genomes_after_contig_qc:
        orientation_summary = {
            'is_consistently_oriented': False,
            'orientation_gene': None,
            'confidence': 0.0,
            'details': "No high-quality genomes available for analysis."
        }
    else:
        # Ensure we only check genomes that have a gene order available
        valid_genomes_for_orient = [gid for gid in genomes_after_contig_qc if gene_order_db.get(gid)]
        
        if not valid_genomes_for_orient:
             orientation_summary = {
                'is_consistently_oriented': False,
                'orientation_gene': None,
                'confidence': 0.0,
                'details': "High-quality genomes have empty gene orders."
            }
        else:
            start_genes = [gene_order_db[gid][0] for gid in valid_genomes_for_orient]
            start_gene_counts = Counter(start_genes)
            
            dominant_start_gene, count = start_gene_counts.most_common(1)[0]
            confidence = count / len(valid_genomes_for_orient)
            is_consistent = confidence >= orientation_threshold
            
            orientation_summary = {
                'is_consistently_oriented': is_consistent,
                'orientation_gene': dominant_start_gene,
                'confidence': round(confidence, 4),
                'details': f"Dominant start gene '{dominant_start_gene}' found in {count}/{len(valid_genomes_for_orient)} ({confidence:.2%}) of high-quality genomes."
            }
        print(f"  Orientation check complete. Consistent: {orientation_summary['is_consistently_oriented']}.")

    # --- Final Report Generation ---
    qc_report = {
        'parameters': { 
            'max_contigs': max_contigs, 
            'orientation_consistency_threshold': orientation_threshold 
        },
        'summary': { 
            'total_genomes_input': len(total_genomes_input), 
            'genomes_passed_qc': len(genomes_after_contig_qc), 
            'genomes_failed_qc': len(excluded_genomes) 
        },
        'failed_genomes': excluded_genomes,
        'fasta_contig_counts': fasta_contig_counts,
        'orientation_analysis': orientation_summary,
        'passed_genomes_list': sorted(genomes_after_contig_qc)
    }
    
    print("Quality control finished.")
    return genomes_after_contig_qc, qc_report
