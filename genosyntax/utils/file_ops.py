# /GenoSynax_test/genosyntax/utils/file_ops.py

import json
import os
import subprocess
from itertools import combinations
import numpy as np
from scipy.stats import spearmanr
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def load_json(file_path):
    """Safely loads a JSON file."""
    try:
        with open(file_path, 'r') as f:
            return json.load(f)
    except FileNotFoundError:
        print(f"[ERROR] JSON file not found: {file_path}")
        return None

def save_json(data, file_path):
    """Saves data to a JSON file with indentation."""
    try:
        with open(file_path, 'w') as f:
            json.dump(data, f, indent=4)
        return True
    except IOError as e:
        print(f"[ERROR] Could not write to JSON file {file_path}: {e}")
        return False

def calculate_asi_for_gene(sequences, temp_dir):
    """
    Calculates the average pairwise identity for a single gene's sequences
    using an external alignment tool like MAFFT.
    """
    if len(sequences) < 2:
        return 1.0

    pid = os.getpid()
    temp_fasta_path = os.path.join(temp_dir, f"temp_asi_{pid}.fasta")
    
    records = [SeqRecord(Seq(seq), id=f"seq_{i}") for i, seq in enumerate(sequences)]
    SeqIO.write(records, temp_fasta_path, "fasta")

    aln_path = temp_fasta_path + ".aln"
    
    try:
        subprocess.run(
            ["mafft", "--auto", "--quiet", temp_fasta_path],
            stdout=open(aln_path, "w"),
            stderr=subprocess.DEVNULL,
            text=True,
            check=True
        )
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        if isinstance(e, FileNotFoundError):
            print("\n[FATAL ERROR] MAFFT command not found. Is it installed and in your PATH?")
            print("  To install: conda install -c bioconda mafft")
            exit(1)
        return None

    try:
        alignment = list(SeqIO.parse(aln_path, "fasta"))
    except Exception:
        alignment = []
        
    pids = []
    if alignment:
        for rec1, rec2 in combinations(alignment, 2):
            seq1_str, seq2_str = str(rec1.seq), str(rec2.seq)
            matches = sum(1 for a, b in zip(seq1_str, seq2_str) if a == b and a != '-')
            len1 = len(seq1_str.replace('-', ''))
            len2 = len(seq2_str.replace('-', ''))
            min_len = min(len1, len2)
            if min_len > 0:
                pids.append(matches / min_len)
    
    if os.path.exists(temp_fasta_path): os.remove(temp_fasta_path)
    if os.path.exists(aln_path): os.remove(aln_path)

    return np.mean(pids) if pids else 1.0

# --- MODIFIED FUNCTION ---
def calculate_gc_content(sequence):
    """
    Calculates GC content of a DNA sequence string.
    This is now the primary function for boundary signal analysis.
    """
    if not isinstance(sequence, str) or len(sequence) == 0:
        return 0.0
    g = sequence.upper().count('G')
    c = sequence.upper().count('C')
    return (g + c) / len(sequence)
