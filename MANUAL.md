# GenoSyntax usage and calculations

## Parameters

Use `python main.py <command> --help` for the full argument list.

| Stage | Parameter | Default | Meaning |
|---|---|---|---|
| Architecture | `--strict_core_threshold` | 1.0 | Minimum prevalence for strict-core genes, calculated before QC filtering |
| Architecture | `--soft_core_threshold` | 0.98 | Minimum prevalence for soft-core genes; lower-prevalence families are accessory |
| Architecture | `--max_contigs` | 1 | Maximum FASTA record count |
| Architecture | `--orientation_consistency_threshold` | 0.9 | Fraction sharing the same first mapped gene for terminal-slot inclusion |
| Architecture | `--tolerance` | 2 | Maximum consecutive non-core genes between successive skeleton genes |
| Architecture | `--coverage_threshold` | 0.9 | Minimum fraction of QC-passing genomes supporting an anchor |
| Architecture | `--min_length` | 2 | Minimum core-gene count in an anchor |
| Slots | `--phi_threshold` | 0.8 | Minimum gene-pair presence/absence correlation |
| Slots | `--p_value_threshold` | 0.01 | Unadjusted one-sided Fisher p-value cutoff for gene pairs |
| Slots | `--adjacency_distance` | 250 | Maximum separation in bp between module members; overlaps are rejected |
| Slots | `--support_threshold` | 0.8 | Minimum adjacency fraction among genomes containing all candidate members |
| Slots | `--spi_resolution` | low | Merge intra-anchor positions by anchor; `high` retains flanking-gene labels |
| Statistics | `--permutations` | 10000 | Number of random assignments |
| Downstream | `--spi_enrichment_threshold` | 1.0 | SPI values above this threshold enter the high-SPI group |

`--seed` defaults to 0 for the integrated pipeline, slot analysis, and statistical tests. It seeds Python and NumPy randomness and Louvain clustering. Reproduction also requires the same inputs, parameters, and dependency versions. Python dependencies are currently unpinned.

## Separate stages

```bash
python main.py define_architecture \
    --panaroo_file /path/to/gene_presence_absence_roary.csv \
    --pharokka_dir /path/to/pharokka_outputs \
    --fasta_dir /path/to/genome_fastas \
    --output_dir ./analysis/1_define_architecture

python main.py analyze_slots \
    --architecture_file ./analysis/1_define_architecture/architecture_definition.json \
    --intermediate_dir_prev ./analysis/1_define_architecture/intermediate \
    --output_dir ./analysis/2_slot_analysis --seed 0

python main.py run_stats \
    --slot_analysis_dir ./analysis/2_slot_analysis \
    --output_dir ./analysis/3_statistical_test --seed 0

python main.py explore_mechanisms \
    --analyses motif evo_dynamics adjacency_analysis extract_layouts extract_cds_layout \
    --def_arch_dir ./analysis/1_define_architecture \
    --slot_analysis_dir ./analysis/2_slot_analysis \
    --stats_dir ./analysis/3_statistical_test \
    --pharokka_dir /path/to/pharokka_outputs \
    --output_dir ./analysis/4_mechanism_analysis
```

The integrated `run` command additionally generates the summary after statistical testing. Use a fresh output directory when changing inputs or parameters to avoid retaining files from an earlier run.

## Calculations

### Anchors and architecture

Candidate skeletons come from the first QC-passing genome in matrix-column order. Candidates are verified in the forward gene order across QC-passing genomes. Failed candidates are discarded as whole blocks.

ASI averages protein identity across anchor gene families after MAFFT alignment. Its identity denominator is the shorter ungapped sequence length. Families with fewer than two supplied sequences receive 1.0.

AOI averages Spearman correlations between gene rank columns across genomes. It measures covariation of absolute ranks. Constant rank columns can yield undefined values.

Anchor boundary GC values use up to 200 bp outside each anchor, standardized against the distribution of whole-genome GC fractions.

Architecture grouping merges observed anchor orders that are ordered subsequences of longer observed orders. Reported group frequencies therefore include genomes with incomplete anchor sets.

### Slots and modules

Recorded slots contain mapped genes. Empty intervals are omitted. Upstream and downstream slots are recorded when the first-gene consistency threshold is met. Soft-core and strict-core genes outside anchors may also occur in slots.

Accessory gene pairs pass a Phi threshold and one-sided Fisher test. Louvain communities become module candidates, followed by physical adjacency checks. Pairwise p-values are unadjusted. Each remaining accessory gene becomes a single-gene module.

A module is observed in a slot when all its members occur in that slot. Frequency is the fraction of QC-passing genomes with such an observation. SPI is the number of distinct observed slot labels. The existing behavioral labels use SPI and a frequency cutoff of 0.5.

### Enrichment tests

The primary tests pool all recorded slot gene positions across genomes. Each permutation samples as many positions as the observed accessory-member count, then counts sampled positions by low-resolution slot type. The second test restricts the counted genes to multi-gene module members.

One-sided empirical p-values use `(extreme permutations + 1) / (permutations + 1)`. Fold enrichment divides the observed count by the simulated mean. Reports contain raw p-values; downstream boundary selection uses a raw 0.05 cutoff.

The separate high-SPI analysis compares module-to-slot label observations using Fisher tests and Benjamini–Hochberg correction. A module can contribute observations to several slots.

### Summary scores

GRI combines anchor coverage of strict-core genes (weight 0.6), mean anchor coverage/ASI/AOI (0.3), and normalized architecture-group entropy conservation (0.1).

The plasticity profile reports frequency-weighted SPI and two enrichment summaries. Each enrichment summary sums `fold_enrichment × -log10(p)` over raw-p significant slot types. The legacy combined GPS uses weights 0.4, 0.3, and 0.3. These scores have no biological calibration in this repository.

## Optional exports

| Analysis | Output |
|---|---|
| `motif` | Boundary FASTA files grouped by enrichment class; high-SPI boundaries require a significant high-SPI enrichment report |
| `evo_dynamics` | `reports/spi_vs_ped_analysis.csv`; `reports/high_spi_slot_enrichment.csv` when testable |
| `adjacency_analysis` | Module annotations, annotation-coverage categories, and member tables for Cytoscape |
| `extract_layouts` | `reports/all_genome_layouts.csv`; uses the sibling `2_slot_analysis` directory |
| `extract_cds_layout` | `reports/all_genomes_cds_layout.csv`; CDS entries mapped to Panaroo families |

Motif extraction uses up to 200 bp of intergenic sequence at slot boundaries and retains sequences at least 20 bp long. Run motif discovery separately. Protein differences use MAFFT alignments with gap-containing columns included except double gaps. Annotation categories use keyword matching.

## Troubleshooting

- Missing MAFFT: activate the environment containing the executable and check `mafft --version`.
- Missing annotations: check genome directory names, `<genome>.gff`, and the correspondence between CDS `ID` attributes and matrix locus IDs.
- No anchors: inspect core-gene prevalence, input orientation, the reference genome, and anchor thresholds. The pipeline retains the architecture report and stops.
- Missing boundary exports: inspect whether selected slots have available nucleotide sequences and intergenic segments meeting the minimum length.

## Development notebooks

The root notebooks are retained development snapshots containing duplicate code and plotting snippets. The Python modules and `main.py` define the supported command-line workflow.

## Local checks

Run `python -m unittest discover -s tests -v` from the repository root. CLI checks use the standard library; analysis checks require NumPy, pandas, and tqdm. The synthetic example additionally requires the full dependencies and MAFFT.
