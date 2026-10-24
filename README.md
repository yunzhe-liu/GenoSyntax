# GenoSyntax

GenoSyntax compares conserved gene blocks and accessory gene locations in related phage genomes. It reads a Panaroo gene presence/absence matrix and Pharokka annotations, identifies conserved core-gene blocks (anchors), and describes accessory genes relative to those blocks (slots).

## Workflow

1. Classify gene families by prevalence and filter genomes by contig count and annotation availability.
2. Identify candidate anchors in the first genome passing QC and verify their order across the cohort.
3. Record nonempty slots between anchors and between core genes within anchors. Terminal slots depend on the cohort's start-gene consistency check.
4. Group accessory genes by co-occurrence and physical adjacency; retain remaining genes as single-gene modules.
5. Count the slot types occupied by each module and test accessory-gene enrichment across slot types.
6. Optionally export boundary sequences, protein sequence differences, annotation tables, and genome coordinates.

## Installation

Python 3.9 or later and MAFFT are required for the full analysis.

```bash
git clone https://github.com/yunzhe-liu/GenoSyntax.git
cd GenoSyntax
conda create -n genosyntax python=3.9
conda activate genosyntax
pip install -r requirements.txt
conda install -c bioconda mafft
```

## Inputs

Genome names and locus IDs must match across the matrix, annotations, and sequences.

| Argument | Expected input |
|---|---|
| `--panaroo_file` | Panaroo `gene_presence_absence_roary.csv`; the first 14 columns are metadata, followed by genome columns |
| `--pharokka_dir` | `<directory>/<genome>/<genome>.gff`, using CDS `ID` attributes that match matrix locus IDs |
| `--fasta_dir` | `<directory>/<genome>.fasta`, `.fa`, or `.fna`; coordinates must match the GFF |
| `--faa_dir` | Optional alternate root containing `<genome>/phanotate.faa`; defaults to `--pharokka_dir` |

The default filter accepts single-contig genomes. Supply consistently oriented genomes with comparable starting positions. The start-gene check reports consistency; sequence rotation and reverse complementation are external preparation steps.

## Usage

Run architecture definition, slot analysis, enrichment tests, and the summary:

```bash
python main.py run \
    --panaroo_file /path/to/gene_presence_absence_roary.csv \
    --pharokka_dir /path/to/pharokka_outputs \
    --fasta_dir /path/to/genome_fastas \
    --output_dir ./analysis \
    --seed 0
```

Add optional analyses with:

```bash
python main.py run \
    --panaroo_file /path/to/gene_presence_absence_roary.csv \
    --pharokka_dir /path/to/pharokka_outputs \
    --fasta_dir /path/to/genome_fastas \
    --output_dir ./analysis \
    --mechanisms motif evo_dynamics adjacency_analysis extract_layouts extract_cds_layout \
    --until-step 4 \
    --seed 0
```

`motif` exports boundary FASTA files for separate motif analysis. `evo_dynamics` aligns proteins with MAFFT and reports sequence differences. These optional analyses read `phanotate.faa` and nucleotide FASTA files from the Pharokka directory; nucleotide sequences must match the annotation coordinates.

`--until-step` stops after the selected stage. `--dry_run` checks required input paths and executable availability without running the analysis; file contents and biological consistency are checked during execution.

Individual commands are `define_architecture`, `analyze_slots`, `run_stats`, and `explore_mechanisms`. See [MANUAL.md](MANUAL.md) for parameters and output definitions, or run `python main.py <command> --help`.

Generate and run a small synthetic cohort:

```bash
python examples/make_demo.py /tmp/genosyntax_demo
python main.py run \
    --panaroo_file /tmp/genosyntax_demo/gene_presence_absence_roary.csv \
    --pharokka_dir /tmp/genosyntax_demo/pharokka \
    --fasta_dir /tmp/genosyntax_demo/fastas \
    --output_dir /tmp/genosyntax_demo/result \
    --permutations 100 --seed 0
```

The synthetic sequences exercise file handling and the workflow. They carry no biological interpretation.

## Outputs

| Location | Contents |
|---|---|
| `1_define_architecture/architecture_definition.json` | QC, gene categories, anchor definitions, and anchor-order groups |
| `1_define_architecture/intermediate/` | Gene order, coordinates, annotations, and locus mappings |
| `2_slot_analysis/intermediate/` | Slot contents and accessory module definitions |
| `2_slot_analysis/reports/module_dynamics_report.csv` | Module frequency within recorded slots, SPI, and observed slot labels |
| `3_statistical_test/reports/` | Enrichment and depletion tests for all accessory genes and for multi-gene module members |
| `comprehensive_summary_report.json` | Composite scores and selected results |
| `4_mechanism_analysis/` | Requested downstream exports; files depend on available data and selected analyses |

When no anchors satisfy the criteria, the pipeline saves the architecture report and stops. This outcome depends on the input cohort, reference genome, and thresholds.

## Interpreting results

- **SPI** counts distinct observed slot labels. Anchor recognition and slot resolution affect this count; it describes positional diversity.
- **Multi-gene modules** are co-occurring, physically adjacent candidates. Functional descriptions come from supplied annotations.
- **Slot enrichment** compares accessory-gene member counts with random assignments across pooled, observed slot gene positions. The null model pools genomes and treats cassette members individually. The primary reports contain unadjusted p-values.
- **PED** is one minus average protein alignment identity, averaged across a module's gene families.
- **GRI and the plasticity profile** are heuristic summaries with fixed weights. Cohort composition, sample size, thresholds, and slot definitions affect their values. Use the component measurements when comparing results.

The implementation uses linear gene order and reference-derived anchors. Circular boundaries, inversions, paralogs, and fragmented assemblies require inspection. The repository contains research analysis code developed for a phage cohort.

## License

[MIT](LICENSE).
