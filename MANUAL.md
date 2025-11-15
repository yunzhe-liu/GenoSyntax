
# GenoSyntax v1.0 - User Manual

## 1. Introduction

### 1.1. Overview
GenoSyntax is a comprehensive bioinformatics framework for the analysis of genomic architecture and plasticity in related microbial genomes, particularly suited for phage and bacterial pangenomes. It deconstructs complex pangenomes into a stable, conserved **backbone** (defined by "Anchor Modules") and dynamic **slots** where mobile genetic elements reside.

The framework is designed as a sequential pipeline that automates the entire analytical process, from raw annotations to quantitative, high-level insights. It culminates in the generation of two novel metrics—the **Genomic Rigidity Index (GRI)** and the **Genomic Plasticity Score (GPS)**—which allow for direct, quantitative comparison of architectural dynamics between different datasets.

### 1.2. Conceptual Workflow
The analysis is performed in four main stages:
1.  **Define Architecture:** Identify the conserved genomic backbone.
2.  **Analyze Slots:** Characterize the plastic regions and the mobile elements within them.
3.  **Statistical Testing:** Determine the significance of mobile element insertion patterns.
4.  **Explore Mechanisms:** Prepare data to investigate the molecular and evolutionary drivers of the observed architecture.

---

## 2. Installation

GenoSyntax requires a Conda environment to manage its dependencies and external tools.

### 2.1. Prerequisites
-   **Conda:** [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or Anaconda must be installed.
-   **Git:** For cloning the repository.

### 2.2. Installation Steps

1.  **Clone the Repository:**
    ```bash
    git clone https://github.com/yunzhe-liu/GenoSynax.git
    cd GenoSyntax
    ```

2.  **Create and Activate Conda Environment:**
    A dedicated environment is highly recommended.
    ```bash
    conda create -n genosyntax python=3.9
    conda activate genosyntax
    ```

3.  **Install Dependencies:**
    Install the required Python packages and external bioinformatics tools.
    ```bash
    # Install Python packages from the requirements file
    pip install -r requirements.txt
    
    # Install essential external tools via Conda
    conda install -c bioconda mafft
    ```

---

## 3. Input Data Requirements

GenoSyntax requires four sources of data, provided as command-line arguments.

| Argument                | Description                                                                                                                                                             | Source Tool |
| ----------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ----------- |
| `--panaroo_file`        | The `gene_presence_absence_roary.csv` file. This provides the pangenome definition and gene family prevalence.                                                           | Panaroo     |
| `--pharokka_dir`        | The main output directory from Pharokka. GenoSyntax expects this directory to contain a subdirectory for each genome, which holds the `.gff`, `.faa`, and reoriented `.fna` files. | Pharokka    |
| `--fasta_dir`           | A directory containing the original, raw FASTA sequence files for each genome (one file per genome). This is used for definitive assembly quality control.                  | User        |
| `--faa_dir`     | A directory containing protein FASTA (`.faa`) files. For Pharokka users, this is the same path as `--pharokka_dir`.                                                        | Pharokka    |

---

## 4. The GenoSyntax Pipeline

### 4.1. Recommended Usage: The `run` Command

For most users, the integrated `run` command is the recommended entrypoint. It automates the entire workflow, manages all intermediate files, and provides a single, unified interface for all parameters.

**Default Usage (Core Analysis):**
This command executes Steps 1, 2, and 3, culminating in the statistical analysis and the generation of the `comprehensive_summary_report.json`. This is sufficient for most standard analyses.

```bash
python main.py run \
    --panaroo_file /path/to/panaroo/gene_presence_absence_roary.csv \
    --pharokka_dir /path/to/pharokka_outputs/ \
    --fasta_dir /path/to/raw_genome_fastas/ \
    --output_dir ./My_GenoSyntax_Analysis/
```

**Full Pipeline with Mechanism Exploration:**
To run all four steps, including the final mechanism analysis, use the `--mechanisms` flag and set `--until-step 4`.

```bash
python main.py run \
    --panaroo_file /path/to/panaroo/gene_presence_absence_roary.csv \
    --pharokka_dir /path/to/pharokka_outputs/ \
    --fasta_dir /path/to/raw_genome_fastas/ \
    --output_dir ./My_Full_Analysis/ \
    --mechanisms motif evo_dynamics adjacency_analysis extract_layouts \
    --until-step 4
```

### 4.2. Detailed Module Breakdown

#### Module 1: `define_architecture`

**Purpose:** This foundational module identifies and characterizes the stable structural backbone of the genome cohort.

**Key Parameters:**
*   `--max_contigs` (Default: `1`): The maximum number of contigs a genome can have to pass assembly QC.
*   `--tolerance` (Default: `2`): The maximum number of accessory genes allowed between two core genes within an anchor.
*   `--coverage_threshold` (Default: `0.9`): The minimum fraction of high-quality genomes in which an anchor must be found.
*   `--min_length` (Default: `2`): The minimum number of core genes required to form a valid anchor.
*   `--strict_core_threshold` (Default: `1.0`): Prevalence threshold for a gene to be 'strict core'.
*   `--soft_core_threshold` (Default: `0.98`): Prevalence threshold for a gene to be 'soft core', also defines the boundary for accessory genes.

**Key Outputs:**
*   `1_define_architecture/architecture_definition.json`: The master blueprint of the genomic architecture. It contains:
    *   `metadata`: QC summary and all analysis parameters.
    *   `anchor_definitions`: A detailed list of each anchor, including its gene skeleton, coverage, sequence conservation (**ASI**), structural rigidity (**AOI**), and functional class.
    *   `macro_architectures`: A summary of the global arrangement of anchors, including the dominant permutation and overall structural entropy.

#### Module 2: `analyze_slots`

**Purpose:** This module analyzes the plastic regions of the genome (the "slot space"), discovers modular accessory elements, and quantifies their structural dynamics.

**Key Parameters:**
*   `--spi_resolution` (Default: `'low'`): Sets the granularity for SPI calculation. `'low'` resolution is recommended as it aggregates related slots and reduces noise. `'high'` resolution treats every unique inter-core-gene position as distinct.
*   `--phi_threshold` (Default: `0.8`): Phi coefficient threshold for identifying statistical co-occurrence between accessory genes.
*   `--p_value_threshold` (Default: `0.01`): P-value threshold for gene co-occurrence.
*   `--adjacency_distance` (Default: `250`): Maximum distance (in bp) between adjacent genes to be considered part of the same physical module.
*   `--support_threshold` (Default: `0.8`): Minimum fraction of genomes where a candidate module's members must be physically adjacent to be validated.

**Key Outputs:**
*   `2_slot_analysis/reports/module_dynamics_report.csv`: The central report for plasticity analysis. Each row represents an Accessory Module (AM), with key columns:
    *   `spi`: The **Structural Plasticity Index**, indicating how many different structural locations the module occupies.
    *   `frequency`: The prevalence of the module in the population.
    *   `behavioral_class`: A classification based on SPI and frequency (e.g., "Multi-Position", "Structurally Conserved").

#### Module 3: `run_stats`

**Purpose:** To apply rigorous statistical testing to determine if the observed distribution of accessory modules is non-random, thereby identifying insertion "hotspots" and "coldspots".

**Key Parameters:**
*   `--permutations` (Default: `10000`): The number of simulations to perform for the permutation test. Higher numbers increase precision but take longer.

**Key Outputs:**
This module performs a **dual-tier analysis** and generates two key reports:
1.  `3_statistical_test/reports/permutation_test_slot_significance.csv`: The primary report testing for hotspots/coldspots based on the insertion of **any accessory gene**. Use this to understand broad patterns of genomic plasticity.
2.  `3_statistical_test/reports/permutation_test_multigene_module_significance.csv`: An additional report testing specifically for hotspots/coldspots of **multi-gene modules (cassettes)**. This is highly recommended for discovering hotspots likely governed by dedicated insertion mechanisms (e.g., site-specific integrases).

#### Module 4: `explore_mechanisms`

**Purpose:** The final analytical stage, designed to prepare high-quality data for specialized downstream tools to explore the mechanisms underlying the observed architectural patterns.

**Key Parameters:**
*   `--analyses`: A space-separated list of which analyses to run. Choices:
    *   `motif`: Prepare FASTA files of slot boundaries for motif discovery.
    *   `evo_dynamics`: Analyze the relationship between module mobility (SPI) and protein evolutionary rate.
    *   `adjacency_analysis`: Perform a deep functional analysis of multi-gene modules.
    *   `extract_layouts`: Generate a simple coordinate map of all anchors and slots for visualization.
    *   `extract_cds_layout`: Generate a detailed table of every Coding Sequence (CDS) in all high-quality genomes, including coordinates and functional category (CDS layout).
*   `--spi_enrichment_threshold` (Default: `1.0`): For `evo_dynamics`, the SPI value above which a module is considered "highly plastic".

**Key Outputs:**
*   **For `motif` analysis:** Creates two directories:
    *   `4_mechanism_analysis/motif_analysis_all_modules/`: FASTA files of boundaries from general hotspots.
    *   `4_mechanism_analysis/motif_analysis_multigene_modules/`: FASTA files from cassette-specific hotspots (recommended for finding integration motifs).
    *   `4_mechanism_analysis/motif_analysis_high_spi_enriched/`: Contains boundary sequences from slots that are statistically enriched for hosting **highly mobile (high-SPI) modules**. This dataset is ideal for exploratory analysis aimed at identifying signals that might facilitate promiscuous or hyper-mobile element integration. *Note: This output is only generated if the `evo_dynamics` analysis is also run and identifies significant enrichment.*
*   **For `evo_dynamics` analysis:**
    *   `reports/spi_vs_ped_analysis.csv`: A report linking each module's structural plasticity (SPI) to its Protein Evolutionary Distance (PED).
    *   `reports/high_spi_slot_enrichment.csv`: A statistical report identifying slots preferentially targeted by highly plastic modules.
*   **For `adjacency_analysis` analysis:**
    *   `adjacency_module_analysis/`: A directory containing a narrative report, a CSV summary, and network files (`.csv`) for visualizing multi-gene module composition in tools like Cytoscape.
*   **For `extract_layouts` analysis:**
    *   `reports/all_genome_layouts.csv`: A coordinate map of all anchors and slots in every high-quality genome, perfect for input into plotting libraries like `ggplot2`'s `geom_rect`.
*   **For `extract_cds_layout` analysis:**
    *   `reports/all_genomes_cds_layout.csv`: A comprehensive table listing every CDS from all high-quality genomes, with columns for start/end coordinates, strand, product annotation, and a high-level functional category. This file is ideal for detailed gene-level visualizations and custom downstream analyses.


---

## 5. The Comprehensive Summary Report

After the core analysis (Steps 1-3) is complete, GenoSyntax generates a final, high-level summary file named `comprehensive_summary_report.json` in the main output directory. This report provides a quantitative snapshot of the cohort's architectural properties.

### 5.1. Genomic Rigidity Index (GRI)

The GRI is a composite score between **0 (highly fluid) and 1 (highly rigid)** that quantifies the stability of the genomic backbone.

-   **Interpretation:** A score near **1.0** indicates a highly conserved architecture where most core genes are locked into stable anchor modules with a consistent global arrangement. A score near **0** indicates a fluid pangenome with few or no identifiable anchors.
-   **Components:** The GRI is a weighted average of `backbone_coverage`, `backbone_stability`, and `macro_conservation`.

### 5.2. Genomic Plasticity Profile

Instead of a single, weighted score, GenoSyntax provides a **Genomic Plasticity Profile** consisting of three independent metrics. This allows for a more nuanced, multi-dimensional understanding of how plasticity is expressed in a given cohort. Comparing these scores between different datasets can reveal distinct evolutionary strategies.

-   **`genetic_mobility_score` (M_genetic):**
    *   **What it measures:** The population-weighted average mobility (SPI) of all accessory modules.
    *   **Interpretation:** A high score indicates that accessory modules in the population are frequently found in many different genomic locations. This points to a "promiscuous" or highly dynamic accessory genome where functional units are readily exchanged and moved.

-   **`general_insertion_activity_score` (H_activity):**
    *   **What it measures:** The cumulative statistical significance of all insertion hotspots for **any** type of accessory gene (single or multi-gene).
    *   **Interpretation:** A high score signifies that gene insertion is highly non-random and is strongly concentrated in specific genomic regions. This indicates the presence of powerful, broadly targeted insertion mechanisms or selection pressures favoring insertions at these loci.

-   **`cassette_driven_plasticity_score` (C_cassette):**
    *   **What it measures:** The cumulative statistical significance of insertion hotspots specifically for validated **multi-gene modules (cassettes)**.
    *   **Interpretation:** This is a highly specific metric. A high score strongly suggests the activity of dedicated molecular machinery (e.g., site-specific integrases) that recognizes and inserts organized, co-evolving functional cassettes. It provides powerful evidence for structured, mechanism-driven plasticity.

### 5.3. Key Findings Summary
This section provides a quick, human-readable overview of the most important results, including the dominant architecture, the top 5 most plastic modules, and the top 3 insertion hotspots.

---

### **6. Advanced Usage: A Guide to Step-by-Step Execution and Custom Workflows**

While the integrated `run` command is recommended for standard analyses, GenoSyntax's true power lies in its modular design. Executing each module as a separate command provides granular control, making it an ideal approach for debugging, parameter tuning, and creating custom analytical workflows.

In this mode, you are responsible for managing the input and output paths between steps. The fundamental rule is that the output of one step serves as the primary input for the next.

**General Help:**
You can get a detailed list of parameters for any command using the `-h` or `--help` flag.
```bash
# General help for all commands
python main.py -h

# Specific help for the define_architecture module
python main.py define_architecture -h
```

---

#### **Step 1: `define_architecture`**

**Purpose:** The foundational step that identifies the conserved structural backbone and produces the master `architecture_definition.json` blueprint.

**Command Template:**
```bash
python main.py define_architecture \
    --panaroo_file /path/to/gene_presence_absence_roary.csv \
    --pharokka_dir /path/to/pharokka_outputs/ \
    --fasta_dir /path/to/raw_fastas/ \
    --output_dir ./output/1_define_architecture/ \
    --coverage_threshold 0.9 \
    --tolerance 2
```

*   **Key Inputs:** Raw data files (`--panaroo_file`, `--pharokka_dir`, `--fasta_dir`).
*   **Key Outputs:**
    *   `./output/1_define_architecture/architecture_definition.json` (The master blueprint).
    *   `./output/1_define_architecture/intermediate/` (Directory containing pre-processed databases like `gene_order_database.json`, which are essential for the next step).

*   **Advanced Use Case:** This module is ideal for iterative parameter tuning. If the initial run with default settings yields no anchors, you can re-run *only this step* with a lower `--coverage_threshold` (e.g., 0.8) to test for a less conserved, "dominant sub-structure" backbone without having to repeat any downstream analysis.

---

#### **Step 2: `analyze_slots`**

**Purpose:** To dissect the "plastic" regions of the genome, discover mobile accessory modules (AMs), and quantify their mobility (SPI).

**Command Template:**
```bash
python main.py analyze_slots \
    --architecture_file ./output/1_define_architecture/architecture_definition.json \
    --intermediate_dir_prev ./output/1_define_architecture/intermediate/ \
    --output_dir ./output/2_slot_analysis/ \
    --spi_resolution low
```

*   **Key Inputs:** The `architecture_definition.json` file and the `intermediate/` directory from the `define_architecture` step.
*   **Key Outputs:**
    *   `./output/2_slot_analysis/reports/module_dynamics_report.csv` (The central plasticity report).
    *   `./output/2_slot_analysis/intermediate/` (Contains the `slot_database.json` and `accessory_module_definitions.json`).

*   **Advanced Use Case:** You can experiment with different SPI calculation granularities. Run this step once with `--spi_resolution low` (default) for a broad overview, then re-run it with `--spi_resolution high` to investigate fine-grained positional shifts of specific modules.

---

#### **Step 3: `run_stats`**

**Purpose:** To apply rigorous permutation testing to determine if the observed distribution of AMs is statistically significant.

**Command Template:**
```bash
python main.py run_stats \
    --slot_analysis_dir ./output/2_slot_analysis/ \
    --output_dir ./output/3_statistical_test/ \
    --permutations 10000
```

*   **Key Inputs:** The entire output directory from the `analyze_slots` step.
*   **Key Outputs:**
    *   `./output/3_statistical_test/reports/permutation_test_slot_significance.csv`
    *   `./output/3_statistical_test/reports/permutation_test_multigene_module_significance.csv`

*   **Advanced Use Case:** For a final, publication-quality analysis, you might want to increase the number of simulations for higher confidence in your p-values. You can easily re-run just this step with `--permutations 100000` after completing the initial, faster exploration.

---

#### **Step 4: `explore_mechanisms`**

**Purpose:** To synthesize all prior findings and prepare high-quality data for specialized downstream tools (MEME, Cytoscape, etc.).

**Command Template:**
```bash
python main.py explore_mechanisms \
    --analyses motif evo_dynamics extract_layouts \
    --def_arch_dir ./output/1_define_architecture/ \
    --slot_analysis_dir ./output/2_slot_analysis/ \
    --stats_dir ./output/3_statistical_test/ \
    --pharokka_dir /path/to/pharokka_outputs/ \
    --output_dir ./output/4_mechanism_analysis/
```

*   **Key Inputs:** The output directories from all three previous steps are required, as this module integrates data from the entire pipeline. The raw `pharokka_dir` is also needed for sequence extraction.
*   **Key Outputs:** A series of directories and files within `./output/4_mechanism_analysis/` tailored to the requested `--analyses` (e.g., `motif_analysis_multigene_modules/`, `reports/all_genome_layouts.csv`).

*   **Advanced Use Case:** This module is designed for selective execution. After a full analysis, you might realize you also want to analyze adjacency modules. You can simply run `python main.py explore_mechanisms --analyses adjacency_analysis ...` without having to re-compute anything from steps 1-3.

---

## 7. Troubleshooting

-   **`mafft: command not found`**: This is the most common error. It means the `mafft` executable is not in your system's PATH. Ensure you have activated the correct conda environment (`conda activate genosyntax`) where `mafft` was installed.
-   **File Not Found Errors**: Double-check all input paths. The `run` command's `--dry_run` flag is excellent for validating paths before starting a long computation.
-   **No Anchors Found**: If the `define_architecture` step completes but finds no anchors, this is a valid biological result indicating a highly fluid genome cohort. The pipeline will gracefully skip downstream steps that depend on anchors. Consider re-running with a lower `--coverage-threshold` (e.g., 0.7 or 0.6) to search for less-conserved "dominant sub-structure backbones."
```
