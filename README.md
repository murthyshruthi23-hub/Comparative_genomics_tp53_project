# Comparative Sequence Analysis and Evolutionary Conservation of TP53 Orthologs Across Mammalian Species

## 1. Project Overview

This project investigates the evolutionary conservation of the TP53 tumour suppressor protein across different mammalian species using computational bioinformatics approaches.

TP53 is a well-characterised tumour suppressor gene that encodes the p53 protein, which plays an important role in regulating cell cycle progression, DNA damage responses, cellular senescence and apoptosis.

The objective of this project is to compare TP53 protein sequences from different mammalian species, identify conserved and variable amino acid positions, calculate pairwise sequence identity, and investigate evolutionary relationships using phylogenetic analysis.

The project uses publicly available sequence data from the National Center for Biotechnology Information (NCBI) and implements a reproducible computational pipeline using Python and established bioinformatics tools.

## 2. Research Question

How conserved is the TP53 protein across different mammalian species, and what can sequence similarity, amino acid conservation and phylogenetic relationships reveal about its evolutionary importance?

## 3. Research Objectives

- Retrieve TP53 orthologous protein sequences from different mammalian species using NCBI.
- Perform FASTA parsing, sequence preprocessing and quality control.
- Select representative protein sequences for comparative analysis.
- Perform multiple sequence alignment using MAFFT.
- Calculate pairwise amino acid sequence identity between species.
- Identify conserved and variable amino acid positions.
- Analyse amino acid composition across species.
- Perform Principal Component Analysis (PCA) to explore compositional differences.
- Construct a protein phylogenetic tree using FastTree.
- Visualise and interpret the results in an evolutionary context.

## 4. Biological Background

TP53 is a tumour suppressor gene located on chromosome 17 in humans. It encodes the p53 protein, commonly described as the "guardian of the genome" because of its role in maintaining genomic stability.

In response to cellular stress, DNA damage and other signals, p53 can regulate the expression of downstream genes involved in:

- Cell cycle arrest
- DNA repair
- Apoptosis
- Cellular senescence

Evolutionary conservation of TP53 across mammals provides an opportunity to investigate how important regions of this protein have been maintained throughout evolution.

Comparative sequence analysis can help identify conserved amino acid positions that may be subject to evolutionary constraints.

## 5. Data Source and Retrieval

### Database

National Center for Biotechnology Information (NCBI)

### Reference Gene

- Gene: TP53
- Human NCBI Gene ID: 7157
- Data type: Protein sequences
- Format: FASTA
- Taxonomic scope: Mammalia

### Data Retrieval Method

TP53 orthologous protein sequences were retrieved using the NCBI Datasets command-line interface.

The following command was used:

```bash
datasets download gene gene-id 7157 \
--ortholog mammals \
--include protein \
--filename data/raw/tp53_mammals.zip
```

The downloaded dataset was extracted, and protein sequences were retrieved from the FASTA file for downstream analysis.

### Species Analysed

- Homo sapiens
- Pan troglodytes
- Macaca mulatta
- Mus musculus
- Rattus norvegicus
- Canis lupus familiaris
- Bos taurus
- Sus scrofa
- Equus caballus


## 6. Tools and Technologies

| Tool / Package | Application |
|---|---|
| Python 3 | Main programming language |
| Biopython | FASTA parsing and biological sequence processing |
| NumPy | Numerical calculations |
| Pandas | Data manipulation and analysis |
| Matplotlib | Data visualisation |
| Scikit-learn | Principal Component Analysis |
| NCBI Datasets CLI | Retrieval of orthologous sequences |
| MAFFT | Multiple sequence alignment |
| FastTree | Phylogenetic tree construction |
| Git and GitHub | Version control and project hosting |

## 7. Methodology and Computational Workflow

### Step 1: Data Retrieval

TP53 orthologous protein sequences were retrieved from NCBI using the gene-specific orthologue download functionality of NCBI Datasets.

The dataset was restricted to protein sequences to minimise storage requirements and computational complexity.

### Step 2: Sequence Preprocessing and Quality Control

Protein FASTA sequences were processed using Biopython.

The following preprocessing operations were performed:

- Conversion of amino acid sequences to uppercase.
- Removal of terminal stop symbols.
- Filtering based on sequence length.
- Validation of amino acid characters.
- Identification of ambiguous amino acid residues represented by X.
- Exclusion of sequences with excessive ambiguity.
- Selection of the longest available protein sequence per species as a representative sequence.

Sequence metadata and quality-control information were exported as CSV files.

### Step 3: Multiple Sequence Alignment

Multiple sequence alignment was performed using MAFFT.

The alignment was used to compare homologous amino acid positions across species while accounting for insertions and deletions.

### Step 4: Pairwise Sequence Identity Analysis

Pairwise amino acid sequence identity was calculated for all combinations of aligned protein sequences.

Sequence identity was calculated as:

Identity (%) = (Number of identical amino acid positions / Number of comparable positions) × 100

Alignment gaps and ambiguous amino acid residues were excluded from the comparison.

A pairwise identity matrix was generated and visualised as a heatmap.

### Step 5: Amino Acid Conservation Analysis

Conservation was assessed at each alignment position by identifying the most frequently occurring amino acid.

The percentage conservation was calculated as:

Conservation (%) = (Frequency of the most common amino acid / Number of valid amino acids at that position) × 100

The resulting values were used to identify highly conserved and variable regions.

### Step 6: Amino Acid Composition Analysis

The relative frequency of each of the 20 standard amino acids was calculated for each selected TP53 protein sequence.

The results were represented as percentage composition values and exported into a CSV file.

### Step 7: Principal Component Analysis

PCA was performed using Scikit-learn on the amino acid composition profiles.

The 20 amino acid composition features were standardised before dimensionality reduction.

The first two principal components were visualised to explore differences in amino acid composition between species.

### Step 8: Phylogenetic Analysis

The multiple sequence alignment was provided to FastTree for protein phylogenetic inference using the LG amino acid substitution model.

The resulting phylogenetic tree was saved in Newick format and visualised using Biopython.


## 8. Results and Visualisations

The computational pipeline generates the following outputs:

### 8.1 Pairwise Sequence Identity Heatmap

File: `results/figures/pairwise_identity_heatmap.png`

This heatmap represents percentage amino acid sequence identity between pairs of mammalian TP53 proteins.

It helps identify species pairs with higher or lower sequence identity.

### 8.2 Amino Acid Conservation Plot

File: `results/figures/conservation_plot.png`

This plot represents amino acid conservation across aligned TP53 protein positions.

It helps identify regions that are highly conserved or variable across the selected species.

### 8.3 PCA of Amino Acid Composition

File: `results/figures/pca_composition.png`

This visualisation represents differences in amino acid composition between species using the first two principal components.

It provides an exploratory representation of compositional similarities and differences.

### 8.4 TP53 Phylogenetic Tree

File: `results/figures/tp53_phylogenetic_tree.png`

This tree illustrates evolutionary relationships inferred from the aligned TP53 protein sequences.

The phylogenetic tree is based on a single protein and should not be interpreted as a complete species phylogeny.

### 8.5 Output Tables

The pipeline also generates:

| Output file | Description |
|---|---|
| `sequence_metadata.csv` | Selected species and sequence metadata |
| `sequence_quality.csv` | Sequence quality-control information |
| `pairwise_identity.csv` | Pairwise identity matrix |
| `conservation_analysis.csv` | Conservation percentages at alignment positions |
| `amino_acid_composition.csv` | Amino acid composition profiles |
| `pca_results.csv` | Principal component coordinates |
| `analysis_summary.txt` | Summary of the analysis |


## 9. Repository Structure

```text
comparative-genomics-tp53_project/
│
├── README.md
├── requirements.txt
├── .gitignore
│
├── data/
│   ├── raw/
│   └── processed/
│
├── scripts/
│   └── run_analysis.py
│
└── results/
    ├── alignment/
    ├── figures/
    ├── tables/
    └── phylogeny/
```

## 10. Installation and Reproducibility

### Clone the repository

```bash
git clone https://github.com/murthyshruthi23-hub/comparative-genomics-tp53_project.git

cd comparative-genomics-tp53_project
```

### Create a Python virtual environment

```bash
python3 -m venv .venv

source .venv/bin/activate
```

### Install Python dependencies

```bash
pip install -r requirements.txt
```

### Install external bioinformatics tools

On Ubuntu/Debian-based Linux:

```bash
sudo apt update

sudo apt install mafft fasttree
```

Install the NCBI Datasets CLI separately using the official NCBI installation instructions.

### Execute the pipeline

```bash
python scripts/run_analysis.py
```

To retrieve a fresh dataset and rerun the analysis:

```bash
python scripts/run_analysis.py --refresh
```

The pipeline retrieves the required sequences, performs preprocessing, executes the analysis and generates the result files.

## 11. Limitations

- The study focuses on a limited number of mammalian TP53 orthologues rather than complete genomes.
- Selection of the longest available protein isoform does not necessarily guarantee selection of the canonical isoform.
- Sequence identity measures exact amino acid matches and does not account for biochemical similarity between different amino acids.
- Conservation does not independently establish the functional importance of a particular amino acid position.
- PCA is based on amino acid composition and does not capture the complete positional arrangement of amino acids.
- The phylogenetic tree represents relationships inferred from a single protein rather than genome-wide evolutionary relationships.
- Phylogenetic inference depends on sequence alignment quality and the selected substitution model.

## 12. Conclusion

This project provides a computational framework for investigating evolutionary conservation of TP53 across mammalian species.

By combining sequence retrieval, preprocessing, multiple sequence alignment, pairwise identity calculations, amino acid conservation analysis, PCA and phylogenetic inference, the project demonstrates the application of computational bioinformatics techniques to comparative genomics.

The workflow also provides a reproducible approach for analysing other conserved genes and exploring evolutionary relationships across species.

## 13. References and Resources

1. NCBI Gene – TP53: https://www.ncbi.nlm.nih.gov/gene/7157
2. NCBI Datasets: https://www.ncbi.nlm.nih.gov/datasets/
3. Biopython Documentation: https://biopython.org/wiki/Documentation
4. MAFFT: https://mafft.cbrc.jp/alignment/software/
5. FastTree: http://www.microbesonline.org/fasttree/
6. Scikit-learn PCA Documentation: https://scikit-learn.org/stable/modules/decomposition.html#pca

---

**Project Type:** Bioinformatics / Comparative Genomics

**Programming Language:** Python

**Data Source:** NCBI

**Repository:** GitHub