
# ============================================================
# PROJECT: Comparative Genomics of TP53 Across Mammals
# PURPOSE: Analyse evolutionary conservation of TP53 proteins
# DATABASE: NCBI Datasets
# TOOLS: Biopython, NumPy, Pandas, Matplotlib, scikit-learn,
#        MAFFT and FastTree
# ============================================================

# -------------------- IMPORT LIBRARIES ----------------------

import re
import zipfile
import subprocess
import shutil
import argparse

from pathlib import Path
from datetime import datetime, timezone
from collections import Counter

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from Bio import SeqIO, AlignIO, Phylo
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA


# -------------------- PROJECT PATHS -------------------------

ROOT = Path(__file__).resolve().parents[1]

RAW_DIR = ROOT / "data" / "raw"
EXTRACT_DIR = RAW_DIR / "extracted"
PROCESSED_DIR = ROOT / "data" / "processed"

FIGURE_DIR = ROOT / "results" / "figures"
TABLE_DIR = ROOT / "results" / "tables"
ALIGNMENT_DIR = ROOT / "results" / "alignment"
TREE_DIR = ROOT / "results" / "phylogeny"

ZIP_FILE = RAW_DIR / "tp53_mammals.zip"

# Create folders if they do not exist
for folder in [
    RAW_DIR, EXTRACT_DIR, PROCESSED_DIR,
    FIGURE_DIR, TABLE_DIR, ALIGNMENT_DIR, TREE_DIR
]:
    folder.mkdir(parents=True, exist_ok=True)


# Standard amino acids
AMINO_ACIDS = list("ACDEFGHIKLMNPQRSTVWY")
VALID_AA = set(AMINO_ACIDS + ["X"])

# Preferred species selection order
PREFERRED_SPECIES = [
    "Homo sapiens",
    "Pan troglodytes",
    "Macaca mulatta",
    "Mus musculus",
    "Rattus norvegicus",
    "Canis lupus familiaris",
    "Bos taurus",
    "Sus scrofa",
    "Equus caballus"
]

MAX_SPECIES = 10


# -------------------- STEP 1: DOWNLOAD DATA -----------------

def download_data(refresh=False):

    if ZIP_FILE.exists() and not refresh:
        print("NCBI dataset already downloaded.")
        return

    if shutil.which("datasets") is None:
        raise RuntimeError(
            "NCBI Datasets CLI not found. Please install datasets."
        )

    print("\nDownloading TP53 mammalian orthologues from NCBI...")

    command = [
        "datasets",
        "download",
        "gene",
        "gene-id",
        "7157",
        "--ortholog",
        "mammals",
        "--include",
        "protein",
        "--filename",
        str(ZIP_FILE)
    ]

    subprocess.run(command, check=True)

    print("NCBI download completed.")


# -------------------- STEP 2: EXTRACT DATA -------------------

def extract_data():

    protein_file = EXTRACT_DIR / "protein.faa"

    if protein_file.exists():
        print("Extracted FASTA file already exists.")
        return protein_file

    print("\nExtracting downloaded ZIP file...")

    with zipfile.ZipFile(ZIP_FILE, "r") as archive:
        archive.extractall(EXTRACT_DIR)

    # Locate the protein FASTA file
    matches = list(EXTRACT_DIR.rglob("protein.faa"))

    if not matches:
        raise FileNotFoundError(
            "protein.faa not found in NCBI dataset package."
        )

    protein_file = matches[0]

    print("Protein FASTA found:", protein_file)

    return protein_file


# -------------------- STEP 3: PARSE FASTA --------------------

def parse_fasta(protein_file):

    print("\nReading protein sequences...")

    records = list(SeqIO.parse(protein_file, "fasta"))

    if not records:
        raise ValueError("No sequences found in FASTA file.")

    print("Total downloaded sequences:", len(records))

    return records


# -------------------- STEP 4: DATA CLEANING -----------------

def clean_sequences(records):

    print("\nPerforming sequence quality checks...")

    species_sequences = {}

    quality_data = []

    for record in records:

        header = record.description

        # Extract species name from FASTA header
        match = re.search(
            r"\[organism=([^\]]+)\]",
            header
        )

        if not match:
            continue

        species = match.group(1)

        # Convert sequence to uppercase
        sequence = str(record.seq).upper()

        # Remove terminal stop symbol
        sequence = sequence.rstrip("*")

        # Check sequence length
        if len(sequence) < 200 or len(sequence) > 600:
            continue

        # Check for unexpected amino acid characters
        if not set(sequence).issubset(VALID_AA):
            continue

        # Calculate ambiguous amino acid percentage
        ambiguous_count = sequence.count("X")

        ambiguous_percentage = (
            ambiguous_count / len(sequence)
        ) * 100

        # Exclude sequences with excessive ambiguity
        if ambiguous_percentage > 5:
            continue

        quality_data.append({
            "species": species,
            "accession": record.id,
            "sequence_length": len(sequence),
            "ambiguous_percentage": ambiguous_percentage
        })

        # Keep the longest protein isoform for each species
        if (
            species not in species_sequences
            or len(sequence) > len(species_sequences[species].seq)
        ):

            species_sequences[species] = SeqRecord(
                Seq(sequence),
                id=record.id,
                description=header
            )

    quality_df = pd.DataFrame(quality_data)

    quality_df.to_csv(
        TABLE_DIR / "sequence_quality.csv",
        index=False
    )

    print("Species passing quality checks:",
          len(species_sequences))

    return species_sequences


# -------------------- STEP 5: SELECT SPECIES ---------------

def select_species(species_sequences):

    print("\nSelecting representative species...")

    selected = []

    # Select preferred species first
    for species in PREFERRED_SPECIES:

        if species in species_sequences:
            selected.append(
                (species, species_sequences[species])
            )

    # If fewer than 4 preferred species are found,
    # supplement with other available mammalian orthologues.
    if len(selected) < 4:

        already_selected = {
            species for species, record in selected
        }

        for species in sorted(species_sequences):

            if species not in already_selected:
                selected.append(
                    (species, species_sequences[species])
                )

    selected = selected[:MAX_SPECIES]

    if len(selected) < 4:
        raise ValueError(
            "Fewer than four usable species found. "
            "Check the NCBI orthologue download."
        )

    cleaned_records = []
    metadata = []

    for species, record in selected:

        accession = record.id

        # Use simple, machine-readable FASTA identifiers
        safe_species = species.replace(" ", "_")

        record.id = f"{safe_species}__{accession}"
        record.name = record.id
        record.description = record.id

        cleaned_records.append(record)

        sequence = str(record.seq)

        metadata.append({
            "species": species,
            "accession": accession,
            "sequence_length": len(sequence),
            "ambiguous_percentage":
                sequence.count("X") / len(sequence) * 100
        })

    # Save cleaned FASTA
    cleaned_fasta = PROCESSED_DIR / "tp53_cleaned.faa"

    SeqIO.write(
        cleaned_records,
        cleaned_fasta,
        "fasta"
    )

    metadata_df = pd.DataFrame(metadata)

    metadata_df.to_csv(
        TABLE_DIR / "sequence_metadata.csv",
        index=False
    )

    print("\nSelected species:")

    print(metadata_df[["species", "sequence_length"]])

    return cleaned_records, metadata_df


# -------------------- STEP 6: MULTIPLE SEQUENCE ALIGNMENT --

def perform_alignment():

    print("\nPerforming multiple sequence alignment...")

    mafft = shutil.which("mafft")

    if mafft is None:
        raise RuntimeError("MAFFT is not installed.")

    input_file = PROCESSED_DIR / "tp53_cleaned.faa"

    output_file = ALIGNMENT_DIR / "tp53_aligned.faa"

    command = [
        mafft,
        "--auto",
        "--thread",
        "2",
        str(input_file)
    ]

    result = subprocess.run(
        command,
        capture_output=True,
        text=True,
        check=True
    )

    output_file.write_text(result.stdout)

    alignment = AlignIO.read(
        output_file,
        "fasta"
    )

    print("Alignment completed.")
    print("Number of sequences:", len(alignment))
    print("Alignment length:", alignment.get_alignment_length())

    return alignment


# -------------------- STEP 7: PAIRWISE IDENTITY --------------

def calculate_identity(alignment):

    print("\nCalculating pairwise sequence identity...")

    sequences = list(alignment)

    species = [
        record.id.split("__")[0]
        for record in sequences
    ]

    n = len(sequences)

    matrix = np.zeros((n, n))

    for i in range(n):

        for j in range(n):

            seq1 = str(sequences[i].seq)
            seq2 = str(sequences[j].seq)

            matches = 0
            compared = 0

            for aa1, aa2 in zip(seq1, seq2):

                # Ignore gaps and ambiguous residues
                if aa1 not in AMINO_ACIDS:
                    continue

                if aa2 not in AMINO_ACIDS:
                    continue

                compared += 1

                if aa1 == aa2:
                    matches += 1

            identity = (
                matches / compared * 100
                if compared > 0 else np.nan
            )

            matrix[i, j] = identity

    identity_df = pd.DataFrame(
        matrix,
        index=species,
        columns=species
    )

    identity_df.to_csv(
        TABLE_DIR / "pairwise_identity.csv"
    )

    return identity_df


# -------------------- STEP 8: CONSERVED POSITIONS -----------

def analyse_conservation(alignment):

    print("\nAnalysing conserved amino acid positions...")

    sequences = [
        str(record.seq)
        for record in alignment
    ]

    alignment_length = len(sequences[0])

    conservation_data = []

    for position in range(alignment_length):

        column = [
            sequence[position]
            for sequence in sequences
            if sequence[position] in AMINO_ACIDS
        ]

        if not column:
            continue

        counts = Counter(column)

        most_common_aa, count = counts.most_common(1)[0]

        conservation = (
            count / len(column)
        ) * 100

        conservation_data.append({
            "position": position + 1,
            "most_common_amino_acid": most_common_aa,
            "conservation_percentage": conservation,
            "number_of_valid_sequences": len(column)
        })

    conservation_df = pd.DataFrame(
        conservation_data
    )

    conservation_df.to_csv(
        TABLE_DIR / "conservation_analysis.csv",
        index=False
    )

    return conservation_df


# -------------------- STEP 9: AMINO ACID COMPOSITION --------

def analyse_composition(records, metadata):

    print("\nCalculating amino acid composition...")

    composition = []

    for record in records:

        sequence = str(record.seq)

        counts = Counter(sequence)

        total = sum(
            counts[aa] for aa in AMINO_ACIDS
        )

        row = {
            "species": record.id.split("__")[0]
                       .replace("_", " ")
        }

        for aa in AMINO_ACIDS:

            row[aa] = (
                counts[aa] / total * 100
                if total > 0 else 0
            )

        composition.append(row)

    composition_df = pd.DataFrame(composition)

    composition_df.to_csv(
        TABLE_DIR / "amino_acid_composition.csv",
        index=False
    )

    return composition_df


# -------------------- STEP 10: EXPLORATORY PCA --------------

def perform_pca(composition_df):

    print("\nPerforming PCA...")

    species = composition_df["species"]

    # Use the 20 amino acid percentage columns
    X = composition_df[AMINO_ACIDS].values

    # Standardise features
    scaler = StandardScaler()

    X_scaled = scaler.fit_transform(X)

    # Reduce data to two principal components
    pca = PCA(n_components=2)

    X_pca = pca.fit_transform(X_scaled)

    pca_df = pd.DataFrame({
        "species": species,
        "PC1": X_pca[:, 0],
        "PC2": X_pca[:, 1]
    })

    pca_df.to_csv(
        TABLE_DIR / "pca_results.csv",
        index=False
    )

    # Plot PCA
    plt.figure(figsize=(9, 6))

    plt.scatter(
        X_pca[:, 0],
        X_pca[:, 1],
        s=100
    )

    for i, name in enumerate(species):

        plt.annotate(
            name.replace(" ", "\n"),
            (X_pca[i, 0], X_pca[i, 1]),
            fontsize=8
        )

    variance = pca.explained_variance_ratio_ * 100

    plt.xlabel(f"PC1 ({variance[0]:.2f}% variance)")
    plt.ylabel(f"PC2 ({variance[1]:.2f}% variance)")

    plt.title("PCA of TP53 Amino Acid Composition")

    plt.tight_layout()

    plt.savefig(
        FIGURE_DIR / "pca_composition.png",
        dpi=300
    )

    plt.close()

    return pca_df


# -------------------- STEP 11: PHYLOGENETIC TREE ------------

def construct_tree():

    print("\nConstructing phylogenetic tree...")

    fasttree = shutil.which("FastTree") or shutil.which("fasttree")

    if fasttree is None:
        raise RuntimeError("FastTree is not installed.")

    alignment_file = ALIGNMENT_DIR / "tp53_aligned.faa"

    tree_file = TREE_DIR / "tp53_tree.nwk"

    command = [
        fasttree,
        "-lg",
        str(alignment_file)
    ]

    result = subprocess.run(
        command,
        capture_output=True,
        text=True,
        check=True
    )

    tree_file.write_text(result.stdout)

    tree = Phylo.read(
        tree_file,
        "newick"
    )

    plt.figure(figsize=(12, 8))

    Phylo.draw(
        tree,
        do_show=False
    )

    plt.title("Phylogenetic Tree of TP53 Mammalian Orthologues")

    plt.tight_layout()

    plt.savefig(
        FIGURE_DIR / "tp53_phylogenetic_tree.png",
        dpi=300,
        bbox_inches="tight"
    )

    plt.close()

    print("Phylogenetic tree generated.")

    return tree


# -------------------- STEP 12: VISUALISATION ----------------

def create_visualisations(identity_df, conservation_df):

    print("\nGenerating visualisations...")

    # Pairwise identity heatmap
    plt.figure(figsize=(10, 8))

    plt.imshow(
        identity_df.values,
        cmap="viridis",
        vmin=0,
        vmax=100
    )

    plt.colorbar(
        label="Sequence identity (%)"
    )

    labels = [
        name.replace("_", " ")
        for name in identity_df.index
    ]

    plt.xticks(
        range(len(labels)),
        labels,
        rotation=90
    )

    plt.yticks(
        range(len(labels)),
        labels
    )

    plt.title("TP53 Pairwise Sequence Identity")

    plt.tight_layout()

    plt.savefig(
        FIGURE_DIR / "pairwise_identity_heatmap.png",
        dpi=300
    )

    plt.close()

    # Conservation plot
    plt.figure(figsize=(12, 5))

    plt.plot(
        conservation_df["position"],
        conservation_df["conservation_percentage"]
    )

    plt.xlabel("Alignment Position")
    plt.ylabel("Conservation (%)")

    plt.title("TP53 Amino Acid Conservation")

    plt.ylim(0, 105)

    plt.tight_layout()

    plt.savefig(
        FIGURE_DIR / "conservation_plot.png",
        dpi=300
    )

    plt.close()


# -------------------- MAIN PIPELINE -------------------------

def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--refresh",
        action="store_true",
        help="Download a fresh NCBI dataset"
    )

    args = parser.parse_args()

    print("=" * 55)
    print("TP53 COMPARATIVE GENOMICS PIPELINE")
    print("=" * 55)

    download_data(refresh=args.refresh)

    protein_file = extract_data()

    records = parse_fasta(protein_file)

    species_sequences = clean_sequences(records)

    cleaned_records, metadata = select_species(
        species_sequences
    )

    alignment = perform_alignment()

    identity_df = calculate_identity(alignment)

    conservation_df = analyse_conservation(alignment)

    composition_df = analyse_composition(
        cleaned_records,
        metadata
    )

    perform_pca(composition_df)

    construct_tree()

    create_visualisations(
        identity_df,
        conservation_df
    )

    # Save project summary
    summary = f"""
TP53 COMPARATIVE GENOMICS PROJECT

Database: NCBI Datasets
NCBI Gene ID: 7157
Taxonomic scope: Mammals

Date of analysis (UTC):
{datetime.now(timezone.utc).isoformat()}

Number of selected species:
{len(cleaned_records)}

Number of aligned sequences:
{len(alignment)}

Alignment length:
{alignment.get_alignment_length()}

Main analysis:
- Sequence quality control
- Pairwise sequence identity
- Amino acid conservation
- Amino acid composition
- Principal Component Analysis
- Phylogenetic tree construction

External tools:
- MAFFT
- FastTree
"""

    (TABLE_DIR / "analysis_summary.txt").write_text(summary)

    print("\n" + "=" * 55)
    print("ANALYSIS COMPLETED SUCCESSFULLY!")
    print("=" * 55)

    print("\nCheck the results/ directory for output files.")


if __name__ == "__main__":
    main()
