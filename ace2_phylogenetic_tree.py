"""
ACE2 Cross-Species Phylogenetic Tree Builder
---------------------------------------------
Compares the ACE2 gene (the human receptor used by SARS-CoV-2 to enter
cells) across multiple species and builds a family tree showing how
similar or different each species' version is.

How it works:
1. Reads a multi-sequence FASTA file (several species' ACE2 sequences).
2. Aligns every pair of sequences and calculates a similarity score.
3. Uses those similarity scores to build a phylogenetic (family) tree
   using the UPGMA clustering method.
4. Saves the tree as both a text diagram and a PNG image.

Usage:
    python ace2_phylogenetic_tree.py ace2_multi_species.fasta

Input:
    A single FASTA file containing MULTIPLE sequences, one per species.
    Each sequence's header (the line starting with ">") should clearly
    name the species, e.g.:
        >Human_ACE2
        ATGTCG...
        >Cat_ACE2
        ATGTCC...

Output:
    - tree.txt   : a simple text-based diagram of the tree
    - tree.png   : an image of the tree
"""

import sys
from Bio import SeqIO, Phylo
from Bio.Align import PairwiseAligner
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor
import matplotlib
matplotlib.use("Agg")  # allows saving images without needing a display
import matplotlib.pyplot as plt


def load_sequences(filepath):
    """Read all sequences from a multi-FASTA file. Returns a list of (name, sequence)."""
    records = list(SeqIO.parse(filepath, "fasta"))
    if len(records) < 3:
        print(f"Found only {len(records)} sequence(s). You need at least 3 species "
              f"to build a meaningful tree.")
        sys.exit(1)
    return records


def calculate_distance_matrix(records):
    """
    Compare every pair of sequences and calculate a 'distance' between them
    (0 = identical, higher = more different). This uses simple pairwise
    alignment scoring, converted into a percentage-difference distance.
    """
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 1
    aligner.mismatch_score = 0
    aligner.open_gap_score = -1
    aligner.extend_gap_score = -0.5

    names = [record.id for record in records]
    n = len(records)

    # Build a lower-triangle distance matrix, the format Biopython expects
    matrix = []
    for i in range(n):
        row = []
        for j in range(i + 1):
            if i == j:
                row.append(0.0)
            else:
                seq1 = str(records[i].seq)
                seq2 = str(records[j].seq)
                score = aligner.score(seq1, seq2)
                max_possible = max(len(seq1), len(seq2))
                similarity = score / max_possible
                distance = round(1 - similarity, 4)  # convert similarity to distance
                row.append(distance)
        matrix.append(row)

    print("Pairwise comparisons complete.\n")
    return DistanceMatrix(names, matrix)


def build_tree(distance_matrix):
    """Build a phylogenetic tree from the distance matrix using UPGMA clustering."""
    constructor = DistanceTreeConstructor()
    tree = constructor.upgma(distance_matrix)
    return tree


def main():
    if len(sys.argv) != 2:
        print("Usage: python ace2_phylogenetic_tree.py <multi_fasta_file>")
        sys.exit(1)

    fasta_path = sys.argv[1]
    print(f"Reading sequences from {fasta_path}...\n")
    records = load_sequences(fasta_path)

    print(f"Loaded {len(records)} species:")
    for r in records:
        print(f"  - {r.id} ({len(r.seq)} bases)")
    print()

    distance_matrix = calculate_distance_matrix(records)

    print("Building phylogenetic tree...\n")
    tree = build_tree(distance_matrix)

    # Save a simple text diagram of the tree
    with open("tree.txt", "w") as f:
        Phylo.draw_ascii(tree, file=f)
    print("Text tree diagram saved to tree.txt")

    # Also print it to the screen
    print("\nYour tree:\n")
    Phylo.draw_ascii(tree)

    # Save a picture of the tree
    fig = plt.figure(figsize=(10, 6))
    axes = fig.add_subplot(1, 1, 1)
    Phylo.draw(tree, axes=axes, do_show=False)
    plt.savefig("tree.png", dpi=150, bbox_inches="tight")
    print("\nTree image saved to tree.png")


if __name__ == "__main__":
    main()
