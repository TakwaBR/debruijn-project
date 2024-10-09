# Assembleur basé sur les graphes de Debruijn

Vous trouverez la description complète du TP [ici]( 
https://docs.google.com/document/d/1P4v3bHbSurD7RXA-ldVwmtNKGvWsnBae51RMGye_KLs/edit?usp=sharing).

## Introduction

The objective of this project is to assemble the genome of Enterovirus A71. This genome is of particular interest because it is very short: 7408 nucleotides, linear, and non-segmented. The FASTQ file we have was generated using the ART program [Huang 2011]. 

In the debruijn-tp/data/ trajectory, we have:
eva71.fna : genome of the virus of interest
eva71_plus_perfect.fq: reads 


## Installation of Dependencies and Requirements

We will use:
- Python 3.x
- Required libraries:
  - `argparse`
  - `os`
  - `sys`
  - `statistics`
  - `textwrap`
  - `networkx`
  - `matplotlib

You can install the required libraries using pip:

```bash
pip install networkx matplotlib
```

## Utilisation

To run the program, use the following command:

```
python3 debruijn.py -i <input.fastq> -k <kmer_size> -o <output.fasta> -f <graph_image.png>
```
## Arguments
- -i : Path to the input FASTQ file (required).
- -k : K-mer size (default is 22).
- -o : Path to the output file for contigs in FASTA format (default is ../results/contigs.fasta).
- -f : Path to save the graph as an image (PNG format).

## Description

The program performs the following steps:

1. **Read FASTQ File**: Extracts reads from the specified FASTQ file.
2. **Build K-mer Dictionary**: Constructs a dictionary of all k-mer occurrences in the FASTQ file.
3. **Build De Bruijn Graph**: Creates a directed graph of all k-mer substrings and their weights.
4. **Simplify Bubbles**: Detects and resolves bubbles in the graph.
5. **Solve Entry and Out Tips**: Removes entry tips and out tips from the graph.
6. **Extract Contigs**: Generates contiguous sequences from the graph.
7. **Save Contigs**: Writes the contigs to a specified output file in FASTA format.
8. **Draw Graph**: (optional) Saves a graphical representation of the De Bruijn graph.

## Example
Here’s an example command to run the script:
```
python3 debruijn.py -i ../data/eva71_hundred_reads.fq -k 100 -o ../results/hundred_k100.fq -f ../results/hundred_k100.png
```
## Contact
For any inquiries, feel free to contact the author via email at takwa.ben-radhia@etu.u-paris.fr
