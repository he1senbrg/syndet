# Synteny Block Detector

A Python tool for detecting synteny blocks (conserved regions) between genomic sequences.

## Overview

This script uses a k-mer based approach for detecting common substrings between two sequences. Here's how it works:

- K-mer indexing: It extracts all k-length substrings (k-mers) from the first sequence and stores their positions
- Matching: It scans the second sequence for the same k-mers.
- Extension: When a matching k-mer is found, it tries to extend the match left and right as far as possible while the characters match.
- Merging: Overlapping or adjacent matches are merged into longer synteny blocks.
- Filtering: It applies heuristics to filter overlapping blocks and keep the best ones.

This is a heuristic method, more memory-efficient and suitable for large genomic sequences, but it does not construct a suffix tree or suffix array, and thus avoids the complexity of Ukkonen’s or McCreight’s algorithms.

## Installation

### Prerequisites

- Python 3.7 or higher

### Setup

1. Clone this repository or download the script:

    ```bash
    git clone https://github.com/he1senbrg/syndet.git
    cd syndet
    ```

## Usage

### Basic Usage

```bash
python syndet.py <genome1.fasta> <genome2.fasta> <output.txt> [min_block_length]
```

Where:
- `<genome1.fasta>`: Path to the first FASTA file
- `<genome2.fasta>`: Path to the second FASTA file
- `[min_block_length]`: Parameter to specify the minimum length of synteny blocks.

### Example Usage

```bash
python syndet.py datasets/ecoli_k12.fasta datasets/ecoli_o157.fasta out_100.txt 200
```
