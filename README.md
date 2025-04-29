# Synteny Block Detector

A Python tool for detecting synteny blocks (conserved regions) between genomic sequences using Ukkonen's suffix tree algorithm.

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