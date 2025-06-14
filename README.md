# Syndet - Synteny Block Detector

A Python tool for detecting synteny blocks (conserved genomic regions) between DNA sequences using Ukkonen's suffix tree algorithm for efficient comparative genomics analysis.

## Overview

Syndet implements a memory-optimized suffix tree construction algorithm to identify common substrings between genomic sequences. The tool processes FASTA-formatted sequences and produces detailed reports of conserved regions with configurable minimum length thresholds.

### Key Features

- **Efficient Algorithm**: Uses Ukkonen's linear-time suffix tree construction
- **Memory Optimization**: Employs `bytearray` and `array.array` for reduced memory footprint
- **Block Merging**: Automatically merges overlapping synteny blocks
- **Comprehensive Output**: Generates CSV reports with GC content analysis

## Installation

### Prerequisites
- Python 3.7 or higher

### Setup
```bash
git clone https://github.com/he1senbrg/syndet.git
cd syndet
```

## Usage

### Basic Command
```bash
python syndet.py <genome1.fasta> <genome2.fasta> <output.csv> [min_block_length]
```

### Parameters
- `genome1.fasta`: First FASTA sequence file
- `genome2.fasta`: Second FASTA sequence file  
- `output.csv`: Output file for synteny blocks
- `min_block_length`: Minimum synteny block length (default: 100)

### Example
```bash
python syndet.py datasets/ecoli_k12.fasta datasets/ecoli_o157.fasta results.csv 200
```

## Algorithm Details

The core algorithm centers on the `SuffixTree` class which implements several optimizations:

- **Generalized Suffix Tree**: Builds a combined tree for both input sequences
- **Source Tracking**: Uses bitmasks to track sequence origins
- **Overlap Filtering**: Removes excessive overlaps while preserving longest blocks

## Output Format

The tool generates CSV files containing:
- Block coordinates in both sequences
- Block length and identity percentage
- GC content analysis
- Sequence previews and summary statistics

## Data Structures

Synteny blocks are represented using a named tuple structure that encapsulates positional information, block characteristics, and quality metrics.

## Test Data

The repository includes sample datasets in the `datasets/` directory for testing with various genomic sequences including E. coli variants and viral genomes.

## Logging

The tool provides comprehensive logging to both console and file (`synteny.log`) for monitoring analysis progress and debugging.
