# Syndet - Synteny Block Detector

A Python tool for detecting synteny blocks (conserved genomic regions) between DNA sequences using Ukkonen's suffix tree algorithm for efficient comparative genomics analysis.

## Overview

Syndet implements a memory-optimized suffix tree construction algorithm to identify common substrings between genomic sequences. The tool processes FASTA-formatted sequences and produces detailed reports of conserved regions with configurable minimum length thresholds.

### Key Features

- **Efficient Algorithm**: Uses Ukkonen's linear-time suffix tree construction
- **Memory Optimization**: Employs `bytearray` and `array.array` for reduced memory footprint
- **Block Merging**: Automatically merges overlapping synteny blocks
- **Comprehensive Output**: Generates CSV reports with GC content analysis
- **Modular Design**: Clean separation of concerns with multiple modules

## Project Structure

```
syndet/
├── syndet.py              # Main CLI entry point
├── models.py              # Data structures and named tuples
├── suffix_tree.py         # Core suffix tree implementation
├── synteny_detector.py    # Synteny detection business logic
├── utils.py               # File I/O and utility functions
├── config.py              # Configuration and constants
├── README.md              # This file
├── .gitignore            # Git ignore patterns
└── datasets/             # Sample data files
```

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
python -m syndet <genome1.fasta> <genome2.fasta> <output.csv> [min_block_length]
```

### Parameters
- `genome1.fasta`: First FASTA sequence file
- `genome2.fasta`: Second FASTA sequence file  
- `output.csv`: Output file for synteny blocks
- `min_block_length`: Minimum synteny block length (default: 100)

### Example
```bash
python -m syndet datasets/ecoli_k12.fasta datasets/ecoli_o157.fasta results.csv 200
```

## Module Overview

### Core Components

- **`models.py`**: Contains `SyntenyBlock` named tuple and `Node` class for suffix tree nodes
- **`suffix_tree.py`**: Implements the core Ukkonen's suffix tree algorithm with memory optimizations
- **`synteny_detector.py`**: High-level synteny detection logic including block merging and filtering
- **`utils.py`**: File I/O operations, FASTA parsing, and utility functions like GC content calculation
- **`config.py`**: Configuration constants and logging setup

### Algorithm Details

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

## Test Data

The repository includes sample datasets in the `datasets/` directory for testing with various genomic sequences including E. coli variants and viral genomes.

## Logging

The tool provides comprehensive logging to both console and file (`synteny.log`) for monitoring analysis progress and debugging.

## Development

To extend the functionality:

1. **Add new algorithms**: Implement in separate modules following the existing pattern
2. **Modify output formats**: Update functions in `utils.py`
3. **Add new data structures**: Define in `models.py`
4. **Adjust parameters**: Modify constants in `config.py`
