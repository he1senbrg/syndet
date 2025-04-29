#!/usr/bin/env python3
"""
Improved Synteny Block Detector
This script detects all synteny blocks between two genome sequences provided as FASTA files.
Usage: python syndet.py <genome1.fasta> <genome2.fasta> <output.txt> [min_block_length]
"""

import sys
import time
from collections import defaultdict, namedtuple
import logging

# Set up logging to both file and console
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[
        logging.FileHandler("synteny.log", mode='w'),
        logging.StreamHandler()  # This sends logs to the console
    ]
)

# Definition of infinity for Ukkonen's algorithm
INF = float('inf')

# For storing comprehensive information about a synteny block
SyntenyBlock = namedtuple("SyntenyBlock", [
    "seq1_start", "seq1_end", "seq2_start", "seq2_end", 
    "length", "sequence", "identity"
])

def read_fasta(file_path):
    """Read a FASTA file and return the sequence and header."""
    sequences = []
    headers = []
    current_seq = []
    current_header = ""
    
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if current_seq:
                    sequences.append(''.join(current_seq))
                    headers.append(current_header)
                    current_seq = []
                current_header = line[1:]  # Remove '>' character
            else:
                current_seq.append(line)
    
    if current_seq:
        sequences.append(''.join(current_seq))
        headers.append(current_header)
    
    if not sequences:
        return "", ""
    
    return sequences[0], headers[0]

def find_common_substrings(text1, text2, min_length=100):
    """Find common substrings using a sliding window approach."""
    common_substrings = []
    
    # Create a dictionary to store all k-mers from text1
    k = min_length  # Use the minimum length as k-mer size
    text1_kmers = {}
    
    # Store all k-mers from text1 with their positions
    for i in range(len(text1) - k + 1):
        kmer = text1[i:i+k]
        if kmer not in text1_kmers:
            text1_kmers[kmer] = []
        text1_kmers[kmer].append(i)
    
    # Look for matches in text2
    for i in range(len(text2) - k + 1):
        kmer = text2[i:i+k]
        if kmer in text1_kmers:
            # For each occurrence in text1
            for pos1 in text1_kmers[kmer]:
                # Extend match as far as possible in both directions
                left_extension = 0
                right_extension = k
                
                # Extend to the left
                while pos1 - left_extension - 1 >= 0 and i - left_extension - 1 >= 0 and \
                      text1[pos1 - left_extension - 1] == text2[i - left_extension - 1]:
                    left_extension += 1
                
                # Extend to the right
                while pos1 + right_extension < len(text1) and i + right_extension < len(text2) and \
                      text1[pos1 + right_extension] == text2[i + right_extension]:
                    right_extension += 1
                
                # Calculate final positions
                seq1_start = pos1 - left_extension
                seq1_end = pos1 + right_extension - 1
                seq2_start = i - left_extension
                seq2_end = i + right_extension - 1
                length = right_extension + left_extension
                
                # Extract the sequence
                sequence = text1[seq1_start:seq1_end + 1]
                
                # Add to results
                if length >= min_length:
                    common_substrings.append(SyntenyBlock(
                        seq1_start, seq1_end,
                        seq2_start, seq2_end,
                        length, sequence, 100.0  # 100% identity for exact matches
                    ))
    
    return common_substrings

def calculate_gc_content(sequence):
    """Calculate the GC content of a DNA sequence."""
    if not sequence:
        return 0.0
    gc_count = sequence.count('G') + sequence.count('g') + sequence.count('C') + sequence.count('c')
    return (gc_count / len(sequence)) * 100.0

def detect_synteny_blocks(seq1, seq2, min_length=100):
    """Detect all synteny blocks between two sequences."""
    logging.info(f"Finding all synteny blocks with minimum length {min_length}...")
    
    # Find all common substrings
    common_substrings = find_common_substrings(seq1, seq2, min_length)
    
    # Sort by length (descending)
    common_substrings.sort(key=lambda x: x.length, reverse=True)
    
    logging.info(f"Found {len(common_substrings)} common substrings with length >= {min_length}")
    
    # OPTION 1: Return all blocks without filtering
    # return common_substrings
    
    # OPTION 2: Filter overlapping blocks (greedy approach) but with a less restrictive filter
    selected_blocks = []
    max_overlap_allowed = min_length // 4  # Allow some overlap (25% of min length)
    
    count = 1
    for block in common_substrings:
        logging.info(f"Processing block {count} with length {block.length} of {len(common_substrings)}")
        count += 1
        # Skip if too much overlap with existing blocks
        should_add = True
        for existing in selected_blocks:
            # Check overlap in seq1
            s1_overlap_start = max(block.seq1_start, existing.seq1_start)
            s1_overlap_end = min(block.seq1_end, existing.seq1_end)
            s1_overlap = max(0, s1_overlap_end - s1_overlap_start + 1)
            
            # Check overlap in seq2
            s2_overlap_start = max(block.seq2_start, existing.seq2_start)
            s2_overlap_end = min(block.seq2_end, existing.seq2_end)
            s2_overlap = max(0, s2_overlap_end - s2_overlap_start + 1)
            
            # If overlap is too large in either sequence, skip this block
            if s1_overlap > max_overlap_allowed or s2_overlap > max_overlap_allowed:
                should_add = False
                break
        
        if should_add:
            selected_blocks.append(block)
    
    # Sort by position in seq1
    selected_blocks.sort(key=lambda x: x.seq1_start)
    
    logging.info(f"Selected {len(selected_blocks)} non-overlapping blocks (allowing {max_overlap_allowed} bp overlap)")
    
    return selected_blocks

def write_synteny_blocks(blocks, output_file, seq1_header, seq2_header, elapsed_time):
    """Write synteny blocks to an output file with detailed information."""
    with open(output_file, 'w') as f:
        f.write(f"# Synteny blocks between:\n")
        f.write(f"# Sequence 1: {seq1_header}\n")
        f.write(f"# Sequence 2: {seq2_header}\n")
        f.write(f"# Total blocks found: {len(blocks)}\n")
        f.write(f"# Analysis date: {time.strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        
        f.write("Block_ID,Seq1_Start,Seq1_End,Seq2_Start,Seq2_End,Length,Identity,GC_Content,Sequence\n")
        
        for i, block in enumerate(blocks):
            gc_content = calculate_gc_content(block.sequence)
            
            f.write(f"{i+1},{block.seq1_start},{block.seq1_end},{block.seq2_start},"
                    f"{block.seq2_end},{block.length},{block.identity:.2f},{gc_content:.2f},"
                    f"{block.sequence[:50]}{'...' if len(block.sequence) > 50 else ''}\n")
        
        # Summary statistics
        if blocks:
            avg_length = sum(b.length for b in blocks) / len(blocks)
            total_synteny = sum(b.length for b in blocks)
            f.write(f"\n# Summary Statistics:\n")
            f.write(f"# Average block length: {avg_length:.2f}\n")
            f.write(f"# Total nucleotides in synteny: {total_synteny}\n")
        
        f.write(f"# Total time taken: {elapsed_time:.2f} seconds\n")


def main():
    if len(sys.argv) < 4:
        logging.info("Usage: python improved_synteny_detector.py <genome1.fasta> <genome2.fasta> <output.txt> [min_block_length]")
        sys.exit(1)
    
    genome1_file = sys.argv[1]
    genome2_file = sys.argv[2]
    output_file = sys.argv[3]
    min_block_length = int(sys.argv[4]) if len(sys.argv) > 4 else 100
    
    logging.info(f"Reading sequences from {genome1_file} and {genome2_file}...")
    seq1, header1 = read_fasta(genome1_file)
    seq2, header2 = read_fasta(genome2_file)
    
    logging.info(f"Sequence 1 length: {len(seq1)}")
    logging.info(f"Sequence 2 length: {len(seq2)}")
    
    # Handle case where sequences are very short
    if len(seq1) < min_block_length or len(seq2) < min_block_length:
        logging.info(f"Warning: Sequences are shorter than the minimum block length ({min_block_length}).")
        logging.info("Adjusting minimum block length to the shorter sequence length.")
        min_block_length = min(len(seq1), len(seq2)) // 2  # Set to half the shorter sequence length
    
    start_time = time.time()
    blocks = detect_synteny_blocks(seq1, seq2, min_block_length)
    end_time = time.time()
    
    logging.info(f"Found {len(blocks)} synteny blocks in {end_time - start_time:.2f} seconds")
    
    write_synteny_blocks(blocks, output_file, header1, header2,end_time - start_time)
    logging.info(f"Synteny blocks written to {output_file}")
    
    # logging.info summary of all blocks to console
    logging.info("\nSummary of synteny blocks:")
    logging.info("Block ID\tSeq1 Range\t\tSeq2 Range\t\tLength")
    logging.info("-" * 70)
    for i, block in enumerate(blocks):
        logging.info(f"{i+1}\t\t{block.seq1_start}-{block.seq1_end}\t\t{block.seq2_start}-{block.seq2_end}\t\t{block.length}")

if __name__ == "__main__":
    main()