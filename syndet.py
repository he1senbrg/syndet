#!/usr/bin/env python3
"""
Synteny Block Detector
This script detects all synteny blocks between two genome sequences provided as FASTA files.
It uses Ukkonen's algorithm for efficient suffix tree construction to identify common substrings.
Usage: python syndet.py <genome1.fasta> <genome2.fasta> <output.txt> [min_block_length]
"""

import sys
import time
from collections import namedtuple
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

# For storing comprehensive information about a synteny block
SyntenyBlock = namedtuple("SyntenyBlock", [
    "seq1_start", "seq1_end", "seq2_start", "seq2_end", 
    "length", "sequence", "identity"
])

# Initialize constants
INF = float('inf')  # Definition of infinity

# ----------------- Ukkonen's Suffix Tree Implementation -----------------

class Node:
    """Node in the suffix tree."""
    def __init__(self, leaf=False):
        self.children = {}  # {char: (node, start_idx, end_idx)}
        self.leaf = leaf
        self.suffix_link = None
        self.id = -1  # Will be used to identify which string this suffix belongs to
        self.suffix_idx = -1  # Index of the suffix in the original string

class SuffixTree:
    """Suffix tree implementation using Ukkonen's algorithm."""
    def __init__(self):
        self.root = Node()
        self.root.suffix_link = self.root
        
    def build(self, text1, text2):
        """Build a generalized suffix tree for two strings using Ukkonen's algorithm."""
        # Combine sequences with a unique separator
        self.separator = "#"  # Must not appear in either string
        self.full_text = text1 + self.separator + text2
        self.text1_len = len(text1)
        self.text2_len = len(text2)
        self.length = len(self.full_text)
        
        # Active point tracking
        self.active_node = self.root
        self.active_edge = -1
        self.active_length = 0
        
        # Counter for remaining suffixes to be inserted
        self.remainder = 0
        self.pos = -1
        
        # For end position of leafs
        self.global_end = [-1]
        
        for i in range(self.length):
            self.pos = i
            self.global_end[0] = i
            self.remainder += 1
            
            # Insert all remaining suffixes
            last_created_node = None
            
            while self.remainder > 0:
                # If active length is zero, active edge must be updated
                if self.active_length == 0:
                    self.active_edge = i
                
                # Check if active edge exists from active node
                if self.active_edge >= 0 and self.full_text[self.active_edge] not in self.active_node.children:
                    # Create a new leaf node
                    leaf = Node(leaf=True)
                    if i < self.text1_len:
                        leaf.id = 1
                        leaf.suffix_idx = i - self.remainder + 1
                    else:
                        leaf.id = 2
                        leaf.suffix_idx = i - self.remainder + 1 - self.text1_len - 1
                    
                    self.active_node.children[self.full_text[self.active_edge]] = (leaf, i, self.global_end)
                    
                    # Set suffix link for previously created internal node
                    if last_created_node is not None:
                        last_created_node.suffix_link = self.active_node
                        last_created_node = None
                    
                else:
                    # Next character in the tree
                    next_node, edge_start, edge_end = self.active_node.children[self.full_text[self.active_edge]]
                    edge_length = self._edge_length(edge_start, edge_end)
                    
                    # If active length is within the current edge
                    if self.active_length >= edge_length:
                        # Move active point down
                        self.active_edge += edge_length
                        self.active_length -= edge_length
                        self.active_node = next_node
                        continue
                    
                    # If next character already exists in the tree
                    if self.full_text[edge_start + self.active_length] == self.full_text[i]:
                        self.active_length += 1
                        
                        # Set suffix link for previously created internal node
                        if last_created_node is not None:
                            last_created_node.suffix_link = self.active_node
                            last_created_node = None
                        
                        # Rule 3: Done for this phase
                        break
                    
                    # Split the edge
                    split_node = Node()
                    self.active_node.children[self.full_text[self.active_edge]] = (split_node, edge_start, [edge_start + self.active_length - 1])
                    
                    # Create new leaf
                    leaf = Node(leaf=True)
                    if i < self.text1_len:
                        leaf.id = 1
                        leaf.suffix_idx = i - self.remainder + 1
                    else:
                        leaf.id = 2
                        leaf.suffix_idx = i - self.remainder + 1 - self.text1_len - 1
                    
                    split_node.children[self.full_text[i]] = (leaf, i, self.global_end)
                    
                    # Update original node in split node's children
                    split_node.children[self.full_text[edge_start + self.active_length]] = (next_node, edge_start + self.active_length, edge_end)
                    
                    # Set suffix link
                    if last_created_node is not None:
                        last_created_node.suffix_link = split_node
                    
                    last_created_node = split_node
                
                # Decrease remaining suffixes and follow suffix link if needed
                self.remainder -= 1
                
                if self.active_node == self.root and self.active_length > 0:
                    self.active_length -= 1
                    self.active_edge = i - self.remainder + 1
                else:
                    self.active_node = self.active_node.suffix_link if self.active_node.suffix_link else self.root
    
    def _edge_length(self, start, end):
        """Calculate edge length."""
        if isinstance(end, list):
            return end[0] - start + 1
        return end - start + 1
    
    def find_common_substrings(self, min_length=100):
        """Find all common substrings between the two texts."""
        common_substrings = []
        
        # DFS to find nodes with both strings
        self._dfs_find_common(self.root, 0, "", {1: False, 2: False}, common_substrings, min_length)
        
        return common_substrings
    
    def _dfs_find_common(self, node, path_len, current_path, string_ids, common_substrings, min_length):
        """DFS to find nodes representing common substrings."""
        # Mark which strings are represented in the current path
        if node.leaf:
            string_ids[node.id] = True
            return {1: node.id == 1, 2: node.id == 2}
        
        # Process each child
        strings_below = {1: False, 2: False}
        
        for char, (child, start, end) in node.children.items():
            # Skip the separator character
            if self.full_text[start] == self.separator:
                continue
            
            # Calculate edge length
            edge_len = self._edge_length(start, end)
            child_path = self.full_text[start:start + edge_len] if isinstance(end, list) else self.full_text[start:end + 1]
            
            # Recursively process the child
            child_strings = self._dfs_find_common(child, path_len + edge_len, current_path + child_path, 
                                                dict(string_ids), common_substrings, min_length)
            
            # Update strings below this node
            strings_below[1] = strings_below[1] or child_strings[1]
            strings_below[2] = strings_below[2] or child_strings[2]
        
        # If this path represents both strings and is long enough, add to common substrings
        if strings_below[1] and strings_below[2] and path_len >= min_length:
            # Find occurrences in both strings
            self._find_occurrences(current_path, common_substrings)
        
        return strings_below
    
    def _find_occurrences(self, pattern, common_substrings):
        """Find all occurrences of pattern in both strings and add to common substrings."""
        # Naive string matching is used here for simplicity
        # In a real implementation, we would use the suffix tree itself to locate the positions
        
        text1 = self.full_text[:self.text1_len]
        text2 = self.full_text[self.text1_len + 1:self.text1_len + 1 + self.text2_len]
        
        for i in range(len(text1) - len(pattern) + 1):
            if text1[i:i + len(pattern)] == pattern:
                for j in range(len(text2) - len(pattern) + 1):
                    if text2[j:j + len(pattern)] == pattern:
                        # Found a match in both texts
                        seq1_start = i
                        seq1_end = i + len(pattern) - 1
                        seq2_start = j
                        seq2_end = j + len(pattern) - 1
                        
                        common_substrings.append(SyntenyBlock(
                            seq1_start, seq1_end,
                            seq2_start, seq2_end,
                            len(pattern), pattern, 100.0  # 100% identity for exact matches
                        ))

def find_common_substrings(seq1, seq2, min_length=100):
    """Find common substrings between two sequences using Ukkonen's algorithm."""
    logging.info(f"Finding common substrings with minimum length {min_length} using Ukkonen's algorithm...")
    
    # Build suffix tree
    suffix_tree = SuffixTree()
    suffix_tree.build(seq1, seq2)
    
    # Find common substrings
    common_substrings = suffix_tree.find_common_substrings(min_length)
    
    logging.info(f"Found {len(common_substrings)} initial common substrings")
    
    # Merge overlapping blocks
    common_substrings.sort(key=lambda x: (x.seq1_start, x.seq2_start))
    merged_blocks = []
    
    if common_substrings:
        current_block = common_substrings[0]
        
        for next_block in common_substrings[1:]:
            # Check if blocks are adjacent or overlapping
            if (next_block.seq1_start <= current_block.seq1_end + 1 and 
                next_block.seq2_start <= current_block.seq2_end + 1):
                # Merge blocks
                current_block = SyntenyBlock(
                    current_block.seq1_start,
                    max(current_block.seq1_end, next_block.seq1_end),
                    current_block.seq2_start,
                    max(current_block.seq2_end, next_block.seq2_end),
                    max(current_block.seq1_end, next_block.seq1_end) - current_block.seq1_start + 1,
                    seq1[current_block.seq1_start:max(current_block.seq1_end, next_block.seq1_end) + 1],
                    100.0
                )
            else:
                merged_blocks.append(current_block)
                current_block = next_block
        
        merged_blocks.append(current_block)
    
    logging.info(f"After merging: {len(merged_blocks)} synteny blocks")
    return merged_blocks

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

def calculate_gc_content(sequence):
    """Calculate the GC content of a DNA sequence."""
    if not sequence:
        return 0.0
    gc_count = sequence.count('G') + sequence.count('g') + sequence.count('C') + sequence.count('c')
    return (gc_count / len(sequence)) * 100.0

def detect_synteny_blocks(seq1, seq2, min_length=100):
    """Detect all synteny blocks between two sequences using Ukkonen's algorithm."""
    logging.info(f"Finding synteny blocks with minimum length {min_length}...")
    
    # Find common substrings using Ukkonen's algorithm
    common_substrings = find_common_substrings(seq1, seq2, min_length)
    
    # Sort by length (descending)
    common_substrings.sort(key=lambda x: x.length, reverse=True)
    
    logging.info(f"Found {len(common_substrings)} common substrings with length >= {min_length}")
    
    # Filter for non-overlapping blocks
    selected_blocks = []
    max_overlap_allowed = min_length // 4  # Allow some overlap (25% of min length)
    
    count = 1
    for block in common_substrings:
        logging.info(f"Processing block {count} with length {len(block)} of {len(common_substrings)} ")
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
        logging.info("Usage: python syndet.py <genome1.fasta> <genome2.fasta> <output.txt> [min_block_length]")
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
    
    write_synteny_blocks(blocks, output_file, header1, header2, end_time - start_time)
    logging.info(f"Synteny blocks written to {output_file}")

if __name__ == "__main__":
    main()