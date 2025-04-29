#!/usr/bin/env python3
"""
Synteny Block Detector
This script detects all synteny blocks between two genome sequences provided as FASTA files.
It uses Ukkonen's algorithm for efficient suffix tree construction to identify common substrings.
Usage: python syndet.py <genome1.fasta> <genome2.fasta> <output.txt> [min_block_length]
"""

import sys
import time
from collections import namedtuple, defaultdict
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

# ----------------- Optimized Ukkonen's Suffix Tree Implementation -----------------

class Edge:
    """Edge in the suffix tree."""
    def __init__(self, first_char, start_node, end_node, start_index, end_index):
        self.first_char = first_char
        self.start_node = start_node
        self.end_node = end_node
        self.start_index = start_index
        self.end_index = end_index
        
    def length(self, current_index):
        """Get the length of the edge."""
        return min(self.end_index, current_index + 1) - self.start_index

class Node:
    """Node in the suffix tree."""
    def __init__(self, suffix_index=-1, string_id=None):
        self.edges = {}  # {first_char: Edge}
        self.suffix_link = None
        self.suffix_index = suffix_index  # -1 for internal nodes
        self.string_id = string_id  # 1 for first string, 2 for second string
        self.end_indices = []  # List of end indices for this node
        self.sources = set()  # Which strings this node represents (1, 2, or both)

class SuffixTree:
    """Optimized suffix tree implementation using Ukkonen's algorithm."""
    def __init__(self):
        self.root = Node()
        self.active_node = self.root
        self.active_edge = -1
        self.active_length = 0
        self.remaining = 0
        self.current_index = -1
        self.nodes = [self.root]  # List of all nodes for easy traversal
        
    def edge_lookup(self, node, first_char):
        """Look up an edge by first character."""
        if first_char in node.edges:
            return node.edges[first_char]
        return None
    
    def walk_down(self, current_node, edge):
        """Walk down the tree to find the right active point."""
        if self.active_length >= edge.length(self.current_index):
            self.active_edge += edge.length(self.current_index)
            self.active_length -= edge.length(self.current_index)
            self.active_node = edge.end_node
            return True
        return False
    
    def build_suffix_tree(self, text, string_id=1):
        """Build a suffix tree for a single text."""
        # Save original string and its length
        start_index = len(getattr(self, 'text', ''))
        if not hasattr(self, 'text'):
            self.text = text
        else:
            self.text += text
        
        # Process one character at a time
        for i in range(len(text)):
            self.current_index = start_index + i
            
            # Add all remaining suffixes
            self.remaining += 1
            last_new_node = None
            
            while self.remaining > 0:
                # If active length is zero, active edge needs to be updated
                if self.active_length == 0:
                    self.active_edge = self.current_index
                
                # Get the active edge
                active_edge_char = self.text[self.active_edge]
                edge = self.edge_lookup(self.active_node, active_edge_char)
                
                # Create a new edge if needed
                if edge is None:
                    # Create new leaf node
                    leaf = Node(self.current_index - self.remaining + 1, string_id)
                    leaf.sources.add(string_id)
                    self.nodes.append(leaf)
                    
                    # Create new edge
                    new_edge = Edge(active_edge_char, self.active_node, leaf, self.active_edge, sys.maxsize)
                    self.active_node.edges[active_edge_char] = new_edge
                    
                    # Set suffix link for the last created node
                    if last_new_node is not None:
                        last_new_node.suffix_link = self.active_node
                        last_new_node = None
                else:
                    # Walk down to find the right active point
                    if self.walk_down(self.active_node, edge):
                        continue
                    
                    # If current character is already on the edge
                    next_char = self.text[edge.start_index + self.active_length]
                    if next_char == self.text[self.current_index]:
                        # Rule 3: already in tree, we're done for this phase
                        self.active_length += 1
                        
                        # Set suffix link for the last created node
                        if last_new_node is not None:
                            last_new_node.suffix_link = self.active_node
                            last_new_node = None
                        break
                    
                    # Create a new internal node
                    split_node = Node()
                    split_node.sources.add(string_id)
                    self.nodes.append(split_node)
                    
                    # Create a new edge for the split node
                    original_edge = edge
                    split_edge = Edge(active_edge_char, self.active_node, split_node, 
                                     edge.start_index, edge.start_index + self.active_length)
                    self.active_node.edges[active_edge_char] = split_edge
                    
                    # Create a new leaf node
                    leaf = Node(self.current_index - self.remaining + 1, string_id)
                    leaf.sources.add(string_id)
                    self.nodes.append(leaf)
                    
                    # Create edges from split node
                    split_node.edges[self.text[self.current_index]] = Edge(
                        self.text[self.current_index], split_node, leaf, self.current_index, sys.maxsize)
                    
                    # Update original edge
                    original_edge.start_index += self.active_length
                    original_edge.start_node = split_node
                    split_node.edges[self.text[original_edge.start_index]] = original_edge
                    
                    # Set suffix link for the last created node
                    if last_new_node is not None:
                        last_new_node.suffix_link = split_node
                    last_new_node = split_node
                
                # Rule 1: decrease remaining and follow suffix link
                self.remaining -= 1
                if self.active_node == self.root and self.active_length > 0:
                    self.active_length -= 1
                    self.active_edge = self.current_index - self.remaining + 1
                else:
                    self.active_node = self.active_node.suffix_link if self.active_node.suffix_link else self.root
        
        # Set correct end indices for leaves
        for node in self.nodes:
            for edge_char, edge in node.edges.items():
                if edge.end_index == sys.maxsize:
                    edge.end_index = self.current_index + 1
        
        # Mark the source for each node
        self._mark_sources(string_id)
        
        return self
    
    def _mark_sources(self, string_id):
        """Mark which string each node represents."""
        # First mark all leaf nodes
        for node in self.nodes:
            if node.string_id == string_id:
                node.sources.add(string_id)
        
        # Propagate up the tree
        self._propagate_sources(self.root)
    
    def _propagate_sources(self, node):
        """Propagate source information up the tree."""
        if not node.edges:  # Leaf node
            return node.sources
        
        # For internal nodes, collect sources from all children
        for edge_char, edge in node.edges.items():
            sources = self._propagate_sources(edge.end_node)
            node.sources.update(sources)
        
        return node.sources
    
    def build_generalized_suffix_tree(self, text1, text2):
        """Build a generalized suffix tree for two texts."""
        # Save text lengths
        self.text1_len = len(text1)
        self.text2_len = len(text2)
        
        # Build tree for first text
        self.build_suffix_tree(text1, 1)
        
        # Add a separator
        self.separator = "#"
        self.build_suffix_tree(self.separator, None)
        
        # Continue with second text
        self.build_suffix_tree(text2, 2)
        
        return self
    
    def _extract_string(self, start_idx, end_idx):
        """Extract a substring from the text."""
        if start_idx < 0 or end_idx > len(self.text):
            return ""
        return self.text[start_idx:end_idx]
    
    def _collect_starting_positions(self, node, positions1, positions2):
        """Collect starting positions of suffixes in both texts."""
        if node.suffix_index != -1:
            if node.string_id == 1:
                positions1.append(node.suffix_index)
            elif node.string_id == 2:
                positions2.append(node.suffix_index - self.text1_len - 1)  # Adjust for separator
            return
        
        # Recursive traversal
        for edge_char, edge in node.edges.items():
            self._collect_starting_positions(edge.end_node, positions1, positions2)
    
    def find_common_substrings(self, min_length=100):
        """Find all common substrings that occur in both texts."""
        common_substrings = []
        
        # Use DFS to find nodes that represent both strings
        self._find_common_substrings_dfs(self.root, "", common_substrings, min_length)
        
        return common_substrings
    
    def _find_common_substrings_dfs(self, node, current_string, common_substrings, min_length):
        """DFS to find common substrings."""
        # Skip nodes that don't represent both strings
        if 1 not in node.sources or 2 not in node.sources:
            return
        
        # For each edge from this node
        for edge_char, edge in node.edges.items():
            # Skip separator
            if edge_char == self.separator:
                continue
            
            # Get substring for this edge
            substring = self._extract_string(edge.start_index, edge.end_index)
            new_string = current_string + substring
            
            # Check if this node represents both strings
            if (1 in edge.end_node.sources and 2 in edge.end_node.sources and 
                len(new_string) >= min_length):
                
                # Collect starting positions
                positions1 = []
                positions2 = []
                self._collect_starting_positions(edge.end_node, positions1, positions2)
                
                # Create synteny blocks for each pair of positions
                for pos1 in positions1:
                    for pos2 in positions2:
                        seq1_start = pos1
                        seq1_end = pos1 + len(new_string) - 1
                        seq2_start = pos2
                        seq2_end = pos2 + len(new_string) - 1
                        
                        if seq1_end < self.text1_len and seq2_end < self.text2_len:
                            common_substrings.append(SyntenyBlock(
                                seq1_start, seq1_end,
                                seq2_start, seq2_end,
                                len(new_string), 
                                self.text[pos1:pos1 + len(new_string)],
                                100.0  # 100% identity for exact matches
                            ))
            
            # Continue DFS
            self._find_common_substrings_dfs(edge.end_node, new_string, common_substrings, min_length)

def find_common_substrings(seq1, seq2, min_length=100):
    """Find common substrings between two sequences using Ukkonen's algorithm."""
    logging.info(f"Finding common substrings with minimum length {min_length} using Ukkonen's algorithm...")
    
    # Build suffix tree
    suffix_tree = SuffixTree()
    suffix_tree.build_generalized_suffix_tree(seq1, seq2)
    
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
    
    for block in common_substrings:
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