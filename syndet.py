import sys
import time
from collections import namedtuple
import logging
import array

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[
        logging.FileHandler("synteny.log", mode='w'),
        logging.StreamHandler()
    ]
)

SyntenyBlock = namedtuple("SyntenyBlock", [
    "seq1_start", "seq1_end", "seq2_start", "seq2_end", 
    "length", "sequence", "identity"
])

class Node:
    """Node in the suffix tree with optimized memory footprint."""
    __slots__ = ['edges', 'suffix_link', 'suffix_index', 'string_id', 'sources']
    
    def __init__(self, suffix_index=-1, string_id=None):
        self.edges = {}
        self.suffix_link = None
        self.suffix_index = suffix_index
        self.string_id = string_id
        self.sources = 0

class SuffixTree:
    """Memory-optimized suffix tree implementation using Ukkonen's algorithm."""
    def __init__(self):
        self.root = Node()
        self.active_node = self.root
        self.active_edge = -1
        self.active_length = 0
        self.remaining = 0
        self.current_index = -1
        self.text = bytearray()
        self.end_indices = array.array('i')
        self.infinity = -1
    
    def edge_lookup(self, node, first_char):
        """Look up an edge by first character."""
        if first_char in node.edges:
            return node.edges[first_char]
        return None
    
    def get_edge_length(self, start_idx, end_idx_ref):
        """Get the current length of an edge."""
        if end_idx_ref == self.infinity:
            return self.current_index + 1 - start_idx
        return self.end_indices[end_idx_ref] - start_idx
    
    def walk_down(self, node, start_idx, end_idx_ref, target_node):
        """Walk down the tree to find the right active point."""
        edge_length = self.get_edge_length(start_idx, end_idx_ref)
        
        if self.active_length >= edge_length:
            self.active_edge += edge_length
            self.active_length -= edge_length
            self.active_node = target_node
            return True
        return False
    
    def build_suffix_tree(self, text_str, string_id=1):
        """Build a suffix tree for a single text."""
        text = text_str.encode('utf-8') if isinstance(text_str, str) else text_str
        
        start_index = len(self.text)
        self.text.extend(text)
        
        for i in range(len(text)):
            self.current_index = start_index + i
            
            self.remaining += 1
            last_new_node = None
            
            while self.remaining > 0:
                if self.active_length == 0:
                    self.active_edge = self.current_index
                
                active_edge_char = self.text[self.active_edge]
                edge_info = self.edge_lookup(self.active_node, active_edge_char)
                
                if edge_info is None:
                    leaf = Node(self.current_index - self.remaining + 1, string_id)
                    leaf.sources |= (1 << (string_id - 1)) if string_id else 0
                    
                    self.active_node.edges[active_edge_char] = (
                        self.active_edge, self.infinity, leaf
                    )
                    
                    if last_new_node is not None:
                        last_new_node.suffix_link = self.active_node
                        last_new_node = None
                else:
                    start_idx, end_idx_ref, target_node = edge_info
                    
                    if self.walk_down(self.active_node, start_idx, end_idx_ref, target_node):
                        continue
                    
                    edge_length = self.get_edge_length(start_idx, end_idx_ref)
                    if self.active_length < edge_length:
                        next_char_idx = start_idx + self.active_length
                        if next_char_idx < len(self.text) and self.text[next_char_idx] == self.text[self.current_index]:
                            self.active_length += 1
                            
                            if last_new_node is not None:
                                last_new_node.suffix_link = self.active_node
                                last_new_node = None
                            break
                    
                    split_node = Node()
                    split_node.sources |= (1 << (string_id - 1)) if string_id else 0
                    
                    new_end_idx = start_idx + self.active_length
                    
                    new_end_idx_ref = len(self.end_indices)
                    self.end_indices.append(new_end_idx)
                    
                    self.active_node.edges[active_edge_char] = (start_idx, new_end_idx_ref, split_node)
                    
                    leaf = Node(self.current_index - self.remaining + 1, string_id)
                    leaf.sources |= (1 << (string_id - 1)) if string_id else 0
                    
                    split_node.edges[self.text[self.current_index]] = (
                        self.current_index, self.infinity, leaf
                    )
                    
                    split_node.edges[self.text[new_end_idx]] = (
                        new_end_idx, end_idx_ref, target_node
                    )
                    
                    if last_new_node is not None:
                        last_new_node.suffix_link = split_node
                    last_new_node = split_node
                
                self.remaining -= 1
                if self.active_node == self.root and self.active_length > 0:
                    self.active_length -= 1
                    self.active_edge = self.current_index - self.remaining + 1
                else:
                    self.active_node = self.active_node.suffix_link if self.active_node.suffix_link else self.root
        
        for i in range(len(self.end_indices)):
            if self.end_indices[i] == -1:
                self.end_indices[i] = len(self.text)
        
        self._propagate_sources()
        
        return self
    
    def _propagate_sources(self):
        """Propagate source information up the tree."""
        self._propagate_sources_helper(self.root)
    
    def _propagate_sources_helper(self, node):
        """Helper function for propagating source information."""
        if not node.edges:
            return node.sources
        
        sources = 0
        for edge_char, (start_idx, end_idx_ref, target_node) in node.edges.items():
            child_sources = self._propagate_sources_helper(target_node)
            sources |= child_sources
        
        node.sources |= sources
        return node.sources
    
    def build_generalized_suffix_tree(self, text1, text2):
        """Build a generalized suffix tree for two texts."""
        self.text1_len = len(text1)
        self.text2_len = len(text2)
        
        self.build_suffix_tree(text1, 1)
        
        self.separator = ord('#')
        self.separator_idx = len(self.text)
        self.build_suffix_tree(bytearray([self.separator]), None)
        
        self.build_suffix_tree(text2, 2)
        
        return self
    
    def _extract_string(self, start_idx, end_idx):
        """Extract a substring from the text."""
        if start_idx < 0 or end_idx > len(self.text):
            return bytearray()
        return self.text[start_idx:end_idx]
    
    def find_common_substrings(self, min_length=100):
        """Find all common substrings that occur in both texts using efficient DFS."""
        common_substrings = []
        
        self._find_common_substrings_dfs(
            self.root, 
            0,
            common_substrings, 
            min_length,
            []
        )
        
        return common_substrings
    
    def _find_common_substrings_dfs(self, node, path_length, common_substrings, min_length, path_stack):
        """DFS to find common substrings with path stacking instead of string concatenation."""
        if node.sources != 3:
            return
        
        for edge_char, (start_idx, end_idx_ref, target_node) in node.edges.items():
            if edge_char == self.separator:
                continue
            
            end_idx = self.end_indices[end_idx_ref] if end_idx_ref != self.infinity else len(self.text)
            edge_length = end_idx - start_idx
            new_path_length = path_length + edge_length
            
            path_stack.append((start_idx, end_idx))
            
            if target_node.sources == 3 and new_path_length >= min_length:
                seq1_positions = []
                seq2_positions = []
                
                if target_node.suffix_index != -1:
                    if target_node.string_id == 1:
                        seq1_positions.append(target_node.suffix_index)
                    elif target_node.string_id == 2:
                        seq2_positions.append(target_node.suffix_index)
                
                if not seq1_positions or not seq2_positions:
                    self._collect_positions(target_node, seq1_positions, seq2_positions)
                
                if seq1_positions and seq2_positions:
                    sequence_bytes = bytearray()
                    for s_idx, e_idx in path_stack:
                        sequence_bytes.extend(self.text[s_idx:e_idx])
                    
                    try:
                        sequence = sequence_bytes.decode('utf-8')
                    except UnicodeDecodeError:
                        sequence = repr(sequence_bytes)
                    
                    for seq1_start in seq1_positions:
                        seq1_end = seq1_start + new_path_length - 1
                        
                        for seq2_start in seq2_positions:
                            if seq2_start > self.separator_idx:
                                seq2_start -= (self.text1_len + 1)
                            seq2_end = seq2_start + new_path_length - 1
                            
                            if seq1_end < self.text1_len and seq2_end < self.text2_len:
                                common_substrings.append(SyntenyBlock(
                                    seq1_start, seq1_end,
                                    seq2_start, seq2_end,
                                    new_path_length, 
                                    sequence,
                                    100.0
                                ))
            
            self._find_common_substrings_dfs(target_node, new_path_length, common_substrings, min_length, path_stack)
            
            path_stack.pop()
    
    def _collect_positions(self, node, seq1_positions, seq2_positions):
        """Collect starting positions from both sequences."""
        if node.suffix_index != -1:
            if node.string_id == 1:
                seq1_positions.append(node.suffix_index)
            elif node.string_id == 2:
                seq2_positions.append(node.suffix_index)
            return
        
        for edge_char, (_, _, target_node) in node.edges.items():
            if edge_char == self.separator:
                continue
                
            sources = target_node.sources
            if sources & 1:
                if not seq1_positions:
                    self._collect_positions(target_node, seq1_positions, seq2_positions)
            if sources & 2:
                if not seq2_positions:
                    self._collect_positions(target_node, seq1_positions, seq2_positions)
            
            if seq1_positions and seq2_positions:
                break

def find_common_substrings(seq1, seq2, min_length=100):
    """Find common substrings between two sequences using optimized Ukkonen's algorithm."""
    logging.info(f"Finding common substrings with minimum length {min_length} using optimized Ukkonen's algorithm...")
    
    suffix_tree = SuffixTree()
    suffix_tree.build_generalized_suffix_tree(seq1, seq2)
    
    common_substrings = suffix_tree.find_common_substrings(min_length)
    
    logging.info(f"Found {len(common_substrings)} initial common substrings")
    
    if common_substrings:
        common_substrings.sort(key=lambda x: (x.seq1_start, x.seq2_start))
        merged_blocks = []
        
        current_block = common_substrings[0]
        
        for next_block in common_substrings[1:]:
            if (next_block.seq1_start <= current_block.seq1_end + 1 and 
                next_block.seq2_start <= current_block.seq2_end + 1):
                seq1_end = max(current_block.seq1_end, next_block.seq1_end)
                seq2_end = max(current_block.seq2_end, next_block.seq2_end)
                new_length = seq1_end - current_block.seq1_start + 1
                
                sequence = (next_block.sequence if next_block.length > current_block.length 
                           else current_block.sequence)
                
                current_block = SyntenyBlock(
                    current_block.seq1_start,
                    seq1_end,
                    current_block.seq2_start,
                    seq2_end,
                    new_length,
                    sequence,
                    100.0
                )
            else:
                merged_blocks.append(current_block)
                current_block = next_block
        
        merged_blocks.append(current_block)
    else:
        merged_blocks = []
    
    logging.info(f"After merging: {len(merged_blocks)} synteny blocks")
    return merged_blocks

def read_fasta(file_path):
    """Read a FASTA file efficiently and return the sequence and header."""
    sequences = []
    headers = []
    current_seq = bytearray()
    current_header = ""
    
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if current_seq:
                    sequences.append(current_seq.decode('utf-8'))
                    headers.append(current_header)
                    current_seq = bytearray()
                current_header = line[1:]
            else:
                current_seq.extend(line.encode('utf-8'))
    
    if current_seq:
        sequences.append(current_seq.decode('utf-8'))
        headers.append(current_header)
    
    if not sequences:
        return "", ""
    
    return sequences[0], headers[0]

def calculate_gc_content(sequence):
    """Calculate the GC content of a DNA sequence efficiently."""
    if not sequence:
        return 0.0
    
    gc_count = 0
    for char in sequence:
        if char in 'GCgc':
            gc_count += 1
            
    return (gc_count / len(sequence)) * 100.0

def detect_synteny_blocks(seq1, seq2, min_length=100):
    """Detect all synteny blocks using optimized algorithm."""
    logging.info(f"Finding synteny blocks with minimum length {min_length}...")
    
    common_substrings = find_common_substrings(seq1, seq2, min_length)
    
    common_substrings.sort(key=lambda x: x.length, reverse=True)
    
    logging.info(f"Found {len(common_substrings)} common substrings with length >= {min_length}")
    
    max_overlap_allowed = min_length // 4
    
    selected_blocks = []
    covered_seq1 = []
    covered_seq2 = []
    
    for block in common_substrings:
        overlap_too_large = False
        
        for start, end in covered_seq1:
            overlap_start = max(block.seq1_start, start)
            overlap_end = min(block.seq1_end, end)
            overlap = max(0, overlap_end - overlap_start + 1)
            
            if overlap > max_overlap_allowed:
                overlap_too_large = True
                break
        
        if not overlap_too_large:
            for start, end in covered_seq2:
                overlap_start = max(block.seq2_start, start)
                overlap_end = min(block.seq2_end, end)
                overlap = max(0, overlap_end - overlap_start + 1)
                
                if overlap > max_overlap_allowed:
                    overlap_too_large = True
                    break
        
        if not overlap_too_large:
            selected_blocks.append(block)
            covered_seq1.append((block.seq1_start, block.seq1_end))
            covered_seq2.append((block.seq2_start, block.seq2_end))
    
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
    
    if len(seq1) < min_block_length or len(seq2) < min_block_length:
        logging.info(f"Warning: Sequences are shorter than the minimum block length ({min_block_length}).")
        logging.info("Adjusting minimum block length to the shorter sequence length.")
        min_block_length = min(len(seq1), len(seq2)) // 2
    
    start_time = time.time()
    blocks = detect_synteny_blocks(seq1, seq2, min_block_length)
    end_time = time.time()
    
    logging.info(f"Found {len(blocks)} synteny blocks in {end_time - start_time:.2f} seconds")
    
    write_synteny_blocks(blocks, output_file, header1, header2, end_time - start_time)
    logging.info(f"Synteny blocks written to {output_file}")

if __name__ == "__main__":
    main()