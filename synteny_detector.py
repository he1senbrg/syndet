import logging
from typing import List
from suffix_tree import SuffixTree
from models import SyntenyBlock
from config import DEFAULT_MAX_OVERLAP_RATIO

def find_common_substrings(seq1: str, seq2: str, min_length: int = 100) -> List[SyntenyBlock]:
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

def detect_synteny_blocks(seq1: str, seq2: str, min_length: int = 100) -> List[SyntenyBlock]:
    logging.info(f"Finding synteny blocks with minimum length {min_length}...")
    
    common_substrings = find_common_substrings(seq1, seq2, min_length)
    
    common_substrings.sort(key=lambda x: x.length, reverse=True)
    
    logging.info(f"Found {len(common_substrings)} common substrings with length >= {min_length}")
    
    max_overlap_allowed = int(min_length * DEFAULT_MAX_OVERLAP_RATIO)
    
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
