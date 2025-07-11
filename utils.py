import time
import os
import logging
from typing import Tuple, List
from models import SyntenyBlock

OUT_DIR = "output"
CSV_DIR = "csv"
INFO_DIR = "info"

os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(OUT_DIR + "/" + CSV_DIR, exist_ok=True)
os.makedirs(OUT_DIR + "/" + INFO_DIR, exist_ok=True)

def read_fasta(file_path: str) -> Tuple[str, str]:
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

def calculate_gc_content(sequence: str) -> float:
    if not sequence:
        return 0.0
    
    gc_count = 0
    for char in sequence:
        if char in 'GCgc':
            gc_count += 1
            
    return (gc_count / len(sequence)) * 100.0

def write_synteny_blocks(blocks: List[SyntenyBlock], output_file: str, seq1_header: str, seq2_header: str, elapsed_time: float) -> None:
    if blocks:
        with open(OUT_DIR + "/" + CSV_DIR + "/" + output_file, 'w') as f:
            f.write("Block_ID,Seq1_Start,Seq1_End,Seq2_Start,Seq2_End,Length,Identity,GC_Content,Sequence\n")
            
            for i, block in enumerate(blocks):
                gc_content = calculate_gc_content(block.sequence)
                
                f.write(f"{i+1},{block.seq1_start},{block.seq1_end},{block.seq2_start},"
                        f"{block.seq2_end},{block.length},{block.identity:.2f},{gc_content:.2f},"
                        f"{block.sequence}\n")

        with open(OUT_DIR + "/" + INFO_DIR + "/" + output_file.replace('.csv', '_info.md'), 'w') as f:
            f.write(f"# Synteny Blocks Information\n")
            f.write(f"## Between:\n")
            f.write(f"- **Sequence 1**: {seq1_header}\n")
            f.write(f"- **Sequence 2**: {seq2_header}\n")
            f.write(f"## Total blocks found: {len(blocks)}\n")
            f.write(f"## Analysis date: {time.strftime('%Y-%m-%d %H:%M:%S')}\n\n")
            
            avg_length = sum(b.length for b in blocks) / len(blocks)
            total_synteny = sum(b.length for b in blocks)

            f.write(f"### Summary Statistics:\n")
            f.write(f"- Average block length: {avg_length:.2f}\n")
            f.write(f"- Total nucleotides in synteny: {total_synteny}\n")
            f.write(f"- Total time taken: {elapsed_time:.2f} seconds\n")

            f.write("\n### Block Details Preview\n")
            f.write("| Block ID | Seq1 Start | Seq1 End | Seq2 Start | Seq2 End | Length | Identity (%) | GC Content (%) | Sequence |\n")
            f.write("|----------|------------|----------|------------|----------|--------|---------------|----------------|----------|\n")
            
            for i, block in enumerate(blocks):
                gc_content = calculate_gc_content(block.sequence)
                f.write(f"| {i+1} | {block.seq1_start} | {block.seq1_end} | {block.seq2_start} | "
                        f"{block.seq2_end} | {block.length} | {block.identity:.2f} | "
                        f"{gc_content:.2f} | {block.sequence[:50]}{'...' if len(block.sequence) > 50 else ''} |\n")
            
        logging.info(f"Synteny blocks written to {output_file}")
    else:
        logging.info("No synteny blocks found to write to the output file.")
