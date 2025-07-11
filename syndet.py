import sys
import time
import logging
from config import setup_logging, DEFAULT_MIN_BLOCK_LENGTH
from utils import read_fasta, write_synteny_blocks
from synteny_detector import detect_synteny_blocks

def main():
    setup_logging()
    
    if len(sys.argv) < 4:
        print("Usage: python syndet.py <genome1.fasta> <genome2.fasta> <output> [min_block_length]")
        print("\tFile names should not include the .csv extension (or any other), it will be added automatically.")
        print("\n\tExample: python syndet.py genome1.fasta genome2.fasta file_name 100")

        sys.exit(1)
    
    genome1_file = sys.argv[1]
    genome2_file = sys.argv[2]
    output_file = sys.argv[3] + ".csv"
    min_block_length = int(sys.argv[4]) if len(sys.argv) > 4 else DEFAULT_MIN_BLOCK_LENGTH
    
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
    
if __name__ == "__main__":
    main()