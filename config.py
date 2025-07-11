import logging

DEFAULT_MIN_BLOCK_LENGTH = 100
DEFAULT_MAX_OVERLAP_RATIO = 0.25

def setup_logging():
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[
            logging.FileHandler("synteny.log", mode='w'),
            logging.StreamHandler()
        ]
    )
