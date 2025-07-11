from collections import namedtuple

SyntenyBlock = namedtuple("SyntenyBlock", [
    "seq1_start", "seq1_end", "seq2_start", "seq2_end", 
    "length", "sequence", "identity"
])

class Node:
    __slots__ = ['edges', 'suffix_link', 'suffix_index', 'string_id', 'sources']
    
    def __init__(self, suffix_index=-1, string_id=None):
        self.edges = {}
        self.suffix_link = None
        self.suffix_index = suffix_index
        self.string_id = string_id
        self.sources = 0
