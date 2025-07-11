import array
from typing import List, Tuple, Optional, Union
from models import Node, SyntenyBlock

class SuffixTree:
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
        self.text1_len = 0
        self.text2_len = 0
        self.separator = ord('#')
        self.separator_idx = 0
    
    def edge_lookup(self, node: Node, first_char: int) -> Optional[Tuple[int, int, Node]]:
        if first_char in node.edges:
            return node.edges[first_char]
        return None
    
    def get_edge_length(self, start_idx: int, end_idx_ref: int) -> int:
        if end_idx_ref == self.infinity:
            return self.current_index + 1 - start_idx
        return self.end_indices[end_idx_ref] - start_idx
    
    def walk_down(self, node: Node, start_idx: int, end_idx_ref: int, target_node: Node) -> bool:
        edge_length = self.get_edge_length(start_idx, end_idx_ref)
        
        if self.active_length >= edge_length:
            self.active_edge += edge_length
            self.active_length -= edge_length
            self.active_node = target_node
            return True
        return False
    
    def build_suffix_tree(self, text_str: Union[str, bytearray], string_id: Optional[int] = 1) -> 'SuffixTree':
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
    
    def _propagate_sources(self) -> None:
        self._propagate_sources_helper(self.root)
    
    def _propagate_sources_helper(self, node: Node) -> int:
        if not node.edges:
            return node.sources
        
        sources = 0
        for edge_char, (start_idx, end_idx_ref, target_node) in node.edges.items():
            child_sources = self._propagate_sources_helper(target_node)
            sources |= child_sources
        
        node.sources |= sources
        return node.sources
    
    def build_generalized_suffix_tree(self, text1: str, text2: str) -> 'SuffixTree':
        self.text1_len = len(text1)
        self.text2_len = len(text2)
        
        self.build_suffix_tree(text1, 1)
        
        self.separator_idx = len(self.text)
        self.build_suffix_tree(bytearray([self.separator]), None)
        
        self.build_suffix_tree(text2, 2)
        
        return self
    
    def find_common_substrings(self, min_length: int = 100) -> List[SyntenyBlock]:
        common_substrings = []
        
        self._find_common_substrings_dfs(
            self.root, 
            0,
            common_substrings, 
            min_length,
            []
        )
        
        return common_substrings
    
    def _find_common_substrings_dfs(self, node: Node, path_length: int, 
                                   common_substrings: List[SyntenyBlock], 
                                   min_length: int, path_stack: List[Tuple[int, int]]) -> None:
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
    
    def _collect_positions(self, node: Node, seq1_positions: List[int], seq2_positions: List[int]) -> None:
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
