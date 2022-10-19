from typing import List, Set
import bte as mat
import pickle
import sys
import graphviz as gv
import copy


# def get_mcp(src, dest, alpha=0.01):
#     if dest.id == "sink":
#         return 0
#     else:
#         return alpha * src.value - src.clade_size + dest.clade_size

class PartitionNode:
    def __init__(self, id, value, clade_size) -> None:
        self.id = id
        self.clade_size = clade_size
        self.value = value # NOTE: Using num mutes for now...

        # tuple w/ pointer to next node, and cost of that edge
        self.edges = {
            "sib": None,
            "up": None,
            "down": None
        }

        # List of children ids in MATree
        self.mat_children_ids = set()

        self.mcp = -1           # maximum cost of all paths from this node to sink node
        self.max_edge = None    # edge that goes to the child that maximizes

        # TODO: Store number of mutations too?

    def __str__(self) -> str:
        return f"[id:{self.id}, val{self.value}]"

    def children(self) -> Set:
        """ Returns the number of this partition nodes children
        """
        return [child[0] for child in self.edges.values() if child is not None]
    
    def is_sink(self) -> bool:
        return self.clade_size == -1

    def __eq__(self, __o: object) -> bool:
        return self.id == __o.id
    
    def __hash__(self) -> int:
        return self.id.__hash__()


class PartitionPath:
    """ Represents a path through the Partition Graph as just the node ids that compose the
    partition (i.e., children are not externally accessible)
    """
    def __init__(self, partition_graph):
        self.pg = partition_graph
        self.partition_list = []
        self.size = 0
        self.cost = 0
        self.size_history = []
        self.val_history = []

    def get_partition_ids(self) -> List[str]:
        """ Removes all subsumed nodes and returns the resulting list of ids
        """
        part_list = copy.copy(self.partition_list)
        
        prev = self.partition_list[0]
        for curr in self.partition_list[1:]:
            if curr.clade_size >= prev.clade_size and len(curr.mat_children_ids) > 0:   # up edge
                for child in curr.mat_children_ids:
                    part_list.remove(self.pg.id2partition_node[child])

        return [node.id for node in part_list[:-1]] # Remove sink node

    def __len__(self):
        return self.size-1

    def append(self, partition_node, val):
        self.partition_list.append(partition_node)
        size = 1 - len(partition_node.mat_children_ids)
        self.size += size
        self.cost += val

        self.size_history.append(size)
        self.val_history.append(val)

    def get_cost(self):
        return self.cost / self.__len__()

    def rewind(self, num_steps=1):
        """ Undo the last `num_steps` appends
        """
        self.size -= sum(self.size_history[-num_steps:])
        self.cost -= sum(self.val_history[-num_steps:])
        self.size_history = self.size_history[:-num_steps]
        self.val_history = self.val_history[:-num_steps]
        self.partition_list = self.partition_list[:-num_steps]

    def copy(self):
        pp = PartitionPath(self.pg)
        pp.partition_list = copy.copy(self.partition_list)
        pp.size = self.size
        pp.cost = self.cost
        pp.size_history = self.size_history
        pp.val_history = self.val_history
        return pp

    def contains(self, node):
        return node in self.partition_list

    def legal_up_edge(self, node):
        """ Returns true if taking the up edge is valid for the given path (i.e., all children
        of the node are contained in the path list so far)
        """
        return all([self.pg.id2partition_node[id] in self.partition_list
                     for id in node.mat_children_ids])

    def beneficial_up_edge(self, node, val):
        curr_cost = self.get_cost()
        self.append(node, val)
        new_cost = self.get_cost()
        self.rewind()
        return curr_cost <= new_cost
        



class PartitionGraph:
    """ Defines a graph representation of MAT that efficiently finds the best partition of the
    leaves subject to some criteria (e.g., restrictions on the number of roots). Each path through
    the graph represents a set of MAT nodes that partitions the MAT into clades.    
    """

    def __init__(self, mattree, max_clade_size=12000) -> None:
        self.max_clade_size = max_clade_size
        self.total_muts = 0
        self.num_edges = 0
        self.num_leaves = 1
        self.num_paths = 1
        self.sink = PartitionNode("sink", 0, -1)

        node_list = mattree.depth_first_expansion(reverse=True)

        id2clade = {}
        id2partition_node = {}

        # Fencepost with arbitrary leaf in MAT
        self.root = PartitionNode(node_list[0].id, len(node_list[0].mutations), 1)
        id2partition_node[node_list[0].id] = self.root
        id2clade[node_list[0].id] = 1
        prev = self.root

        to_next_leaf = set()

        for i, node in enumerate(node_list[1:]):
            self.total_muts += len(node.mutations)
            if node.is_leaf():
                clade = 1
                self.num_leaves += 1
            else:
                clade = 0
                for child in node.children:
                    clade += id2clade[child.id]
            
            id2clade[node.id] = clade

            if clade > max_clade_size:
                continue

            curr = PartitionNode(node.id, len(node.mutations), clade)
            curr.mat_children_ids.update([node.id for node in node.children])
            id2partition_node[node.id] = curr

            num_in_edges = 1
            if node.is_leaf() and len(to_next_leaf) > 0:    # Node to leaf edge
                for pnode in to_next_leaf:
                    assert pnode.edges["down"] is None
                    pnode.edges["down"] = (curr, curr.value)
                    self.num_edges += 1
                    num_in_edges += 1
                to_next_leaf = set()
            
            if prev.clade_size <= curr.clade_size and len(curr.mat_children_ids) != 0:   # Up edge
                assert prev.edges["up"] is None
                value = len(node.mutations)
                for child in node.children:
                    value -= id2partition_node[child.id].value
                prev.edges["up"] = (curr, value)    # Accumulate node costs through partition path
                to_next_leaf.add(prev)

            else:   # leaf to leaf edge
                assert prev.edges["sib"] is None
                prev.edges["sib"] = (curr, len(node.mutations))
            
            self.num_edges += 1
            self.num_paths *= num_in_edges
            prev = curr

        num_in_edges = 0
        # Add nodes to sink
        for pnode in to_next_leaf:
            assert pnode.edges["down"] is None
            pnode.edges["down"] = (self.sink, 0)
            self.num_edges += 1
            num_in_edges += 1
        
        # Don't forget to add the root!
        assert prev.edges["down"] is None
        prev.edges["down"] = (self.sink, 0)
        self.num_edges += 1
        
        self.num_paths *= num_in_edges
        self.id2partition_node = id2partition_node

        self.__backwards_annotate(mattree)
        
    def __backwards_annotate(self, mattree, alpha=0.1):
        """ Annotate each edge (u, v) in the Partition Graph with maximum achievable path value
        (i.e., weights - alpha*|partition|) from node v to the sink node.
        """

        nl = mattree.depth_first_expansion(reverse=False)
        for node in nl:
            partition_node = self.id2partition_node[node.id]

            max_val = -sys.maxsize()
            max_edge = None
            for dir, edge in partition_node.edges.items():
                if edge is None:
                    continue
                if edge[0].is_sink():
                    edge_val = 0
                else:
                    if dir == "down":
                        partition_size = edge[0].num_mat_children - 1
                    else:
                        partition_size = 1 - edge[0].num_mat_children
                    edge_val = edge[1] - alpha * partition_size

                if max_val < edge_val + edge[0].mcp:
                    max_val = edge_val + edge[0].mcp
                    max_edge = edge

            partition_node.mcp = max_val
            # partition_node.max_edge = max_edge

    def postorder(self):
        """Recursive postorder traversal of the Partion Graph.

        Returns:
            Generator on nodes
        """
        visited = set()
        def traverse(node):
            visited.add(node.id)
            if not node.is_sink():
                for child in node.children():
                    if not child.id in visited:
                        yield from traverse(child)
            yield node
        yield from traverse(self.root)

    def to_graphviz(self) -> gv.Digraph:
        r"""Converts a partition graph to graphviz (dot format) Digraph object.
        """
        G = gv.Digraph("Partition Graph", node_attr={"shape": "record"})
        for node in self.postorder():
            G.node(node.id)
            for edge_info in node.edges.values():
                if edge_info is not None:
                    child, val = edge_info
                    G.edge(node.id, child.id, label=str(val))
        return G

    # TODO: Rewrite this to greedily take the edge that maximizes MCP!

    def find_max_cost_path(self) -> PartitionPath:
        ...
    
    # def find_max_cost_path(self) -> PartitionPath:
    #     # Stack stores tuples of node and their level in this DFS tree.
    #     stack = []
    #     stack.append((self.root, self.root.value, 0))
    #     path = PartitionPath(self)
    #     best_path = None

    #     i = 0
    
    #     # loop till queue is empty
    #     while len(stack) != 0:
    #         curr, val, level = stack.pop()
    #         path.append(curr, val)

    #         ## DEBUG ######################################
    #         if len(path) > self.num_leaves:
    #             print(f"{len(path)} > {self.num_leaves}")
    #             print(len(path.get_partition_ids()))
    #         assert len(path) <= self.num_leaves
    #         ###############################################

    #         # print(val, [node.id for node in path.partition_list])
    #         if i % 1000 == 0:
    #             print(f"{i} / {self.num_edges}\tlen:", len(path), "val:", val, "id:",curr.id)
    #             if best_path is not None:
    #                 print(f"\tbest", best_path.get_cost())
    #         i += 1

    #         if curr.is_sink():
    #             # if best_path is not None:
    #             #     print(f"=> Valid partition\tval:{path.get_cost()}\tlen:{len(path)}\t{best_path.get_partition_ids()}")
    #             # Full partition and best path candidate
    #             if len(path) < self.max_clade_size and (best_path is None or path.cost > best_path.cost):
    #                 best_path = path.copy()
    #                 # print(f"\t current best: val {val}, len {len(path)} \t{[node.id for node in best_path.partition_list]}")

    #             if len(stack) > 0:
    #                 # print("\tRewinding", 1+level-stack[-1][2])
    #                 path.rewind(1+level-stack[-1][2])
            
    #         else:
    #             for edge_type, edge in curr.edges.items():
    #                 if edge is None:
    #                     continue
    #                 child, val = edge
    #                 if edge_type == "up" and not path.legal_up_edge(child):
    #                     continue

    #                 if edge_type == "up":
    #                     ... # TODO: Take this edge

    #                 # if edge_type == "up" and not path.beneficial_up_edge(child, val): # Greedy heuristic
    #                 #     continue
                
    #                 stack.append((child, val, level+1))

    #     return best_path

    # def find_max_cost_path(self) -> PartitionPath:
    #     # Stack stores tuples of node and their level in this DFS tree.
    #     stack = []
    #     stack.append((self.root, self.root.value, 0))
    #     path = PartitionPath(self)
    #     best_path = None

    #     i = 0
    
    #     # loop till queue is empty
    #     while len(stack) != 0:
    #         curr, val, level = stack.pop()
    #         path.append(curr, val)

    #         ## DEBUG ######################################
    #         if len(path) > self.num_leaves:
    #             print(f"{len(path)} > {self.num_leaves}")
    #             print(len(path.get_partition_ids()))
    #         assert len(path) <= self.num_leaves
    #         ###############################################

    #         # print(val, [node.id for node in path.partition_list])
    #         if i % 1000 == 0:
    #             print(f"{i} / {self.num_edges}\tlen:", len(path), "val:", val, "id:",curr.id)
    #             if best_path is not None:
    #                 print(f"\tbest", best_path.get_cost())
    #         i += 1

    #         if curr.is_sink():
    #             # if best_path is not None:
    #             #     print(f"=> Valid partition\tval:{path.get_cost()}\tlen:{len(path)}\t{best_path.get_partition_ids()}")
    #             # Full partition and best path candidate
    #             if len(path) < self.max_clade_size and (best_path is None or path.cost > best_path.cost):
    #                 best_path = path.copy()
    #                 # print(f"\t current best: val {val}, len {len(path)} \t{[node.id for node in best_path.partition_list]}")

    #             if len(stack) > 0:
    #                 # print("\tRewinding", 1+level-stack[-1][2])
    #                 path.rewind(1+level-stack[-1][2])
            
    #         else:
    #             for edge_type, edge in curr.edges.items():
    #                 if edge is None:
    #                     continue
    #                 child, val = edge
    #                 if edge_type == "up" and not path.legal_up_edge(child):
    #                     continue

    #                 if edge_type == "up":
    #                     ... # TODO: Take this edge

    #                 # if edge_type == "up" and not path.beneficial_up_edge(child, val): # Greedy heuristic
    #                 #     continue
                
    #                 stack.append((child, val, level+1))

    #     return best_path




