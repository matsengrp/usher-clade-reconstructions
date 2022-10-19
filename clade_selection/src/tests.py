import bte as mat
from partition_graph import PartitionPath, PartitionGraph
import sys

class recursion_depth:
    def __init__(self, limit):
        self.limit = limit
        self.default_limit = sys.getrecursionlimit()
    def __enter__(self):
        sys.setrecursionlimit(self.limit)
    def __exit__(self, type, value, traceback):
        sys.setrecursionlimit(self.default_limit)

def build_tree(nw, muts):
    t = mat.MATree()
    t.from_newick_string(nw)
    node_list = t.depth_first_expansion(reverse=True)
    mut_map = {}
    for node in node_list:
        mut_map[node.id] = muts
    t.apply_mutations(mut_map)
    return t

def test_dfs():
    # Single mutation means that the value score should always equal the path length
    muts= ["A0A"]
    nws = [
        "(A, B, C, D)1;",
        "((A, B, C)2, D)1;",
        "((A, B)2, (C, D)3)1;",
    ]
    partitions = {}
    for clade_restriction in [5, 4, 3, 2, 1]:
        partitions[clade_restriction] = []
        for nw in nws:
            t = build_tree(nw, muts)
            print("***", clade_restriction, nw)
            try:
                pg = PartitionGraph(t, max_clade_size=clade_restriction)
                print(pg.to_graphviz())
                best_path = pg.find_max_cost_path()
                if best_path is not None:
                    partitions[clade_restriction].append((best_path.get_cost(), len(best_path), best_path.get_partition_ids(), nw))
                    assert best_path.get_cost() == len(best_path)
            except:
                print("\nError!\n")

    for k, v in partitions.items():
        print(k, v)
    return partitions


def test_graph_building():
    pg = build_small_graph()
    G = pg.to_graphviz()
    print(G)

    pp = pg.find_max_cost_path()
    partition = pp.get_partition_ids()
    print("partition:", partition)

    pgs = build_graphs()
    for pg in pgs:
        G = pg.to_graphviz()
        print("\n\n\n")
        print(G)
        pp = pg.find_max_cost_path()
        partition = pp.get_partition_ids()
        print("partition:", partition)

def build_graphs():
    # Multifurcation + unifurcation
    muts= ["A0A", "A0A"]
    t = mat.MATree()

    nw = "((((A,B)1)2)3, C, D, (E, F))4;"
    t.from_newick_string(nw)
    node_list = t.depth_first_expansion(reverse=True)
    mut_map = {}
    for node in node_list:
        mut_map[node.id] = muts
    t.apply_mutations(mut_map)

    return [PartitionGraph(t, max_clade_size=6)]


def build_small_graph():
    muts= ["A0A", "A0A"]
    t = mat.MATree()

    nw = "(B,(C,D)A)Root;"
    t.from_newick_string(nw)
    node_list = t.depth_first_expansion(reverse=True)
    mut_map = {}
    for node in node_list:
        mut_map[node.id] = muts
    t.apply_mutations(mut_map)

    return PartitionGraph(t)


def build_clade_graph():
    clade_path = "clades/20F/subset_mat.pb"
    matree = mat.MATree(clade_path)

    print("Building PG on MAT...")
    pg = PartitionGraph(matree, 10000)
    with recursion_depth(20000):
        print(pg.to_graphviz())
    print(f"Graph contains {pg.num_edges} edges")
    print(f"Graph contains {pg.num_paths} paths")
    print(f"Graph contains {pg.num_leaves} leaves")

    print("Searching for best path...")
    best_path = pg.find_max_cost_path()
    print(f"\tcost:{best_path.get_cost()}\tlen:{len(best_path)}")



def main():
    # test_graph_building()
    # test_dfs()
    build_clade_graph()

    

if __name__ == "__main__":
    main()