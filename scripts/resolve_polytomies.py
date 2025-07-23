import sys
import treeswift

def resolve_polytomies(ref_tree, res_tree):
    tree = treeswift.read_tree_newick(ref_tree)
    treeswift.Tree.resolve_polytomies(tree)
    treeswift.Tree.write_tree_newick(tree, filename=res_tree)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py ref_tree res_tree")
        sys.exit(1)
    
    ref_tree = sys.argv[1]
    res_tree = sys.argv[2]

    resolve_polytomies(ref_tree, res_tree)
