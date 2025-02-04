from old.utils import adjacency_matrix
import numpy as np

class Graph:
    def __init__(self, adj_mat, node_positions, fill_diag=True):
        self.adj_mat = adj_mat

        if (fill_diag):
            np.fill_diagonal(self.adj_mat, 0)

        self.adj_list = self.make_adj_list()
        self.positions = node_positions

    def make_adj_list(self):
        return [np.where(row)[0] for row in self.adj_mat == 1]

    def make_subgraph(self, nodes):
        nodes = np.array(nodes)
        sub_adj_mat = self.adj_mat[np.ix_(nodes, nodes)] # make an adjacency matrix for the nodes of the subgraph by extracting elements i,j where both i and j are in the subgraph
        np.fill_diagonal(sub_adj_mat, 0)    # set the diagonal to 0. Avoids some issues when extracting bonds

        return Graph(sub_adj_mat, self.positions[nodes])
    
    def longest_chain_from_node(self, starting_node):
        max_lenght = 0
        discovered = []
        stack = []

        stack.append((starting_node, 0))
        while (len(stack) > 0):
            v, l = stack[-1]
            stack = stack[:-1]

            if (v not in discovered):
                discovered.append(v)
                max_lenght = max(max_lenght, l)

                for w in self.adj_list[v]:
                    if (w not in discovered):
                        stack.append((w, l+1))

        return max_lenght

    @property
    def longest_chain(self):
        max_lenght = 0
        for i,_ in enumerate(self.adj_list):
            max_lenght = max(self.longest_chain_from_node(i), max_lenght)

        return max_lenght

    @property
    def unique_bonds(self):
        pairs = np.argwhere(self.adj_mat == 1)   # extract (i,j) where a_ij=1
        sorted_pairs = np.sort(pairs, axis=1)   # sort each pair of indices
        unique_pairs = np.unique(sorted_pairs, axis=0)  # remove duplicates

        return unique_pairs

    @property
    def number_of_unique_bonds(self):
        return len(self.unique_bonds)

    @property
    def number_of_nodes(self):
        return self.adj_mat.shape[0]

def get_cns_subgraphs(graph):
    subgraphs = [graph.make_subgraph(np.intersect1d(graph.adj_list[i], graph.adj_list[j])) for i,j in graph.unique_bonds]
    return subgraphs

def compute_signatures(graph):
    subgraphs = get_cns_subgraphs(graph)

    signatures = []
    for subgraph in subgraphs:
        i = subgraph.number_of_nodes
        j = subgraph.number_of_unique_bonds
        k = subgraph.longest_chain
        signatures.append((i,j,k))

    return signatures

def get_occurrences(signatures):
    unique_sigs, counts = np.unique(signatures, axis=0, return_counts=True)

    return {tuple(sig): count for sig, count in zip(unique_sigs, counts)}





# doesn't work
def get_feature_vectors(signatures):
    feature_vectors = {}

    for bond, signature in zip(signatures.keys(), signatures.values()):
        i,j = bond

        # make sure atoms i,j are in the dict
        if (i not in feature_vectors.keys()):
            feature_vectors[i] = {}

        if (j not in feature_vectors.keys()):
            feature_vectors[j] = {}

        # if signature was already encountered for atom i, add 1 to the number of occurences, otherwise set it to 1
        if (signature not in feature_vectors[i].keys()):
            feature_vectors[i][signature] = 1
        else:
            feature_vectors[i][signature] += 1

        # if signature was already encountered for atom j, add 1 to the number of occurences, otherwise set it to 1
        if (signature not in feature_vectors[j].keys()):
            feature_vectors[j][signature] = 1
        else:
            feature_vectors[j][signature] += 1

    return feature_vectors

