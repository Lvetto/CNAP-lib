from old.utils import adjacency_matrix
import numpy as np

def get_all_bonds(adj_mat):
    """
    Make a list of all unique bonded atom indices

    Parameters:
    adj_mat:    np.array
        the adjacency matrix of the system
    
    Returns:
    unique_pairs:   np.array
        an array of unique pairs of bonded atoms
    
    """

    pairs = np.argwhere(adj_mat == 1)   # extract (i,j) where a_ij=1
    sorted_pairs = np.sort(pairs, axis=1)   # sort each pair of indices
    unique_pairs = np.unique(sorted_pairs, axis=0)  # remove duplicates

    return unique_pairs

def get_adj_list(adj_mat):
    """
    Make an adjacency list from an adjacency matrix

    Parameters:
    adj_mat:    np.array
        the adjacency matrix of the system
    
    Returns:
    adj list:   np.array
        a list of lists, (i don't know how to word the rest)
    
    """

    return [np.where(row)[0] for row in adj_mat == 1]

def get_common_neighbors(adj_list, bonds):
    """
    Make a dictionary of lists of common neighbours
    Element [(i, j)] would be the common neighbours for atoms i,j

    Parameters:
    adj_list:   np.array
        the ajacency list for the system
    bonds:  np.array
        a list of all bonded atom pairs
    
    Returns:
    CNs: dict
        a dict of common neighbours
    """

    CNs = {(i,j) : np.intersect1d(adj_list[i], adj_list[j]) for i,j in bonds}   # common neighbours for the two atoms are the intersection of their adjacency list
    return CNs

def get_common_bonds(adj_mat, cns):
    """
    Make a list of unique bonded atom pairs with atoms from a local environment

    Parameters:
    adj_list:   np.array
        the ajacency list for the system
    cns:  np.array
        the common neighbors to consider
    
    Returns:
    CBs: list
        a list of bonds between common neighbours
    """

    cns = np.array(cns)     # cns must be a  numpy array for this to work
    sub_adj_mat = get_adj_sub_mat(adj_mat, cns) # get syb matrix describing the common neighbors
    bonds = np.array(get_all_bonds(sub_adj_mat))  # extract nodes from the matrix

    # change the indices to refer to the whole system
    bonds = cns[bonds]

    return bonds

def get_adj_sub_mat(adj_mat, cns):
    """
    Create and adjacency matrix for a subset of nodes starting from a matrix describing a bigger system

    Args:
        adj_mat (np.array): adjacency matrix for the bigger system
        cns (list): indices for the subset of nodes to consider

    Returns:
        np.array: adjacency matrix for the subset of nodes
    """

    cns = np.array(cns)     # cns must be a  numpy array for this to work
    sub_adj_mat = adj_mat[np.ix_(cns, cns)] # make an adjacency matrix for the cns by extracting elements i,j where both i and j are in cns
    np.fill_diagonal(sub_adj_mat, 0)    # set the diagonal to 0. Avoids some issues when extracting bonds

    return sub_adj_mat


# VERY experimental!

# It should be something similar to this: https://en.wikipedia.org/wiki/Depth-first_search
# I pretty much copied the pseudocode and turned it into python, then added a few lines to keep track of the maximum lenght (depth?)
def explore_graph(adj_list, starting_node):
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
            for w in adj_list[v]:
                if (w not in discovered):
                    stack.append((w, l+1))
    
    return max_lenght

def longest_chain_lenght(adj_list):
    max_lenght = 0
    for i,_ in enumerate(adj_list):
        max_lenght = max(explore_graph(adj_list, i), max_lenght)
    return max_lenght

def get_signature(adj_mat):
    adj_list = get_adj_list(adj_mat)
    bonds = get_all_bonds(adj_mat)
    bonds = [tuple(i) for i in bonds]   # make the bonds hashable

    common_neighbors = get_common_neighbors(adj_list, bonds)
    common_bonds = {bond: get_common_bonds(adj_mat, cns) for bond, cns in zip(common_neighbors.keys(), common_neighbors.values())}
 
    adj_mats = {bond: get_adj_sub_mat(adj_mat, common_neighbors[bond]) for bond in bonds}
    adj_lists = {bond: get_adj_list(adj_mats[bond]) for bond in bonds}
    longest_chains = {bond: longest_chain_lenght(adj_lists[bond]) for bond in bonds}

    signatures = {bond: (len(common_neighbors[bond]), len(common_bonds[bond]), longest_chains[bond]) for bond in bonds}

    return signatures
