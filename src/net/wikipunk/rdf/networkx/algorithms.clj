(ns net.wikipunk.rdf.networkx.algorithms
  {:rdf/type :owl/Ontology}
  (:require   
   [net.wikipunk.rdf.py])
  (:refer-clojure :exclude [ancestors complement descendants load reverse]))

(def adamic_adar_index
  "Compute the Adamic-Adar index of all node pairs in ebunch."
  {:db/ident :networkx.algorithms/adamic_adar_index,
   :rdf/type :py/Function})

(def all
  "Operations on many graphs."
  {:db/ident :networkx.algorithms/all,
   :rdf/type :py/Function})

(def all_node_cuts
  "Returns all minimum k cutsets of an undirected graph G."
  {:db/ident :networkx.algorithms/all_node_cuts,
   :rdf/type :py/Function})

(def all_pairs_bellman_ford_path
  "Compute shortest paths between all nodes in a weighted graph."
  {:db/ident :networkx.algorithms/all_pairs_bellman_ford_path,
   :rdf/type :py/Function})

(def all_pairs_bellman_ford_path_length
  "Compute shortest path lengths between all nodes in a weighted graph."
  {:db/ident :networkx.algorithms/all_pairs_bellman_ford_path_length,
   :rdf/type :py/Function})

(def all_pairs_dijkstra
  "Find shortest weighted paths and lengths between all nodes."
  {:db/ident :networkx.algorithms/all_pairs_dijkstra,
   :rdf/type :py/Function})

(def all_pairs_dijkstra_path
  "Compute shortest paths between all nodes in a weighted graph."
  {:db/ident :networkx.algorithms/all_pairs_dijkstra_path,
   :rdf/type :py/Function})

(def all_pairs_dijkstra_path_length
  "Compute shortest path lengths between all nodes in a weighted graph."
  {:db/ident :networkx.algorithms/all_pairs_dijkstra_path_length,
   :rdf/type :py/Function})

(def all_pairs_lowest_common_ancestor
  "Return the lowest common ancestor of all pairs or the provided pairs"
  {:db/ident :networkx.algorithms/all_pairs_lowest_common_ancestor,
   :rdf/type :py/Function})

(def all_pairs_node_connectivity
  "Compute node connectivity between all pairs of nodes of G."
  {:db/ident :networkx.algorithms/all_pairs_node_connectivity,
   :rdf/type :py/Function})

(def all_pairs_shortest_path
  "Compute shortest paths between all nodes."
  {:db/ident :networkx.algorithms/all_pairs_shortest_path,
   :rdf/type :py/Function})

(def all_pairs_shortest_path_length
  "Computes the shortest path lengths between all nodes in `G`."
  {:db/ident :networkx.algorithms/all_pairs_shortest_path_length,
   :rdf/type :py/Function})

(def all_shortest_paths
  "Compute all shortest simple paths in the graph."
  {:db/ident :networkx.algorithms/all_shortest_paths,
   :rdf/type :py/Function})

(def all_simple_edge_paths
  "Generate lists of edges for all simple paths in G from source to target."
  {:db/ident :networkx.algorithms/all_simple_edge_paths,
   :rdf/type :py/Function})

(def all_simple_paths
  "Generate all simple paths in the graph G from source to target."
  {:db/ident :networkx.algorithms/all_simple_paths,
   :rdf/type :py/Function})

(def all_topological_sorts
  "Returns a generator of _all_ topological sorts of the directed graph G."
  {:db/ident :networkx.algorithms/all_topological_sorts,
   :rdf/type :py/Function})

(def all_triads
  "A generator of all possible triads in G."
  {:db/ident :networkx.algorithms/all_triads,
   :rdf/type :py/Function})

(def all_triplets
  "Returns a generator of all possible sets of 3 nodes in a DiGraph."
  {:db/ident :networkx.algorithms/all_triplets,
   :rdf/type :py/Function})

(def ancestors
  "Returns all nodes having a path to `source` in `G`."
  {:db/ident :networkx.algorithms/ancestors,
   :rdf/type :py/Function})

(def antichains
  "Generates antichains from a directed acyclic graph (DAG)."
  {:db/ident :networkx.algorithms/antichains,
   :rdf/type :py/Function})

(def approximate_current_flow_betweenness_centrality
  "Compute the approximate current-flow betweenness centrality for nodes."
  {:db/ident
   :networkx.algorithms/approximate_current_flow_betweenness_centrality,
   :rdf/type :py/Function})

(def approximation
  "Approximations of graph properties and Heuristic methods for optimization."
  {:db/ident :networkx.algorithms/approximation,
   :rdf/type :py/Function})

(def articulation_points
  "Yield the articulation points, or cut vertices, of a graph."
  {:db/ident :networkx.algorithms/articulation_points,
   :rdf/type :py/Function})

(def assortativity
  "assortativity"
  {:db/ident :networkx.algorithms/assortativity,
   :rdf/type :py/Function})

(def astar
  "Shortest paths and path lengths using the A* (\"A star\") algorithm."
  {:db/ident :networkx.algorithms/astar,
   :rdf/type :py/Function})

(def astar_path
  "Returns a list of nodes in a shortest path between source and target using the A* (\"A-star\") algorithm."
  {:db/ident :networkx.algorithms/astar_path,
   :rdf/type :py/Function})

(def astar_path_length
  "Returns the length of the shortest path between source and target using the A* (\"A-star\") algorithm."
  {:db/ident :networkx.algorithms/astar_path_length,
   :rdf/type :py/Function})

(def asteroidal
  "Algorithms for asteroidal triples and asteroidal numbers in graphs."
  {:db/ident :networkx.algorithms/asteroidal,
   :rdf/type :py/Function})

(def attracting
  "Attracting components."
  {:db/ident :networkx.algorithms/attracting,
   :rdf/type :py/Function})

(def attracting_components
  "Generates the attracting components in `G`."
  {:db/ident :networkx.algorithms/attracting_components,
   :rdf/type :py/Function})

(def attribute_assortativity_coefficient
  "Compute assortativity for node attributes."
  {:db/ident :networkx.algorithms/attribute_assortativity_coefficient,
   :rdf/type :py/Function})

(def attribute_mixing_dict
  "Returns dictionary representation of mixing matrix for attribute."
  {:db/ident :networkx.algorithms/attribute_mixing_dict,
   :rdf/type :py/Function})

(def attribute_mixing_matrix
  "Returns mixing matrix for attribute."
  {:db/ident :networkx.algorithms/attribute_mixing_matrix,
   :rdf/type :py/Function})

(def average_clustering
  "Compute the average clustering coefficient for the graph G."
  {:db/ident :networkx.algorithms/average_clustering,
   :rdf/type :py/Function})

(def average_degree_connectivity
  "Compute the average degree connectivity of graph."
  {:db/ident :networkx.algorithms/average_degree_connectivity,
   :rdf/type :py/Function})

(def average_neighbor_degree
  "Returns the average degree of the neighborhood of each node."
  {:db/ident :networkx.algorithms/average_neighbor_degree,
   :rdf/type :py/Function})

(def average_node_connectivity
  "Returns the average connectivity of a graph G."
  {:db/ident :networkx.algorithms/average_node_connectivity,
   :rdf/type :py/Function})

(def average_shortest_path_length
  "Returns the average shortest path length."
  {:db/ident :networkx.algorithms/average_shortest_path_length,
   :rdf/type :py/Function})

(def barycenter
  "Calculate barycenter of a connected graph, optionally with edge weights."
  {:db/ident :networkx.algorithms/barycenter,
   :rdf/type :py/Function})

(def beamsearch
  "Basic algorithms for breadth-first searching the nodes of a graph."
  {:db/ident :networkx.algorithms/beamsearch,
   :rdf/type :py/Function})

(def bellman_ford_path
  "Returns the shortest path from source to target in a weighted graph G."
  {:db/ident :networkx.algorithms/bellman_ford_path,
   :rdf/type :py/Function})

(def bellman_ford_path_length
  "Returns the shortest path length from source to target in a weighted graph."
  {:db/ident :networkx.algorithms/bellman_ford_path_length,
   :rdf/type :py/Function})

(def bellman_ford_predecessor_and_distance
  "Compute shortest path lengths and predecessors on shortest paths in weighted graphs."
  {:db/ident :networkx.algorithms/bellman_ford_predecessor_and_distance,
   :rdf/type :py/Function})

(def betweenness
  "Betweenness centrality measures."
  {:db/ident :networkx.algorithms/betweenness,
   :rdf/type :py/Function})

(def betweenness_centrality
  "Compute the shortest-path betweenness centrality for nodes."
  {:db/ident :networkx.algorithms/betweenness_centrality,
   :rdf/type :py/Function})

(def betweenness_centrality_subset
  "Compute betweenness centrality for a subset of nodes."
  {:db/ident :networkx.algorithms/betweenness_centrality_subset,
   :rdf/type :py/Function})

(def betweenness_subset
  "Betweenness centrality measures for subsets of nodes."
  {:db/ident :networkx.algorithms/betweenness_subset,
   :rdf/type :py/Function})

(def bfs_beam_edges
  "Iterates over edges in a beam search."
  {:db/ident :networkx.algorithms/bfs_beam_edges,
   :rdf/type :py/Function})

(def bfs_edges
  "Iterate over edges in a breadth-first-search starting at source."
  {:db/ident :networkx.algorithms/bfs_edges,
   :rdf/type :py/Function})

(def bfs_layers
  "Returns an iterator of all the layers in breadth-first search traversal."
  {:db/ident :networkx.algorithms/bfs_layers,
   :rdf/type :py/Function})

(def bfs_predecessors
  "Returns an iterator of predecessors in breadth-first-search from source."
  {:db/ident :networkx.algorithms/bfs_predecessors,
   :rdf/type :py/Function})

(def bfs_successors
  "Returns an iterator of successors in breadth-first-search from source."
  {:db/ident :networkx.algorithms/bfs_successors,
   :rdf/type :py/Function})

(def bfs_tree
  "Returns an oriented tree constructed from of a breadth-first-search starting at source."
  {:db/ident :networkx.algorithms/bfs_tree,
   :rdf/type :py/Function})

(def biconnected
  "Biconnected components and articulation points."
  {:db/ident :networkx.algorithms/biconnected,
   :rdf/type :py/Function})

(def biconnected_component_edges
  "Returns a generator of lists of edges, one list for each biconnected component of the input graph."
  {:db/ident :networkx.algorithms/biconnected_component_edges,
   :rdf/type :py/Function})

(def biconnected_components
  "Returns a generator of sets of nodes, one set for each biconnected component of the graph"
  {:db/ident :networkx.algorithms/biconnected_components,
   :rdf/type :py/Function})

(def bidirectional_dijkstra
  "Dijkstra's algorithm for shortest paths using bidirectional search."
  {:db/ident :networkx.algorithms/bidirectional_dijkstra,
   :rdf/type :py/Function})

(def bidirectional_shortest_path
  "Returns a list of nodes in a shortest path between source and target."
  {:db/ident :networkx.algorithms/bidirectional_shortest_path,
   :rdf/type :py/Function})

(def binary
  "Operations on graphs including union, intersection, difference."
  {:db/ident :networkx.algorithms/binary,
   :rdf/type :py/Function})

(def bipartite
  "This module provides functions and operations for bipartite graphs.  Bipartite graphs `B = (U, V, E)` have two node sets `U,V` and edges in `E` that only connect nodes from opposite sets. It is common in the literature to use an spatial analogy referring to the two node sets as top and bottom nodes."
  {:db/ident :networkx.algorithms/bipartite,
   :rdf/type :py/Function})

(def boundary
  "Routines to find the boundary of a set of nodes."
  {:db/ident :networkx.algorithms/boundary,
   :rdf/type :py/Function})

(def boundary_expansion
  "Returns the boundary expansion of the set `S`."
  {:db/ident :networkx.algorithms/boundary_expansion,
   :rdf/type :py/Function})

(def breadth_first_search
  "Basic algorithms for breadth-first searching the nodes of a graph."
  {:db/ident :networkx.algorithms/breadth_first_search,
   :rdf/type :py/Function})

(def bridges
  "Generate all bridges in a graph."
  {:db/ident :networkx.algorithms/bridges,
   :rdf/type :py/Function})

(def capacity_scaling
  "Find a minimum cost flow satisfying all demands in digraph G."
  {:db/ident :networkx.algorithms/capacity_scaling,
   :rdf/type :py/Function})

(def cartesian_product
  "Returns the Cartesian product of G and H."
  {:db/ident :networkx.algorithms/cartesian_product,
   :rdf/type :py/Function})

(def center
  "Returns the center of the graph G."
  {:db/ident :networkx.algorithms/center,
   :rdf/type :py/Function})

(def centrality
  "centrality"
  {:db/ident :networkx.algorithms/centrality,
   :rdf/type :py/Function})

(def chain_decomposition
  "Returns the chain decomposition of a graph."
  {:db/ident :networkx.algorithms/chain_decomposition,
   :rdf/type :py/Function})

(def chains
  "Functions for finding chains in a graph."
  {:db/ident :networkx.algorithms/chains,
   :rdf/type :py/Function})

(def check_planarity
  "Check if a graph is planar and return a counterexample or an embedding."
  {:db/ident :networkx.algorithms/check_planarity,
   :rdf/type :py/Function})

(def chordal
  "Algorithms for chordal graphs."
  {:db/ident :networkx.algorithms/chordal,
   :rdf/type :py/Function})

(def chordal_graph_cliques
  "Returns all maximal cliques of a chordal graph."
  {:db/ident :networkx.algorithms/chordal_graph_cliques,
   :rdf/type :py/Function})

(def chordal_graph_treewidth
  "Returns the treewidth of the chordal graph G."
  {:db/ident :networkx.algorithms/chordal_graph_treewidth,
   :rdf/type :py/Function})

(def chromatic_polynomial
  "Returns the chromatic polynomial of `G`"
  {:db/ident :networkx.algorithms/chromatic_polynomial,
   :rdf/type :py/Function})

(def clique
  "Functions for finding and manipulating cliques."
  {:db/ident :networkx.algorithms/clique,
   :rdf/type :py/Function})

(def cliques_containing_node
  "Returns a list of cliques containing the given node."
  {:db/ident :networkx.algorithms/cliques_containing_node,
   :rdf/type :py/Function})

(def closeness
  "Closeness centrality measures."
  {:db/ident :networkx.algorithms/closeness,
   :rdf/type :py/Function})

(def closeness_centrality
  "Compute closeness centrality for nodes."
  {:db/ident :networkx.algorithms/closeness_centrality,
   :rdf/type :py/Function})

(def closeness_vitality
  "Returns the closeness vitality for nodes in the graph."
  {:db/ident :networkx.algorithms/closeness_vitality,
   :rdf/type :py/Function})

(def cluster
  "Algorithms to characterize the number of triangles in a graph."
  {:db/ident :networkx.algorithms/cluster,
   :rdf/type :py/Function})

(def clustering
  "Compute the clustering coefficient for nodes."
  {:db/ident :networkx.algorithms/clustering,
   :rdf/type :py/Function})

(def cn_soundarajan_hopcroft
  "Count the number of common neighbors of all node pairs in ebunch using community information."
  {:db/ident :networkx.algorithms/cn_soundarajan_hopcroft,
   :rdf/type :py/Function})

(def coloring
  "coloring"
  {:db/ident :networkx.algorithms/coloring,
   :rdf/type :py/Function})

(def combinatorial_embedding_to_pos
  "Assigns every node a (x, y) position based on the given embedding"
  {:db/ident :networkx.algorithms/combinatorial_embedding_to_pos,
   :rdf/type :py/Function})

(def common_neighbor_centrality
  "Return the CCPA score for each pair of nodes."
  {:db/ident :networkx.algorithms/common_neighbor_centrality,
   :rdf/type :py/Function})

(def communicability
  "Returns communicability between all pairs of nodes in G."
  {:db/ident :networkx.algorithms/communicability,
   :rdf/type :py/Function})

(def communicability_alg
  "Communicability."
  {:db/ident :networkx.algorithms/communicability_alg,
   :rdf/type :py/Function})

(def communicability_betweenness_centrality
  "Returns subgraph communicability for all pairs of nodes in G."
  {:db/ident :networkx.algorithms/communicability_betweenness_centrality,
   :rdf/type :py/Function})

(def communicability_exp
  "Returns communicability between all pairs of nodes in G."
  {:db/ident :networkx.algorithms/communicability_exp,
   :rdf/type :py/Function})

(def community
  "Functions for computing and measuring community structure."
  {:db/ident :networkx.algorithms/community,
   :rdf/type :py/Function})

(def complement
  "Returns the graph complement of G."
  {:db/ident :networkx.algorithms/complement,
   :rdf/type :py/Function})

(def complete_bipartite_graph
  "Returns the complete bipartite graph `K_{n_1,n_2}`."
  {:db/ident :networkx.algorithms/complete_bipartite_graph,
   :rdf/type :py/Function})

(def complete_to_chordal_graph
  "Return a copy of G completed to a chordal graph"
  {:db/ident :networkx.algorithms/complete_to_chordal_graph,
   :rdf/type :py/Function})

(def components
  "components"
  {:db/ident :networkx.algorithms/components,
   :rdf/type :py/Function})

(def compose
  "Compose graph G with H by combining nodes and edges into a single graph."
  {:db/ident :networkx.algorithms/compose,
   :rdf/type :py/Function})

(def compose_all
  "Returns the composition of all graphs."
  {:db/ident :networkx.algorithms/compose_all,
   :rdf/type :py/Function})

(def compute_v_structures
  "Iterate through the graph to compute all v-structures."
  {:db/ident :networkx.algorithms/compute_v_structures,
   :rdf/type :py/Function})

(def condensation
  "Returns the condensation of G."
  {:db/ident :networkx.algorithms/condensation,
   :rdf/type :py/Function})

(def conductance
  "Returns the conductance of two sets of nodes."
  {:db/ident :networkx.algorithms/conductance,
   :rdf/type :py/Function})

(def connected
  "Connected components."
  {:db/ident :networkx.algorithms/connected,
   :rdf/type :py/Function})

(def connected_components
  "Generate connected components."
  {:db/ident :networkx.algorithms/connected_components,
   :rdf/type :py/Function})

(def connected_double_edge_swap
  "Attempts the specified number of double-edge swaps in the graph `G`."
  {:db/ident :networkx.algorithms/connected_double_edge_swap,
   :rdf/type :py/Function})

(def connectivity
  "Connectivity and cut algorithms"
  {:db/ident :networkx.algorithms/connectivity,
   :rdf/type :py/Function})

(def constraint
  "Returns the constraint on all nodes in the graph ``G``."
  {:db/ident :networkx.algorithms/constraint,
   :rdf/type :py/Function})

(def contracted_edge
  "Returns the graph that results from contracting the specified edge."
  {:db/ident :networkx.algorithms/contracted_edge,
   :rdf/type :py/Function})

(def contracted_nodes
  "Returns the graph that results from contracting `u` and `v`."
  {:db/ident :networkx.algorithms/contracted_nodes,
   :rdf/type :py/Function})

(def core
  "Find the k-cores of a graph."
  {:db/ident :networkx.algorithms/core,
   :rdf/type :py/Function})

(def core_number
  "Returns the core number for each vertex."
  {:db/ident :networkx.algorithms/core_number,
   :rdf/type :py/Function})

(def corona_product
  "Returns the Corona product of G and H."
  {:db/ident :networkx.algorithms/corona_product,
   :rdf/type :py/Function})

(def correlation
  "Node assortativity coefficients and correlation measures."
  {:db/ident :networkx.algorithms/correlation,
   :rdf/type :py/Function})

(def cost_of_flow
  "Compute the cost of the flow given by flowDict on graph G."
  {:db/ident :networkx.algorithms/cost_of_flow,
   :rdf/type :py/Function})

(def could_be_isomorphic
  "Returns False if graphs are definitely not isomorphic. True does NOT guarantee isomorphism."
  {:db/ident :networkx.algorithms/could_be_isomorphic,
   :rdf/type :py/Function})

(def covering
  "Functions related to graph covers."
  {:db/ident :networkx.algorithms/covering,
   :rdf/type :py/Function})

(def current_flow_betweenness
  "Current-flow betweenness centrality measures."
  {:db/ident :networkx.algorithms/current_flow_betweenness,
   :rdf/type :py/Function})

(def current_flow_betweenness_centrality
  "Compute current-flow betweenness centrality for nodes."
  {:db/ident :networkx.algorithms/current_flow_betweenness_centrality,
   :rdf/type :py/Function})

(def current_flow_betweenness_centrality_subset
  "Compute current-flow betweenness centrality for subsets of nodes."
  {:db/ident :networkx.algorithms/current_flow_betweenness_centrality_subset,
   :rdf/type :py/Function})

(def current_flow_betweenness_subset
  "Current-flow betweenness centrality measures for subsets of nodes."
  {:db/ident :networkx.algorithms/current_flow_betweenness_subset,
   :rdf/type :py/Function})

(def current_flow_closeness
  "Current-flow closeness centrality measures."
  {:db/ident :networkx.algorithms/current_flow_closeness,
   :rdf/type :py/Function})

(def current_flow_closeness_centrality
  "Compute current-flow closeness centrality for nodes."
  {:db/ident :networkx.algorithms/current_flow_closeness_centrality,
   :rdf/type :py/Function})

(def cut_size
  "Returns the size of the cut between two sets of nodes."
  {:db/ident :networkx.algorithms/cut_size,
   :rdf/type :py/Function})

(def cuts
  "Functions for finding and evaluating cuts in a graph."
  {:db/ident :networkx.algorithms/cuts,
   :rdf/type :py/Function})

(def cycle_basis
  "Returns a list of cycles which form a basis for cycles of G."
  {:db/ident :networkx.algorithms/cycle_basis,
   :rdf/type :py/Function})

(def cycles
  "======================== Cycle finding algorithms ========================"
  {:db/ident :networkx.algorithms/cycles,
   :rdf/type :py/Function})

(def d_separated
  "Return whether node sets ``x`` and ``y`` are d-separated by ``z``."
  {:db/ident :networkx.algorithms/d_separated,
   :rdf/type :py/Function})

(def d_separation
  "Algorithm for testing d-separation in DAGs."
  {:db/ident :networkx.algorithms/d_separation,
   :rdf/type :py/Function})

(def dag
  "Algorithms for directed acyclic graphs (DAGs)."
  {:db/ident :networkx.algorithms/dag,
   :rdf/type :py/Function})

(def dag_longest_path
  "Returns the longest path in a directed acyclic graph (DAG)."
  {:db/ident :networkx.algorithms/dag_longest_path,
   :rdf/type :py/Function})

(def dag_longest_path_length
  "Returns the longest path length in a DAG"
  {:db/ident :networkx.algorithms/dag_longest_path_length,
   :rdf/type :py/Function})

(def dag_to_branching
  "Returns a branching representing all (overlapping) paths from root nodes to leaf nodes in the given directed acyclic graph."
  {:db/ident :networkx.algorithms/dag_to_branching,
   :rdf/type :py/Function})

(def dedensify
  "Compresses neighborhoods around high-degree nodes"
  {:db/ident :networkx.algorithms/dedensify,
   :rdf/type :py/Function})

(def degree_alg
  "Degree centrality measures."
  {:db/ident :networkx.algorithms/degree_alg,
   :rdf/type :py/Function})

(def degree_assortativity_coefficient
  "Compute degree assortativity of graph."
  {:db/ident :networkx.algorithms/degree_assortativity_coefficient,
   :rdf/type :py/Function})

(def degree_centrality
  "Compute the degree centrality for nodes."
  {:db/ident :networkx.algorithms/degree_centrality,
   :rdf/type :py/Function})

(def degree_mixing_dict
  "Returns dictionary representation of mixing matrix for degree."
  {:db/ident :networkx.algorithms/degree_mixing_dict,
   :rdf/type :py/Function})

(def degree_mixing_matrix
  "Returns mixing matrix for attribute."
  {:db/ident :networkx.algorithms/degree_mixing_matrix,
   :rdf/type :py/Function})

(def degree_pearson_correlation_coefficient
  "Compute degree assortativity of graph."
  {:db/ident :networkx.algorithms/degree_pearson_correlation_coefficient,
   :rdf/type :py/Function})

(def dense
  "Floyd-Warshall algorithm for shortest paths."
  {:db/ident :networkx.algorithms/dense,
   :rdf/type :py/Function})

(def depth_first_search
  "Basic algorithms for depth-first searching the nodes of a graph."
  {:db/ident :networkx.algorithms/depth_first_search,
   :rdf/type :py/Function})

(def descendants
  "Returns all nodes reachable from `source` in `G`."
  {:db/ident :networkx.algorithms/descendants,
   :rdf/type :py/Function})

(def descendants_at_distance
  "Returns all nodes at a fixed `distance` from `source` in `G`."
  {:db/ident :networkx.algorithms/descendants_at_distance,
   :rdf/type :py/Function})

(def dfs_edges
  "Iterate over edges in a depth-first-search (DFS)."
  {:db/ident :networkx.algorithms/dfs_edges,
   :rdf/type :py/Function})

(def dfs_labeled_edges
  "Iterate over edges in a depth-first-search (DFS) labeled by type."
  {:db/ident :networkx.algorithms/dfs_labeled_edges,
   :rdf/type :py/Function})

(def dfs_postorder_nodes
  "Generate nodes in a depth-first-search post-ordering starting at source."
  {:db/ident :networkx.algorithms/dfs_postorder_nodes,
   :rdf/type :py/Function})

(def dfs_predecessors
  "Returns dictionary of predecessors in depth-first-search from source."
  {:db/ident :networkx.algorithms/dfs_predecessors,
   :rdf/type :py/Function})

(def dfs_preorder_nodes
  "Generate nodes in a depth-first-search pre-ordering starting at source."
  {:db/ident :networkx.algorithms/dfs_preorder_nodes,
   :rdf/type :py/Function})

(def dfs_successors
  "Returns dictionary of successors in depth-first-search from source."
  {:db/ident :networkx.algorithms/dfs_successors,
   :rdf/type :py/Function})

(def dfs_tree
  "Returns oriented tree constructed from a depth-first-search from source."
  {:db/ident :networkx.algorithms/dfs_tree,
   :rdf/type :py/Function})

(def diameter
  "Returns the diameter of the graph G."
  {:db/ident :networkx.algorithms/diameter,
   :rdf/type :py/Function})

(def difference
  "Returns a new graph that contains the edges that exist in G but not in H."
  {:db/ident :networkx.algorithms/difference,
   :rdf/type :py/Function})

(def dijkstra_path
  "Returns the shortest weighted path from source to target in G."
  {:db/ident :networkx.algorithms/dijkstra_path,
   :rdf/type :py/Function})

(def dijkstra_path_length
  "Returns the shortest weighted path length in G from source to target."
  {:db/ident :networkx.algorithms/dijkstra_path_length,
   :rdf/type :py/Function})

(def dijkstra_predecessor_and_distance
  "Compute weighted shortest path length and predecessors."
  {:db/ident :networkx.algorithms/dijkstra_predecessor_and_distance,
   :rdf/type :py/Function})

(def directed_edge_swap
  "Swap three edges in a directed graph while keeping the node degrees fixed."
  {:db/ident :networkx.algorithms/directed_edge_swap,
   :rdf/type :py/Function})

(def disjoint_union
  "Combine graphs G and H. The nodes are assumed to be unique (disjoint)."
  {:db/ident :networkx.algorithms/disjoint_union,
   :rdf/type :py/Function})

(def disjoint_union_all
  "Returns the disjoint union of all graphs."
  {:db/ident :networkx.algorithms/disjoint_union_all,
   :rdf/type :py/Function})

(def dispersion
  "Calculate dispersion between `u` and `v` in `G`."
  {:db/ident :networkx.algorithms/dispersion,
   :rdf/type :py/Function})

(def distance_measures
  "Graph diameter, radius, eccentricity and other properties."
  {:db/ident :networkx.algorithms/distance_measures,
   :rdf/type :py/Function})

(def distance_regular
  "======================= Distance-regular graphs ======================="
  {:db/ident :networkx.algorithms/distance_regular,
   :rdf/type :py/Function})

(def dominance
  "Dominance algorithms."
  {:db/ident :networkx.algorithms/dominance,
   :rdf/type :py/Function})

(def dominance_frontiers
  "Returns the dominance frontiers of all nodes of a directed graph."
  {:db/ident :networkx.algorithms/dominance_frontiers,
   :rdf/type :py/Function})

(def dominating
  "Functions for computing dominating sets in a graph."
  {:db/ident :networkx.algorithms/dominating,
   :rdf/type :py/Function})

(def dominating_set
  "Finds a dominating set for the graph G."
  {:db/ident :networkx.algorithms/dominating_set,
   :rdf/type :py/Function})

(def double_edge_swap
  "Swap two edges in the graph while keeping the node degrees fixed."
  {:db/ident :networkx.algorithms/double_edge_swap,
   :rdf/type :py/Function})

(def eccentricity
  "Returns the eccentricity of nodes in G."
  {:db/ident :networkx.algorithms/eccentricity,
   :rdf/type :py/Function})

(def edge_betweenness_centrality
  "Compute betweenness centrality for edges."
  {:db/ident :networkx.algorithms/edge_betweenness_centrality,
   :rdf/type :py/Function})

(def edge_betweenness_centrality_subset
  "Compute betweenness centrality for edges for a subset of nodes."
  {:db/ident :networkx.algorithms/edge_betweenness_centrality_subset,
   :rdf/type :py/Function})

(def edge_bfs
  "A directed, breadth-first-search of edges in `G`, beginning at `source`."
  {:db/ident :networkx.algorithms/edge_bfs,
   :rdf/type :py/Function})

(def edge_boundary
  "Returns the edge boundary of `nbunch1`."
  {:db/ident :networkx.algorithms/edge_boundary,
   :rdf/type :py/Function})

(def edge_connectivity
  "Returns the edge connectivity of the graph or digraph G."
  {:db/ident :networkx.algorithms/edge_connectivity,
   :rdf/type :py/Function})

(def edge_current_flow_betweenness_centrality
  "Compute current-flow betweenness centrality for edges."
  {:db/ident :networkx.algorithms/edge_current_flow_betweenness_centrality,
   :rdf/type :py/Function})

(def edge_current_flow_betweenness_centrality_subset
  "Compute current-flow betweenness centrality for edges using subsets of nodes."
  {:db/ident
   :networkx.algorithms/edge_current_flow_betweenness_centrality_subset,
   :rdf/type :py/Function})

(def edge_dfs
  "A directed, depth-first-search of edges in `G`, beginning at `source`."
  {:db/ident :networkx.algorithms/edge_dfs,
   :rdf/type :py/Function})

(def edge_disjoint_paths
  "Returns the edges disjoint paths between source and target."
  {:db/ident :networkx.algorithms/edge_disjoint_paths,
   :rdf/type :py/Function})

(def edge_expansion
  "Returns the edge expansion between two node sets."
  {:db/ident :networkx.algorithms/edge_expansion,
   :rdf/type :py/Function})

(def edge_load_centrality
  "Compute edge load."
  {:db/ident :networkx.algorithms/edge_load_centrality,
   :rdf/type :py/Function})

(def edgebfs
  "============================= Breadth First Search on Edges ============================="
  {:db/ident :networkx.algorithms/edgebfs,
   :rdf/type :py/Function})

(def edgedfs
  "=========================== Depth First Search on Edges ==========================="
  {:db/ident :networkx.algorithms/edgedfs,
   :rdf/type :py/Function})

(def effective_size
  "Returns the effective size of all nodes in the graph ``G``."
  {:db/ident :networkx.algorithms/effective_size,
   :rdf/type :py/Function})

(def efficiency
  "Returns the efficiency of a pair of nodes in a graph."
  {:db/ident :networkx.algorithms/efficiency,
   :rdf/type :py/Function})

(def efficiency_measures
  "Provides functions for computing the efficiency of nodes and graphs."
  {:db/ident :networkx.algorithms/efficiency_measures,
   :rdf/type :py/Function})

(def eigenvector
  "Functions for computing eigenvector centrality."
  {:db/ident :networkx.algorithms/eigenvector,
   :rdf/type :py/Function})

(def eigenvector_centrality
  "Compute the eigenvector centrality for the graph `G`."
  {:db/ident :networkx.algorithms/eigenvector_centrality,
   :rdf/type :py/Function})

(def eigenvector_centrality_numpy
  "Compute the eigenvector centrality for the graph G."
  {:db/ident :networkx.algorithms/eigenvector_centrality_numpy,
   :rdf/type :py/Function})

(def enumerate_all_cliques
  "Returns all cliques in an undirected graph."
  {:db/ident :networkx.algorithms/enumerate_all_cliques,
   :rdf/type :py/Function})

(def equitable_color
  "Provides equitable (r + 1)-coloring for nodes of G in O(r * n^2) time if deg(G) <= r. The algorithm is described in [1]_."
  {:db/ident :networkx.algorithms/equitable_color,
   :rdf/type :py/Function})

(def equivalence_classes
  "Returns equivalence classes of `relation` when applied to `iterable`."
  {:db/ident :networkx.algorithms/equivalence_classes,
   :rdf/type :py/Function})

(def estrada_index
  "Returns the Estrada index of a the graph G."
  {:db/ident :networkx.algorithms/estrada_index,
   :rdf/type :py/Function})

(def euler
  "Eulerian circuits and graphs."
  {:db/ident :networkx.algorithms/euler,
   :rdf/type :py/Function})

(def eulerian_circuit
  "Returns an iterator over the edges of an Eulerian circuit in `G`."
  {:db/ident :networkx.algorithms/eulerian_circuit,
   :rdf/type :py/Function})

(def eulerian_path
  "Return an iterator over the edges of an Eulerian path in `G`."
  {:db/ident :networkx.algorithms/eulerian_path,
   :rdf/type :py/Function})

(def eulerize
  "Transforms a graph into an Eulerian graph."
  {:db/ident :networkx.algorithms/eulerize,
   :rdf/type :py/Function})

(def fast_could_be_isomorphic
  "Returns False if graphs are definitely not isomorphic."
  {:db/ident :networkx.algorithms/fast_could_be_isomorphic,
   :rdf/type :py/Function})

(def faster_could_be_isomorphic
  "Returns False if graphs are definitely not isomorphic."
  {:db/ident :networkx.algorithms/faster_could_be_isomorphic,
   :rdf/type :py/Function})

(def find_asteroidal_triple
  "Find an asteroidal triple in the given graph."
  {:db/ident :networkx.algorithms/find_asteroidal_triple,
   :rdf/type :py/Function})

(def find_cliques
  "Returns all maximal cliques in an undirected graph."
  {:db/ident :networkx.algorithms/find_cliques,
   :rdf/type :py/Function})

(def find_cliques_recursive
  "Returns all maximal cliques in a graph."
  {:db/ident :networkx.algorithms/find_cliques_recursive,
   :rdf/type :py/Function})

(def find_cycle
  "Returns a cycle found via depth-first traversal."
  {:db/ident :networkx.algorithms/find_cycle,
   :rdf/type :py/Function})

(def find_induced_nodes
  "Returns the set of induced nodes in the path from s to t."
  {:db/ident :networkx.algorithms/find_induced_nodes,
   :rdf/type :py/Function})

(def find_negative_cycle
  "Returns a cycle with negative total weight if it exists."
  {:db/ident :networkx.algorithms/find_negative_cycle,
   :rdf/type :py/Function})

(def flow
  "flow"
  {:db/ident :networkx.algorithms/flow,
   :rdf/type :py/Function})

(def flow_hierarchy
  "Returns the flow hierarchy of a directed network."
  {:db/ident :networkx.algorithms/flow_hierarchy,
   :rdf/type :py/Function})

(def flow_matrix
  "flow_matrix"
  {:db/ident :networkx.algorithms/flow_matrix,
   :rdf/type :py/Function})

(def floyd_warshall
  "Find all-pairs shortest path lengths using Floyd's algorithm."
  {:db/ident :networkx.algorithms/floyd_warshall,
   :rdf/type :py/Function})

(def floyd_warshall_numpy
  "Find all-pairs shortest path lengths using Floyd's algorithm."
  {:db/ident :networkx.algorithms/floyd_warshall_numpy,
   :rdf/type :py/Function})

(def floyd_warshall_predecessor_and_distance
  "Find all-pairs shortest path lengths using Floyd's algorithm."
  {:db/ident :networkx.algorithms/floyd_warshall_predecessor_and_distance,
   :rdf/type :py/Function})

(def from_nested_tuple
  "Returns the rooted tree corresponding to the given nested tuple."
  {:db/ident :networkx.algorithms/from_nested_tuple,
   :rdf/type :py/Function})

(def from_prufer_sequence
  "Returns the tree corresponding to the given Prüfer sequence."
  {:db/ident :networkx.algorithms/from_prufer_sequence,
   :rdf/type :py/Function})

(def full_join
  "Returns the full join of graphs G and H."
  {:db/ident :networkx.algorithms/full_join,
   :rdf/type :py/Function})

(def generalized_degree
  "Compute the generalized degree for nodes."
  {:db/ident :networkx.algorithms/generalized_degree,
   :rdf/type :py/Function})

(def generate_random_paths
  "Randomly generate `sample_size` paths of length `path_length`."
  {:db/ident :networkx.algorithms/generate_random_paths,
   :rdf/type :py/Function})

(def generic
  "Compute the shortest paths and path lengths between nodes in the graph."
  {:db/ident :networkx.algorithms/generic,
   :rdf/type :py/Function})

(def global_efficiency
  "Returns the average global efficiency of the graph."
  {:db/ident :networkx.algorithms/global_efficiency,
   :rdf/type :py/Function})

(def global_parameters
  "Returns global parameters for a given intersection array."
  {:db/ident :networkx.algorithms/global_parameters,
   :rdf/type :py/Function})

(def global_reaching_centrality
  "Returns the global reaching centrality of a directed graph."
  {:db/ident :networkx.algorithms/global_reaching_centrality,
   :rdf/type :py/Function})

(def goldberg_radzik
  "Compute shortest path lengths and predecessors on shortest paths in weighted graphs."
  {:db/ident :networkx.algorithms/goldberg_radzik,
   :rdf/type :py/Function})

(def gomory_hu_tree
  "Returns the Gomory-Hu tree of an undirected graph G."
  {:db/ident :networkx.algorithms/gomory_hu_tree,
   :rdf/type :py/Function})

(def google_matrix
  "Returns the Google matrix of the graph."
  {:db/ident :networkx.algorithms/google_matrix,
   :rdf/type :py/Function})

(def graph_clique_number
  "Returns the clique number of the graph."
  {:db/ident :networkx.algorithms/graph_clique_number,
   :rdf/type :py/Function})

(def graph_edit_distance
  "Returns GED (graph edit distance) between graphs G1 and G2."
  {:db/ident :networkx.algorithms/graph_edit_distance,
   :rdf/type :py/Function})

(def graph_hashing
  "Functions for hashing graphs to strings. Isomorphic graphs should be assigned identical hashes. For now, only Weisfeiler-Lehman hashing is implemented."
  {:db/ident :networkx.algorithms/graph_hashing,
   :rdf/type :py/Function})

(def graph_number_of_cliques
  "Returns the number of maximal cliques in the graph."
  {:db/ident :networkx.algorithms/graph_number_of_cliques,
   :rdf/type :py/Function})

(def graphical
  "Test sequences for graphiness."
  {:db/ident :networkx.algorithms/graphical,
   :rdf/type :py/Function})

(def greedy_color
  "Color a graph using various strategies of greedy graph coloring."
  {:db/ident :networkx.algorithms/greedy_color,
   :rdf/type :py/Function})

(def group
  "Group centrality measures."
  {:db/ident :networkx.algorithms/group,
   :rdf/type :py/Function})

(def group_betweenness_centrality
  "Compute the group betweenness centrality for a group of nodes."
  {:db/ident :networkx.algorithms/group_betweenness_centrality,
   :rdf/type :py/Function})

(def group_closeness_centrality
  "Compute the group closeness centrality for a group of nodes."
  {:db/ident :networkx.algorithms/group_closeness_centrality,
   :rdf/type :py/Function})

(def group_degree_centrality
  "Compute the group degree centrality for a group of nodes."
  {:db/ident :networkx.algorithms/group_degree_centrality,
   :rdf/type :py/Function})

(def group_in_degree_centrality
  "Compute the group in-degree centrality for a group of nodes."
  {:db/ident :networkx.algorithms/group_in_degree_centrality,
   :rdf/type :py/Function})

(def group_out_degree_centrality
  "Compute the group out-degree centrality for a group of nodes."
  {:db/ident :networkx.algorithms/group_out_degree_centrality,
   :rdf/type :py/Function})

(def harmonic
  "Functions for computing the harmonic centrality of a graph."
  {:db/ident :networkx.algorithms/harmonic,
   :rdf/type :py/Function})

(def harmonic_centrality
  "Compute harmonic centrality for nodes."
  {:db/ident :networkx.algorithms/harmonic_centrality,
   :rdf/type :py/Function})

(def has_bridges
  "Decide whether a graph has any bridges."
  {:db/ident :networkx.algorithms/has_bridges,
   :rdf/type :py/Function})

(def has_eulerian_path
  "Return True iff `G` has an Eulerian path."
  {:db/ident :networkx.algorithms/has_eulerian_path,
   :rdf/type :py/Function})

(def has_path
  "Returns *True* if *G* has a path from *source* to *target*."
  {:db/ident :networkx.algorithms/has_path,
   :rdf/type :py/Function})

(def hierarchy
  "Flow Hierarchy."
  {:db/ident :networkx.algorithms/hierarchy,
   :rdf/type :py/Function})

(def hits
  "Returns HITS hubs and authorities values for nodes."
  {:db/ident :networkx.algorithms/hits,
   :rdf/type :py/Function})

(def hits_alg
  "Hubs and authorities analysis of graph structure."
  {:db/ident :networkx.algorithms/hits_alg,
   :rdf/type :py/Function})

(def hybrid
  "Provides functions for finding and testing for locally `(k, l)`-connected graphs."
  {:db/ident :networkx.algorithms/hybrid,
   :rdf/type :py/Function})

(def identified_nodes
  "Returns the graph that results from contracting `u` and `v`."
  {:db/ident :networkx.algorithms/identified_nodes,
   :rdf/type :py/Function})

(def immediate_dominators
  "Returns the immediate dominators of all nodes of a directed graph."
  {:db/ident :networkx.algorithms/immediate_dominators,
   :rdf/type :py/Function})

(def in_degree_centrality
  "Compute the in-degree centrality for nodes."
  {:db/ident :networkx.algorithms/in_degree_centrality,
   :rdf/type :py/Function})

(def incremental_closeness_centrality
  "Incremental closeness centrality for nodes."
  {:db/ident :networkx.algorithms/incremental_closeness_centrality,
   :rdf/type :py/Function})

(def information_centrality
  "Compute current-flow closeness centrality for nodes."
  {:db/ident :networkx.algorithms/information_centrality,
   :rdf/type :py/Function})

(def intersection
  "Returns a new graph that contains only the nodes and the edges that exist in both G and H."
  {:db/ident :networkx.algorithms/intersection,
   :rdf/type :py/Function})

(def intersection_all
  "Returns a new graph that contains only the nodes and the edges that exist in all graphs."
  {:db/ident :networkx.algorithms/intersection_all,
   :rdf/type :py/Function})

(def intersection_array
  "Returns the intersection array of a distance-regular graph."
  {:db/ident :networkx.algorithms/intersection_array,
   :rdf/type :py/Function})

(def is_aperiodic
  "Returns True if `G` is aperiodic."
  {:db/ident :networkx.algorithms/is_aperiodic,
   :rdf/type :py/Function})

(def is_arborescence
  "Returns True if `G` is an arborescence."
  {:db/ident :networkx.algorithms/is_arborescence,
   :rdf/type :py/Function})

(def is_at_free
  "Check if a graph is AT-free."
  {:db/ident :networkx.algorithms/is_at_free,
   :rdf/type :py/Function})

(def is_attracting_component
  "Returns True if `G` consists of a single attracting component."
  {:db/ident :networkx.algorithms/is_attracting_component,
   :rdf/type :py/Function})

(def is_biconnected
  "Returns True if the graph is biconnected, False otherwise."
  {:db/ident :networkx.algorithms/is_biconnected,
   :rdf/type :py/Function})

(def is_bipartite
  "Returns True if graph G is bipartite, False if not."
  {:db/ident :networkx.algorithms/is_bipartite,
   :rdf/type :py/Function})

(def is_branching
  "Returns True if `G` is a branching."
  {:db/ident :networkx.algorithms/is_branching,
   :rdf/type :py/Function})

(def is_chordal
  "Checks whether G is a chordal graph."
  {:db/ident :networkx.algorithms/is_chordal,
   :rdf/type :py/Function})

(def is_connected
  "Returns True if the graph is connected, False otherwise."
  {:db/ident :networkx.algorithms/is_connected,
   :rdf/type :py/Function})

(def is_digraphical
  "Returns True if some directed graph can realize the in- and out-degree sequences."
  {:db/ident :networkx.algorithms/is_digraphical,
   :rdf/type :py/Function})

(def is_directed_acyclic_graph
  "Returns True if the graph `G` is a directed acyclic graph (DAG) or False if not."
  {:db/ident :networkx.algorithms/is_directed_acyclic_graph,
   :rdf/type :py/Function})

(def is_distance_regular
  "Returns True if the graph is distance regular, False otherwise."
  {:db/ident :networkx.algorithms/is_distance_regular,
   :rdf/type :py/Function})

(def is_dominating_set
  "Checks if `nbunch` is a dominating set for `G`."
  {:db/ident :networkx.algorithms/is_dominating_set,
   :rdf/type :py/Function})

(def is_edge_cover
  "Decides whether a set of edges is a valid edge cover of the graph."
  {:db/ident :networkx.algorithms/is_edge_cover,
   :rdf/type :py/Function})

(def is_eulerian
  "Returns True if and only if `G` is Eulerian."
  {:db/ident :networkx.algorithms/is_eulerian,
   :rdf/type :py/Function})

(def is_forest
  "Returns True if `G` is a forest."
  {:db/ident :networkx.algorithms/is_forest,
   :rdf/type :py/Function})

(def is_graphical
  "Returns True if sequence is a valid degree sequence."
  {:db/ident :networkx.algorithms/is_graphical,
   :rdf/type :py/Function})

(def is_isolate
  "Determines whether a node is an isolate."
  {:db/ident :networkx.algorithms/is_isolate,
   :rdf/type :py/Function})

(def is_isomorphic
  "Returns True if the graphs G1 and G2 are isomorphic and False otherwise."
  {:db/ident :networkx.algorithms/is_isomorphic,
   :rdf/type :py/Function})

(def is_k_edge_connected
  "Tests to see if a graph is k-edge-connected."
  {:db/ident :networkx.algorithms/is_k_edge_connected,
   :rdf/type :py/Function})

(def is_k_regular
  "Determines whether the graph ``G`` is a k-regular graph."
  {:db/ident :networkx.algorithms/is_k_regular,
   :rdf/type :py/Function})

(def is_kl_connected
  "Returns True if and only if `G` is locally `(k, l)`-connected."
  {:db/ident :networkx.algorithms/is_kl_connected,
   :rdf/type :py/Function})

(def is_matching
  "Return True if ``matching`` is a valid matching of ``G``"
  {:db/ident :networkx.algorithms/is_matching,
   :rdf/type :py/Function})

(def is_maximal_matching
  "Return True if ``matching`` is a maximal matching of ``G``"
  {:db/ident :networkx.algorithms/is_maximal_matching,
   :rdf/type :py/Function})

(def is_minimal_d_separator
  "Determine if a d-separating set is minimal."
  {:db/ident :networkx.algorithms/is_minimal_d_separator,
   :rdf/type :py/Function})

(def is_multigraphical
  "Returns True if some multigraph can realize the sequence."
  {:db/ident :networkx.algorithms/is_multigraphical,
   :rdf/type :py/Function})

(def is_perfect_matching
  "Return True if ``matching`` is a perfect matching for ``G``"
  {:db/ident :networkx.algorithms/is_perfect_matching,
   :rdf/type :py/Function})

(def is_planar
  "Returns True if and only if `G` is planar."
  {:db/ident :networkx.algorithms/is_planar,
   :rdf/type :py/Function})

(def is_pseudographical
  "Returns True if some pseudograph can realize the sequence."
  {:db/ident :networkx.algorithms/is_pseudographical,
   :rdf/type :py/Function})

(def is_regular
  "Determines whether the graph ``G`` is a regular graph."
  {:db/ident :networkx.algorithms/is_regular,
   :rdf/type :py/Function})

(def is_semiconnected
  "Returns True if the graph is semiconnected, False otherwise."
  {:db/ident :networkx.algorithms/is_semiconnected,
   :rdf/type :py/Function})

(def is_semieulerian
  "Return True iff `G` is semi-Eulerian."
  {:db/ident :networkx.algorithms/is_semieulerian,
   :rdf/type :py/Function})

(def is_simple_path
  "Returns True if and only if `nodes` form a simple path in `G`."
  {:db/ident :networkx.algorithms/is_simple_path,
   :rdf/type :py/Function})

(def is_strongly_connected
  "Test directed graph for strong connectivity."
  {:db/ident :networkx.algorithms/is_strongly_connected,
   :rdf/type :py/Function})

(def is_strongly_regular
  "Returns True if and only if the given graph is strongly regular."
  {:db/ident :networkx.algorithms/is_strongly_regular,
   :rdf/type :py/Function})

(def is_tree
  "Returns True if `G` is a tree."
  {:db/ident :networkx.algorithms/is_tree,
   :rdf/type :py/Function})

(def is_triad
  "Returns True if the graph G is a triad, else False."
  {:db/ident :networkx.algorithms/is_triad,
   :rdf/type :py/Function})

(def is_valid_degree_sequence_erdos_gallai
  "Returns True if deg_sequence can be realized by a simple graph."
  {:db/ident :networkx.algorithms/is_valid_degree_sequence_erdos_gallai,
   :rdf/type :py/Function})

(def is_valid_degree_sequence_havel_hakimi
  "Returns True if deg_sequence can be realized by a simple graph."
  {:db/ident :networkx.algorithms/is_valid_degree_sequence_havel_hakimi,
   :rdf/type :py/Function})

(def is_weakly_connected
  "Test directed graph for weak connectivity."
  {:db/ident :networkx.algorithms/is_weakly_connected,
   :rdf/type :py/Function})

(def isolate
  "Functions for identifying isolate (degree zero) nodes."
  {:db/ident :networkx.algorithms/isolate,
   :rdf/type :py/Function})

(def isolates
  "Iterator over isolates in the graph."
  {:db/ident :networkx.algorithms/isolates,
   :rdf/type :py/Function})

(def isomorphism
  "isomorphism"
  {:db/ident :networkx.algorithms/isomorphism,
   :rdf/type :py/Function})

(def jaccard_coefficient
  "Compute the Jaccard coefficient of all node pairs in ebunch."
  {:db/ident :networkx.algorithms/jaccard_coefficient,
   :rdf/type :py/Function})

(def johnson
  "Uses Johnson's Algorithm to compute shortest paths."
  {:db/ident :networkx.algorithms/johnson,
   :rdf/type :py/Function})

(def join
  "Returns a new rooted tree with a root node joined with the roots of each of the given rooted trees."
  {:db/ident :networkx.algorithms/join,
   :rdf/type :py/Function})

(def junction_tree
  "Returns a junction tree of a given graph."
  {:db/ident :networkx.algorithms/junction_tree,
   :rdf/type :py/Function})

(def k_components
  "Returns the k-component structure of a graph G."
  {:db/ident :networkx.algorithms/k_components,
   :rdf/type :py/Function})

(def k_core
  "Returns the k-core of G."
  {:db/ident :networkx.algorithms/k_core,
   :rdf/type :py/Function})

(def k_corona
  "Returns the k-corona of G."
  {:db/ident :networkx.algorithms/k_corona,
   :rdf/type :py/Function})

(def k_crust
  "Returns the k-crust of G."
  {:db/ident :networkx.algorithms/k_crust,
   :rdf/type :py/Function})

(def k_edge_augmentation
  "Finds set of edges to k-edge-connect G."
  {:db/ident :networkx.algorithms/k_edge_augmentation,
   :rdf/type :py/Function})

(def k_edge_components
  "Generates nodes in each maximal k-edge-connected component in G."
  {:db/ident :networkx.algorithms/k_edge_components,
   :rdf/type :py/Function})

(def k_edge_subgraphs
  "Generates nodes in each maximal k-edge-connected subgraph in G."
  {:db/ident :networkx.algorithms/k_edge_subgraphs,
   :rdf/type :py/Function})

(def k_factor
  "Compute a k-factor of G"
  {:db/ident :networkx.algorithms/k_factor,
   :rdf/type :py/Function})

(def k_shell
  "Returns the k-shell of G."
  {:db/ident :networkx.algorithms/k_shell,
   :rdf/type :py/Function})

(def k_truss
  "Returns the k-truss of `G`."
  {:db/ident :networkx.algorithms/k_truss,
   :rdf/type :py/Function})

(def katz
  "Katz centrality."
  {:db/ident :networkx.algorithms/katz,
   :rdf/type :py/Function})

(def katz_centrality
  "Compute the Katz centrality for the nodes of the graph G."
  {:db/ident :networkx.algorithms/katz_centrality,
   :rdf/type :py/Function})

(def katz_centrality_numpy
  "Compute the Katz centrality for the graph G."
  {:db/ident :networkx.algorithms/katz_centrality_numpy,
   :rdf/type :py/Function})

(def kl_connected_subgraph
  "Returns the maximum locally `(k, l)`-connected subgraph of `G`."
  {:db/ident :networkx.algorithms/kl_connected_subgraph,
   :rdf/type :py/Function})

(def kosaraju_strongly_connected_components
  "Generate nodes in strongly connected components of graph."
  {:db/ident :networkx.algorithms/kosaraju_strongly_connected_components,
   :rdf/type :py/Function})

(def lattice_reference
  "Latticize the given graph by swapping edges."
  {:db/ident :networkx.algorithms/lattice_reference,
   :rdf/type :py/Function})

(def lexicographic_product
  "Returns the lexicographic product of G and H."
  {:db/ident :networkx.algorithms/lexicographic_product,
   :rdf/type :py/Function})

(def lexicographical_topological_sort
  "Generate the nodes in the unique lexicographical topological sort order."
  {:db/ident :networkx.algorithms/lexicographical_topological_sort,
   :rdf/type :py/Function})

(def link_analysis
  "link_analysis"
  {:db/ident :networkx.algorithms/link_analysis,
   :rdf/type :py/Function})

(def link_prediction
  "Link prediction algorithms."
  {:db/ident :networkx.algorithms/link_prediction,
   :rdf/type :py/Function})

(def load
  "Load centrality."
  {:db/ident :networkx.algorithms/load,
   :rdf/type :py/Function})

(def load_centrality
  "Compute load centrality for nodes."
  {:db/ident :networkx.algorithms/load_centrality,
   :rdf/type :py/Function})

(def local_bridges
  "Iterate over local bridges of `G` optionally computing the span"
  {:db/ident :networkx.algorithms/local_bridges,
   :rdf/type :py/Function})

(def local_constraint
  "Returns the local constraint on the node ``u`` with respect to the node ``v`` in the graph ``G``."
  {:db/ident :networkx.algorithms/local_constraint,
   :rdf/type :py/Function})

(def local_efficiency
  "Returns the average local efficiency of the graph."
  {:db/ident :networkx.algorithms/local_efficiency,
   :rdf/type :py/Function})

(def local_reaching_centrality
  "Returns the local reaching centrality of a node in a directed graph."
  {:db/ident :networkx.algorithms/local_reaching_centrality,
   :rdf/type :py/Function})

(def lowest_common_ancestor
  "Compute the lowest common ancestor of the given pair of nodes."
  {:db/ident :networkx.algorithms/lowest_common_ancestor,
   :rdf/type :py/Function})

(def lowest_common_ancestors
  "Algorithms for finding the lowest common ancestor of trees and DAGs."
  {:db/ident :networkx.algorithms/lowest_common_ancestors,
   :rdf/type :py/Function})

(def make_clique_bipartite
  "Returns the bipartite clique graph corresponding to `G`."
  {:db/ident :networkx.algorithms/make_clique_bipartite,
   :rdf/type :py/Function})

(def make_max_clique_graph
  "Returns the maximal clique graph of the given graph."
  {:db/ident :networkx.algorithms/make_max_clique_graph,
   :rdf/type :py/Function})

(def matching
  "Functions for computing and verifying matchings in a graph."
  {:db/ident :networkx.algorithms/matching,
   :rdf/type :py/Function})

(def max_flow_min_cost
  "Returns a maximum (s, t)-flow of minimum cost."
  {:db/ident :networkx.algorithms/max_flow_min_cost,
   :rdf/type :py/Function})

(def max_weight_clique
  "Find a maximum weight clique in G."
  {:db/ident :networkx.algorithms/max_weight_clique,
   :rdf/type :py/Function})

(def max_weight_matching
  "Compute a maximum-weighted matching of G."
  {:db/ident :networkx.algorithms/max_weight_matching,
   :rdf/type :py/Function})

(def maximal_independent_set
  "Returns a random maximal independent set guaranteed to contain a given set of nodes."
  {:db/ident :networkx.algorithms/maximal_independent_set,
   :rdf/type :py/Function})

(def maximal_matching
  "Find a maximal matching in the graph."
  {:db/ident :networkx.algorithms/maximal_matching,
   :rdf/type :py/Function})

(def maximum_branching
  "Returns a maximum branching from G."
  {:db/ident :networkx.algorithms/maximum_branching,
   :rdf/type :py/Function})

(def maximum_flow
  "Find a maximum single-commodity flow."
  {:db/ident :networkx.algorithms/maximum_flow,
   :rdf/type :py/Function})

(def maximum_flow_value
  "Find the value of maximum single-commodity flow."
  {:db/ident :networkx.algorithms/maximum_flow_value,
   :rdf/type :py/Function})

(def maximum_spanning_arborescence
  "Returns a maximum spanning arborescence from G."
  {:db/ident :networkx.algorithms/maximum_spanning_arborescence,
   :rdf/type :py/Function})

(def maximum_spanning_edges
  "Generate edges in a maximum spanning forest of an undirected weighted graph."
  {:db/ident :networkx.algorithms/maximum_spanning_edges,
   :rdf/type :py/Function})

(def maximum_spanning_tree
  "Returns a maximum spanning tree or forest on an undirected graph `G`."
  {:db/ident :networkx.algorithms/maximum_spanning_tree,
   :rdf/type :py/Function})

(def min_cost_flow
  "Returns a minimum cost flow satisfying all demands in digraph G."
  {:db/ident :networkx.algorithms/min_cost_flow,
   :rdf/type :py/Function})

(def min_cost_flow_cost
  "Find the cost of a minimum cost flow satisfying all demands in digraph G."
  {:db/ident :networkx.algorithms/min_cost_flow_cost,
   :rdf/type :py/Function})

(def min_edge_cover
  "Returns the min cardinality edge cover of the graph as a set of edges."
  {:db/ident :networkx.algorithms/min_edge_cover,
   :rdf/type :py/Function})

(def min_weight_matching
  "Computing a minimum-weight maximal matching of G."
  {:db/ident :networkx.algorithms/min_weight_matching,
   :rdf/type :py/Function})

(def minimal_d_separator
  "Compute a minimal d-separating set between 'u' and 'v'."
  {:db/ident :networkx.algorithms/minimal_d_separator,
   :rdf/type :py/Function})

(def minimum_branching
  "Returns a minimum branching from G."
  {:db/ident :networkx.algorithms/minimum_branching,
   :rdf/type :py/Function})

(def minimum_cut
  "Compute the value and the node partition of a minimum (s, t)-cut."
  {:db/ident :networkx.algorithms/minimum_cut,
   :rdf/type :py/Function})

(def minimum_cut_value
  "Compute the value of a minimum (s, t)-cut."
  {:db/ident :networkx.algorithms/minimum_cut_value,
   :rdf/type :py/Function})

(def minimum_cycle_basis
  "Returns a minimum weight cycle basis for G"
  {:db/ident :networkx.algorithms/minimum_cycle_basis,
   :rdf/type :py/Function})

(def minimum_edge_cut
  "Returns a set of edges of minimum cardinality that disconnects G."
  {:db/ident :networkx.algorithms/minimum_edge_cut,
   :rdf/type :py/Function})

(def minimum_node_cut
  "Returns a set of nodes of minimum cardinality that disconnects G."
  {:db/ident :networkx.algorithms/minimum_node_cut,
   :rdf/type :py/Function})

(def minimum_spanning_arborescence
  "Returns a minimum spanning arborescence from G."
  {:db/ident :networkx.algorithms/minimum_spanning_arborescence,
   :rdf/type :py/Function})

(def minimum_spanning_edges
  "Generate edges in a minimum spanning forest of an undirected weighted graph."
  {:db/ident :networkx.algorithms/minimum_spanning_edges,
   :rdf/type :py/Function})

(def minimum_spanning_tree
  "Returns a minimum spanning tree or forest on an undirected graph `G`."
  {:db/ident :networkx.algorithms/minimum_spanning_tree,
   :rdf/type :py/Function})

(def minors
  "Subpackages related to graph-minor problems."
  {:db/ident :networkx.algorithms/minors,
   :rdf/type :py/Function})

(def mis
  "Algorithm to find a maximal (not maximum) independent set."
  {:db/ident :networkx.algorithms/mis,
   :rdf/type :py/Function})

(def mixing
  "Mixing matrices for node attributes and degree."
  {:db/ident :networkx.algorithms/mixing,
   :rdf/type :py/Function})

(def mixing_dict
  "Returns a dictionary representation of mixing matrix."
  {:db/ident :networkx.algorithms/mixing_dict,
   :rdf/type :py/Function})

(def mixing_expansion
  "Returns the mixing expansion between two node sets."
  {:db/ident :networkx.algorithms/mixing_expansion,
   :rdf/type :py/Function})

(def moral
  "Function for computing the moral graph of a directed graph."
  {:db/ident :networkx.algorithms/moral,
   :rdf/type :py/Function})

(def moral_graph
  "Return the Moral Graph"
  {:db/ident :networkx.algorithms/moral_graph,
   :rdf/type :py/Function})

(def multi_source_dijkstra
  "Find shortest weighted paths and lengths from a given set of source nodes."
  {:db/ident :networkx.algorithms/multi_source_dijkstra,
   :rdf/type :py/Function})

(def multi_source_dijkstra_path
  "Find shortest weighted paths in G from a given set of source nodes."
  {:db/ident :networkx.algorithms/multi_source_dijkstra_path,
   :rdf/type :py/Function})

(def multi_source_dijkstra_path_length
  "Find shortest weighted path lengths in G from a given set of source nodes."
  {:db/ident :networkx.algorithms/multi_source_dijkstra_path_length,
   :rdf/type :py/Function})

(def negative_edge_cycle
  "Returns True if there exists a negative edge cycle anywhere in G."
  {:db/ident :networkx.algorithms/negative_edge_cycle,
   :rdf/type :py/Function})

(def neighbor_degree
  "neighbor_degree"
  {:db/ident :networkx.algorithms/neighbor_degree,
   :rdf/type :py/Function})

(def network_simplex
  "Find a minimum cost flow satisfying all demands in digraph G."
  {:db/ident :networkx.algorithms/network_simplex,
   :rdf/type :py/Function})

(def node_attribute_xy
  "Returns iterator of node-attribute pairs for all edges in G."
  {:db/ident :networkx.algorithms/node_attribute_xy,
   :rdf/type :py/Function})

(def node_boundary
  "Returns the node boundary of `nbunch1`."
  {:db/ident :networkx.algorithms/node_boundary,
   :rdf/type :py/Function})

(def node_classification
  "This module provides the functions for node classification problem."
  {:db/ident :networkx.algorithms/node_classification,
   :rdf/type :py/Function})

(def node_clique_number
  "Returns the size of the largest maximal clique containing each given node."
  {:db/ident :networkx.algorithms/node_clique_number,
   :rdf/type :py/Function})

(def node_connected_component
  "Returns the set of nodes in the component of graph containing node n."
  {:db/ident :networkx.algorithms/node_connected_component,
   :rdf/type :py/Function})

(def node_connectivity
  "Returns node connectivity for a graph or digraph G."
  {:db/ident :networkx.algorithms/node_connectivity,
   :rdf/type :py/Function})

(def node_degree_xy
  "Generate node degree-degree pairs for edges in G."
  {:db/ident :networkx.algorithms/node_degree_xy,
   :rdf/type :py/Function})

(def node_disjoint_paths
  "Computes node disjoint paths between source and target."
  {:db/ident :networkx.algorithms/node_disjoint_paths,
   :rdf/type :py/Function})

(def node_expansion
  "Returns the node expansion of the set `S`."
  {:db/ident :networkx.algorithms/node_expansion,
   :rdf/type :py/Function})

(def non_randomness
  "Compute the non-randomness of graph G."
  {:db/ident :networkx.algorithms/non_randomness,
   :rdf/type :py/Function})

(def normalized_cut_size
  "Returns the normalized size of the cut between two sets of nodes."
  {:db/ident :networkx.algorithms/normalized_cut_size,
   :rdf/type :py/Function})

(def number_attracting_components
  "Returns the number of attracting components in `G`."
  {:db/ident :networkx.algorithms/number_attracting_components,
   :rdf/type :py/Function})

(def number_connected_components
  "Returns the number of connected components."
  {:db/ident :networkx.algorithms/number_connected_components,
   :rdf/type :py/Function})

(def number_of_cliques
  "Returns the number of maximal cliques for each node."
  {:db/ident :networkx.algorithms/number_of_cliques,
   :rdf/type :py/Function})

(def number_of_isolates
  "Returns the number of isolates in the graph."
  {:db/ident :networkx.algorithms/number_of_isolates,
   :rdf/type :py/Function})

(def number_strongly_connected_components
  "Returns number of strongly connected components in graph."
  {:db/ident :networkx.algorithms/number_strongly_connected_components,
   :rdf/type :py/Function})

(def number_weakly_connected_components
  "Returns the number of weakly connected components in G."
  {:db/ident :networkx.algorithms/number_weakly_connected_components,
   :rdf/type :py/Function})

(def numeric_assortativity_coefficient
  "Compute assortativity for numerical node attributes."
  {:db/ident :networkx.algorithms/numeric_assortativity_coefficient,
   :rdf/type :py/Function})

(def omega
  "Returns the small-world coefficient (omega) of a graph"
  {:db/ident :networkx.algorithms/omega,
   :rdf/type :py/Function})

(def onion_layers
  "Returns the layer of each vertex in an onion decomposition of the graph."
  {:db/ident :networkx.algorithms/onion_layers,
   :rdf/type :py/Function})

(def operators
  "operators"
  {:db/ident :networkx.algorithms/operators,
   :rdf/type :py/Function})

(def optimal_edit_paths
  "Returns all minimum-cost edit paths transforming G1 to G2."
  {:db/ident :networkx.algorithms/optimal_edit_paths,
   :rdf/type :py/Function})

(def optimize_edit_paths
  "GED (graph edit distance) calculation: advanced interface."
  {:db/ident :networkx.algorithms/optimize_edit_paths,
   :rdf/type :py/Function})

(def optimize_graph_edit_distance
  "Returns consecutive approximations of GED (graph edit distance) between graphs G1 and G2."
  {:db/ident :networkx.algorithms/optimize_graph_edit_distance,
   :rdf/type :py/Function})

(def out_degree_centrality
  "Compute the out-degree centrality for nodes."
  {:db/ident :networkx.algorithms/out_degree_centrality,
   :rdf/type :py/Function})

(def overall_reciprocity
  "Compute the reciprocity for the whole graph."
  {:db/ident :networkx.algorithms/overall_reciprocity,
   :rdf/type :py/Function})

(def pagerank
  "Returns the PageRank of the nodes in the graph."
  {:db/ident :networkx.algorithms/pagerank,
   :rdf/type :py/Function})

(def pagerank_alg
  "PageRank analysis of graph structure."
  {:db/ident :networkx.algorithms/pagerank_alg,
   :rdf/type :py/Function})

(def pairs
  "Generators of  x-y pairs of node data."
  {:db/ident :networkx.algorithms/pairs,
   :rdf/type :py/Function})

(def panther_similarity
  "Returns the Panther similarity of nodes in the graph `G` to node ``v``."
  {:db/ident :networkx.algorithms/panther_similarity,
   :rdf/type :py/Function})

(def partition_spanning_tree
  "Find a spanning tree while respecting a partition of edges."
  {:db/ident :networkx.algorithms/partition_spanning_tree,
   :rdf/type :py/Function})

(def percolation
  "Percolation centrality measures."
  {:db/ident :networkx.algorithms/percolation,
   :rdf/type :py/Function})

(def percolation_centrality
  "Compute the percolation centrality for nodes."
  {:db/ident :networkx.algorithms/percolation_centrality,
   :rdf/type :py/Function})

(def periphery
  "Returns the periphery of the graph G."
  {:db/ident :networkx.algorithms/periphery,
   :rdf/type :py/Function})

(def planar_drawing
  "planar_drawing"
  {:db/ident :networkx.algorithms/planar_drawing,
   :rdf/type :py/Function})

(def planarity
  "planarity"
  {:db/ident :networkx.algorithms/planarity,
   :rdf/type :py/Function})

(def polynomials
  "Provides algorithms supporting the computation of graph polynomials."
  {:db/ident :networkx.algorithms/polynomials,
   :rdf/type :py/Function})

(def power
  "Returns the specified power of a graph."
  {:db/ident :networkx.algorithms/power,
   :rdf/type :py/Function})

(def predecessor
  "Returns dict of predecessors for the path from source to all nodes in G."
  {:db/ident :networkx.algorithms/predecessor,
   :rdf/type :py/Function})

(def preferential_attachment
  "Compute the preferential attachment score of all node pairs in ebunch."
  {:db/ident :networkx.algorithms/preferential_attachment,
   :rdf/type :py/Function})

(def product
  "Graph products."
  {:db/ident :networkx.algorithms/product,
   :rdf/type :py/Function})

(def projected_graph
  "Returns the projection of B onto one of its node sets."
  {:db/ident :networkx.algorithms/projected_graph,
   :rdf/type :py/Function})

(def prominent_group
  "Find the prominent group of size $k$ in graph $G$. The prominence of the group is evaluated by the group betweenness centrality."
  {:db/ident :networkx.algorithms/prominent_group,
   :rdf/type :py/Function})

(def quotient_graph
  "Returns the quotient graph of `G` under the specified equivalence relation on nodes."
  {:db/ident :networkx.algorithms/quotient_graph,
   :rdf/type :py/Function})

(def ra_index_soundarajan_hopcroft
  "Compute the resource allocation index of all node pairs in ebunch using community information."
  {:db/ident :networkx.algorithms/ra_index_soundarajan_hopcroft,
   :rdf/type :py/Function})

(def radius
  "Returns the radius of the graph G."
  {:db/ident :networkx.algorithms/radius,
   :rdf/type :py/Function})

(def random_reference
  "Compute a random graph by swapping edges of a given graph."
  {:db/ident :networkx.algorithms/random_reference,
   :rdf/type :py/Function})

(def random_spanning_tree
  "Sample a random spanning tree using the edges weights of `G`."
  {:db/ident :networkx.algorithms/random_spanning_tree,
   :rdf/type :py/Function})

(def random_triad
  "Returns a random triad from a directed graph."
  {:db/ident :networkx.algorithms/random_triad,
   :rdf/type :py/Function})

(def reaching
  "Functions for computing reaching centrality of a node or a graph."
  {:db/ident :networkx.algorithms/reaching,
   :rdf/type :py/Function})

(def reciprocity
  "Compute the reciprocity in a directed graph."
  {:db/ident :networkx.algorithms/reciprocity,
   :rdf/type :py/Function})

(def reconstruct_path
  "Reconstruct a path from source to target using the predecessors dict as returned by floyd_warshall_predecessor_and_distance"
  {:db/ident :networkx.algorithms/reconstruct_path,
   :rdf/type :py/Function})

(def recursive_simple_cycles
  "Find simple cycles (elementary circuits) of a directed graph."
  {:db/ident :networkx.algorithms/recursive_simple_cycles,
   :rdf/type :py/Function})

(def regular
  "Functions for computing and verifying regular graphs."
  {:db/ident :networkx.algorithms/regular,
   :rdf/type :py/Function})

(def resistance_distance
  "Returns the resistance distance between node A and node B on graph G."
  {:db/ident :networkx.algorithms/resistance_distance,
   :rdf/type :py/Function})

(def resource_allocation_index
  "Compute the resource allocation index of all node pairs in ebunch."
  {:db/ident :networkx.algorithms/resource_allocation_index,
   :rdf/type :py/Function})

(def reverse
  "Returns the reverse directed graph of G."
  {:db/ident :networkx.algorithms/reverse,
   :rdf/type :py/Function})

(def rich_club_coefficient
  "Returns the rich-club coefficient of the graph `G`."
  {:db/ident :networkx.algorithms/rich_club_coefficient,
   :rdf/type :py/Function})

(def richclub
  "Functions for computing rich-club coefficients."
  {:db/ident :networkx.algorithms/richclub,
   :rdf/type :py/Function})

(def rooted_product
  "Return the rooted product of graphs G and H rooted at root in H."
  {:db/ident :networkx.algorithms/rooted_product,
   :rdf/type :py/Function})

(def s_metric
  "Returns the s-metric of graph."
  {:db/ident :networkx.algorithms/s_metric,
   :rdf/type :py/Function})

(def second_order
  "Copyright (c) 2015 – Thomson Licensing, SAS"
  {:db/ident :networkx.algorithms/second_order,
   :rdf/type :py/Function})

(def second_order_centrality
  "Compute the second order centrality for nodes of G."
  {:db/ident :networkx.algorithms/second_order_centrality,
   :rdf/type :py/Function})

(def semiconnected
  "Semiconnectedness."
  {:db/ident :networkx.algorithms/semiconnected,
   :rdf/type :py/Function})

(def shortest_path
  "Compute shortest paths in the graph."
  {:db/ident :networkx.algorithms/shortest_path,
   :rdf/type :py/Function})

(def shortest_path_length
  "Compute shortest path lengths in the graph."
  {:db/ident :networkx.algorithms/shortest_path_length,
   :rdf/type :py/Function})

(def shortest_paths
  "shortest_paths"
  {:db/ident :networkx.algorithms/shortest_paths,
   :rdf/type :py/Function})

(def shortest_simple_paths
  "Generate all simple paths in the graph G from source to target, starting from shortest ones."
  {:db/ident :networkx.algorithms/shortest_simple_paths,
   :rdf/type :py/Function})

(def sigma
  "Returns the small-world coefficient (sigma) of the given graph."
  {:db/ident :networkx.algorithms/sigma,
   :rdf/type :py/Function})

(def similarity
  "Functions measuring similarity using graph edit distance."
  {:db/ident :networkx.algorithms/similarity,
   :rdf/type :py/Function})

(def simple_cycles
  "Find simple cycles (elementary circuits) of a directed graph."
  {:db/ident :networkx.algorithms/simple_cycles,
   :rdf/type :py/Function})

(def simple_paths
  "simple_paths"
  {:db/ident :networkx.algorithms/simple_paths,
   :rdf/type :py/Function})

(def simrank_similarity
  "Returns the SimRank similarity of nodes in the graph ``G``."
  {:db/ident :networkx.algorithms/simrank_similarity,
   :rdf/type :py/Function})

(def single_source_bellman_ford
  "Compute shortest paths and lengths in a weighted graph G."
  {:db/ident :networkx.algorithms/single_source_bellman_ford,
   :rdf/type :py/Function})

(def single_source_bellman_ford_path
  "Compute shortest path between source and all other reachable nodes for a weighted graph."
  {:db/ident :networkx.algorithms/single_source_bellman_ford_path,
   :rdf/type :py/Function})

(def single_source_bellman_ford_path_length
  "Compute the shortest path length between source and all other reachable nodes for a weighted graph."
  {:db/ident :networkx.algorithms/single_source_bellman_ford_path_length,
   :rdf/type :py/Function})

(def single_source_dijkstra
  "Find shortest weighted paths and lengths from a source node."
  {:db/ident :networkx.algorithms/single_source_dijkstra,
   :rdf/type :py/Function})

(def single_source_dijkstra_path
  "Find shortest weighted paths in G from a source node."
  {:db/ident :networkx.algorithms/single_source_dijkstra_path,
   :rdf/type :py/Function})

(def single_source_dijkstra_path_length
  "Find shortest weighted path lengths in G from a source node."
  {:db/ident :networkx.algorithms/single_source_dijkstra_path_length,
   :rdf/type :py/Function})

(def single_source_shortest_path
  "Compute shortest path between source and all other nodes reachable from source."
  {:db/ident :networkx.algorithms/single_source_shortest_path,
   :rdf/type :py/Function})

(def single_source_shortest_path_length
  "Compute the shortest path lengths from source to all reachable nodes."
  {:db/ident :networkx.algorithms/single_source_shortest_path_length,
   :rdf/type :py/Function})

(def single_target_shortest_path
  "Compute shortest path to target from all nodes that reach target."
  {:db/ident :networkx.algorithms/single_target_shortest_path,
   :rdf/type :py/Function})

(def single_target_shortest_path_length
  "Compute the shortest path lengths to target from all reachable nodes."
  {:db/ident :networkx.algorithms/single_target_shortest_path_length,
   :rdf/type :py/Function})

(def smallworld
  "Functions for estimating the small-world-ness of graphs."
  {:db/ident :networkx.algorithms/smallworld,
   :rdf/type :py/Function})

(def smetric
  "smetric"
  {:db/ident :networkx.algorithms/smetric,
   :rdf/type :py/Function})

(def snap_aggregation
  "Creates a summary graph based on attributes and connectivity."
  {:db/ident :networkx.algorithms/snap_aggregation,
   :rdf/type :py/Function})

(def spanner
  "Returns a spanner of the given graph with the given stretch."
  {:db/ident :networkx.algorithms/spanner,
   :rdf/type :py/Function})

(def sparsifiers
  "Functions for computing sparsifiers of graphs."
  {:db/ident :networkx.algorithms/sparsifiers,
   :rdf/type :py/Function})

(def square_clustering
  "Compute the squares clustering coefficient for nodes."
  {:db/ident :networkx.algorithms/square_clustering,
   :rdf/type :py/Function})

(def stoer_wagner
  "Returns the weighted minimum edge cut using the Stoer-Wagner algorithm."
  {:db/ident :networkx.algorithms/stoer_wagner,
   :rdf/type :py/Function})

(def strong_product
  "Returns the strong product of G and H."
  {:db/ident :networkx.algorithms/strong_product,
   :rdf/type :py/Function})

(def strongly_connected
  "Strongly connected components."
  {:db/ident :networkx.algorithms/strongly_connected,
   :rdf/type :py/Function})

(def strongly_connected_components
  "Generate nodes in strongly connected components of graph."
  {:db/ident :networkx.algorithms/strongly_connected_components,
   :rdf/type :py/Function})

(def strongly_connected_components_recursive
  "Generate nodes in strongly connected components of graph."
  {:db/ident :networkx.algorithms/strongly_connected_components_recursive,
   :rdf/type :py/Function})

(def structuralholes
  "Functions for computing measures of structural holes."
  {:db/ident :networkx.algorithms/structuralholes,
   :rdf/type :py/Function})

(def subgraph_alg
  "Subraph centrality and communicability betweenness."
  {:db/ident :networkx.algorithms/subgraph_alg,
   :rdf/type :py/Function})

(def subgraph_centrality
  "Returns subgraph centrality for each node in G."
  {:db/ident :networkx.algorithms/subgraph_centrality,
   :rdf/type :py/Function})

(def subgraph_centrality_exp
  "Returns the subgraph centrality for each node of G."
  {:db/ident :networkx.algorithms/subgraph_centrality_exp,
   :rdf/type :py/Function})

(def summarization
  "Graph summarization finds smaller representations of graphs resulting in faster runtime of algorithms, reduced storage needs, and noise reduction. Summarization has applications in areas such as visualization, pattern mining, clustering and community detection, and more.  Core graph summarization techniques are grouping/aggregation, bit-compression, simplification/sparsification, and influence based. Graph summarization algorithms often produce either summary graphs in the form of supergraphs or sparsified graphs, or a list of independent structures. Supergraphs are the most common product, which consist of supernodes and original nodes and are connected by edges and superedges, which represent aggregate edges between nodes and supernodes."
  {:db/ident :networkx.algorithms/summarization,
   :rdf/type :py/Function})

(def swap
  "Swap edges in a graph."
  {:db/ident :networkx.algorithms/swap,
   :rdf/type :py/Function})

(def symmetric_difference
  "Returns new graph with edges that exist in either G or H but not both."
  {:db/ident :networkx.algorithms/symmetric_difference,
   :rdf/type :py/Function})

(def tensor_product
  "Returns the tensor product of G and H."
  {:db/ident :networkx.algorithms/tensor_product,
   :rdf/type :py/Function})

(def to_nested_tuple
  "Returns a nested tuple representation of the given tree."
  {:db/ident :networkx.algorithms/to_nested_tuple,
   :rdf/type :py/Function})

(def to_prufer_sequence
  "Returns the Prüfer sequence of the given tree."
  {:db/ident :networkx.algorithms/to_prufer_sequence,
   :rdf/type :py/Function})

(def topological_generations
  "Stratifies a DAG into generations."
  {:db/ident :networkx.algorithms/topological_generations,
   :rdf/type :py/Function})

(def topological_sort
  "Returns a generator of nodes in topologically sorted order."
  {:db/ident :networkx.algorithms/topological_sort,
   :rdf/type :py/Function})

(def tournament
  "Functions concerning tournament graphs."
  {:db/ident :networkx.algorithms/tournament,
   :rdf/type :py/Function})

(def transitive_closure
  "Returns transitive closure of a graph"
  {:db/ident :networkx.algorithms/transitive_closure,
   :rdf/type :py/Function})

(def transitive_closure_dag
  "Returns the transitive closure of a directed acyclic graph."
  {:db/ident :networkx.algorithms/transitive_closure_dag,
   :rdf/type :py/Function})

(def transitive_reduction
  "Returns transitive reduction of a directed graph"
  {:db/ident :networkx.algorithms/transitive_reduction,
   :rdf/type :py/Function})

(def transitivity
  "Compute graph transitivity, the fraction of all possible triangles present in G."
  {:db/ident :networkx.algorithms/transitivity,
   :rdf/type :py/Function})

(def traversal
  "traversal"
  {:db/ident :networkx.algorithms/traversal,
   :rdf/type :py/Function})

(def tree
  "tree"
  {:db/ident :networkx.algorithms/tree,
   :rdf/type :py/Function})

(def tree_all_pairs_lowest_common_ancestor
  "Yield the lowest common ancestor for sets of pairs in a tree."
  {:db/ident :networkx.algorithms/tree_all_pairs_lowest_common_ancestor,
   :rdf/type :py/Function})

(def triad_type
  "Returns the sociological triad type for a triad."
  {:db/ident :networkx.algorithms/triad_type,
   :rdf/type :py/Function})

(def triadic_census
  "Determines the triadic census of a directed graph."
  {:db/ident :networkx.algorithms/triadic_census,
   :rdf/type :py/Function})

(def triads
  "Functions for analyzing triads of a graph."
  {:db/ident :networkx.algorithms/triads,
   :rdf/type :py/Function})

(def triads_by_type
  "Returns a list of all triads for each triad type in a directed graph. There are exactly 16 different types of triads possible. Suppose 1, 2, 3 are three nodes, they will be classified as a particular triad type if their connections are as follows:"
  {:db/ident :networkx.algorithms/triads_by_type,
   :rdf/type :py/Function})

(def triangles
  "Compute the number of triangles."
  {:db/ident :networkx.algorithms/triangles,
   :rdf/type :py/Function})

(def trophic
  "Trophic levels"
  {:db/ident :networkx.algorithms/trophic,
   :rdf/type :py/Function})

(def trophic_differences
  "Compute the trophic differences of the edges of a directed graph."
  {:db/ident :networkx.algorithms/trophic_differences,
   :rdf/type :py/Function})

(def trophic_incoherence_parameter
  "Compute the trophic incoherence parameter of a graph."
  {:db/ident :networkx.algorithms/trophic_incoherence_parameter,
   :rdf/type :py/Function})

(def trophic_levels
  "Compute the trophic levels of nodes."
  {:db/ident :networkx.algorithms/trophic_levels,
   :rdf/type :py/Function})

(def tutte_polynomial
  "Returns the Tutte polynomial of `G`"
  {:db/ident :networkx.algorithms/tutte_polynomial,
   :rdf/type :py/Function})

(def unary
  "Unary operations on graphs"
  {:db/ident :networkx.algorithms/unary,
   :rdf/type :py/Function})

(def union
  "Combine graphs G and H. The names of nodes must be unique."
  {:db/ident :networkx.algorithms/union,
   :rdf/type :py/Function})

(def union_all
  "Returns the union of all graphs."
  {:db/ident :networkx.algorithms/union_all,
   :rdf/type :py/Function})

(def unweighted
  "Shortest path algorithms for unweighted graphs."
  {:db/ident :networkx.algorithms/unweighted,
   :rdf/type :py/Function})

(def vf2pp_all_isomorphisms
  "Yields all the possible mappings between G1 and G2."
  {:db/ident :networkx.algorithms/vf2pp_all_isomorphisms,
   :rdf/type :py/Function})

(def vf2pp_is_isomorphic
  "Examines whether G1 and G2 are isomorphic."
  {:db/ident :networkx.algorithms/vf2pp_is_isomorphic,
   :rdf/type :py/Function})

(def vf2pp_isomorphism
  "Return an isomorphic mapping between `G1` and `G2` if it exists."
  {:db/ident :networkx.algorithms/vf2pp_isomorphism,
   :rdf/type :py/Function})

(def vitality
  "Vitality measures."
  {:db/ident :networkx.algorithms/vitality,
   :rdf/type :py/Function})

(def volume
  "Returns the volume of a set of nodes."
  {:db/ident :networkx.algorithms/volume,
   :rdf/type :py/Function})

(def voronoi
  "Functions for computing the Voronoi cells of a graph."
  {:db/ident :networkx.algorithms/voronoi,
   :rdf/type :py/Function})

(def voronoi_cells
  "Returns the Voronoi cells centered at `center_nodes` with respect to the shortest-path distance metric."
  {:db/ident :networkx.algorithms/voronoi_cells,
   :rdf/type :py/Function})

(def voterank
  "Select a list of influential nodes in a graph using VoteRank algorithm"
  {:db/ident :networkx.algorithms/voterank,
   :rdf/type :py/Function})

(def voterank_alg
  "Algorithm to select influential nodes in a graph using VoteRank."
  {:db/ident :networkx.algorithms/voterank_alg,
   :rdf/type :py/Function})

(def weakly_connected
  "Weakly connected components."
  {:db/ident :networkx.algorithms/weakly_connected,
   :rdf/type :py/Function})

(def weakly_connected_components
  "Generate weakly connected components of G."
  {:db/ident :networkx.algorithms/weakly_connected_components,
   :rdf/type :py/Function})

(def weighted
  "Shortest path algorithms for weighted graphs."
  {:db/ident :networkx.algorithms/weighted,
   :rdf/type :py/Function})

(def weisfeiler_lehman_graph_hash
  "Return Weisfeiler Lehman (WL) graph hash."
  {:db/ident :networkx.algorithms/weisfeiler_lehman_graph_hash,
   :rdf/type :py/Function})

(def weisfeiler_lehman_subgraph_hashes
  "Return a dictionary of subgraph hashes by node."
  {:db/ident :networkx.algorithms/weisfeiler_lehman_subgraph_hashes,
   :rdf/type :py/Function})

(def wiener
  "Functions related to the Wiener index of a graph."
  {:db/ident :networkx.algorithms/wiener,
   :rdf/type :py/Function})

(def wiener_index
  "Returns the Wiener index of the given graph."
  {:db/ident :networkx.algorithms/wiener_index,
   :rdf/type :py/Function})

(def within_inter_cluster
  "Compute the ratio of within- and inter-cluster common neighbors of all node pairs in ebunch."
  {:db/ident :networkx.algorithms/within_inter_cluster,
   :rdf/type :py/Function})
