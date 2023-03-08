(ns net.wikipunk.rdf.networkx.generators
  {:rdf/type :owl/Ontology}
  (:require   
   [net.wikipunk.rdf.py]))

(def LCF_graph
  "Return the cubic graph specified in LCF notation."
  {:db/ident :networkx.generators/LCF_graph,
   :rdf/type :py/Function})

(def LFR_benchmark_graph
  "Returns the LFR benchmark graph."
  {:db/ident :networkx.generators/LFR_benchmark_graph,
   :rdf/type :py/Function})

(def atlas
  "Generators for the small graph atlas."
  {:db/ident :networkx.generators/atlas,
   :rdf/type :py/Function})

(def balanced_tree
  "Returns the perfectly balanced `r`-ary tree of height `h`."
  {:db/ident :networkx.generators/balanced_tree,
   :rdf/type :py/Function})

(def barabasi_albert_graph
  "Returns a random graph using Barabási–Albert preferential attachment"
  {:db/ident :networkx.generators/barabasi_albert_graph,
   :rdf/type :py/Function})

(def barbell_graph
  "Returns the Barbell Graph: two complete graphs connected by a path."
  {:db/ident :networkx.generators/barbell_graph,
   :rdf/type :py/Function})

(def binomial_graph
  "Returns a $G_{n,p}$ random graph, also known as an Erdős-Rényi graph or a binomial graph."
  {:db/ident :networkx.generators/binomial_graph,
   :rdf/type :py/Function})

(def binomial_tree
  "Returns the Binomial Tree of order n."
  {:db/ident :networkx.generators/binomial_tree,
   :rdf/type :py/Function})

(def bull_graph
  "Returns the Bull Graph"
  {:db/ident :networkx.generators/bull_graph,
   :rdf/type :py/Function})

(def caveman_graph
  "Returns a caveman graph of `l` cliques of size `k`."
  {:db/ident :networkx.generators/caveman_graph,
   :rdf/type :py/Function})

(def chordal_cycle_graph
  "Returns the chordal cycle graph on `p` nodes."
  {:db/ident :networkx.generators/chordal_cycle_graph,
   :rdf/type :py/Function})

(def chvatal_graph
  "Returns the Chvátal Graph"
  {:db/ident :networkx.generators/chvatal_graph,
   :rdf/type :py/Function})

(def circulant_graph
  "Returns the circulant graph $Ci_n(x_1, x_2, ..., x_m)$ with $n$ nodes."
  {:db/ident :networkx.generators/circulant_graph,
   :rdf/type :py/Function})

(def circular_ladder_graph
  "Returns the circular ladder graph $CL_n$ of length n."
  {:db/ident :networkx.generators/circular_ladder_graph,
   :rdf/type :py/Function})

(def classic
  "Generators for some classic graphs."
  {:db/ident :networkx.generators/classic,
   :rdf/type :py/Function})

(def cographs
  "Generators for cographs"
  {:db/ident :networkx.generators/cographs,
   :rdf/type :py/Function})

(def community
  "Generators for classes of graphs used in studying social networks."
  {:db/ident :networkx.generators/community,
   :rdf/type :py/Function})

(def complete_graph
  "Return the complete graph `K_n` with n nodes."
  {:db/ident :networkx.generators/complete_graph,
   :rdf/type :py/Function})

(def complete_multipartite_graph
  "Returns the complete multipartite graph with the specified subset sizes."
  {:db/ident :networkx.generators/complete_multipartite_graph,
   :rdf/type :py/Function})

(def configuration_model
  "Returns a random graph with the given degree sequence."
  {:db/ident :networkx.generators/configuration_model,
   :rdf/type :py/Function})

(def connected_caveman_graph
  "Returns a connected caveman graph of `l` cliques of size `k`."
  {:db/ident :networkx.generators/connected_caveman_graph,
   :rdf/type :py/Function})

(def connected_watts_strogatz_graph
  "Returns a connected Watts–Strogatz small-world graph."
  {:db/ident :networkx.generators/connected_watts_strogatz_graph,
   :rdf/type :py/Function})

(def cubical_graph
  "Returns the 3-regular Platonic Cubical Graph"
  {:db/ident :networkx.generators/cubical_graph,
   :rdf/type :py/Function})

(def cycle_graph
  "Returns the cycle graph $C_n$ of cyclically connected nodes."
  {:db/ident :networkx.generators/cycle_graph,
   :rdf/type :py/Function})

(def davis_southern_women_graph
  "Returns Davis Southern women social network."
  {:db/ident :networkx.generators/davis_southern_women_graph,
   :rdf/type :py/Function})

(def degree_seq
  "Generate graphs with a given degree sequence or expected degree sequence."
  {:db/ident :networkx.generators/degree_seq,
   :rdf/type :py/Function})

(def degree_sequence_tree
  "Make a tree for the given degree sequence."
  {:db/ident :networkx.generators/degree_sequence_tree,
   :rdf/type :py/Function})

(def dense_gnm_random_graph
  "Returns a $G_{n,m}$ random graph."
  {:db/ident :networkx.generators/dense_gnm_random_graph,
   :rdf/type :py/Function})

(def desargues_graph
  "Returns the Desargues Graph"
  {:db/ident :networkx.generators/desargues_graph,
   :rdf/type :py/Function})

(def diamond_graph
  "Returns the Diamond graph"
  {:db/ident :networkx.generators/diamond_graph,
   :rdf/type :py/Function})

(def directed
  "Generators for some directed graphs, including growing network (GN) graphs and scale-free graphs."
  {:db/ident :networkx.generators/directed,
   :rdf/type :py/Function})

(def directed_configuration_model
  "Returns a directed_random graph with the given degree sequences."
  {:db/ident :networkx.generators/directed_configuration_model,
   :rdf/type :py/Function})

(def directed_havel_hakimi_graph
  "Returns a directed graph with the given degree sequences."
  {:db/ident :networkx.generators/directed_havel_hakimi_graph,
   :rdf/type :py/Function})

(def directed_joint_degree_graph
  "Generates a random simple directed graph with the joint degree."
  {:db/ident :networkx.generators/directed_joint_degree_graph,
   :rdf/type :py/Function})

(def dodecahedral_graph
  "Returns the Platonic Dodecahedral graph."
  {:db/ident :networkx.generators/dodecahedral_graph,
   :rdf/type :py/Function})

(def dorogovtsev_goltsev_mendes_graph
  "Returns the hierarchically constructed Dorogovtsev-Goltsev-Mendes graph."
  {:db/ident :networkx.generators/dorogovtsev_goltsev_mendes_graph,
   :rdf/type :py/Function})

(def dual_barabasi_albert_graph
  "Returns a random graph using dual Barabási–Albert preferential attachment"
  {:db/ident :networkx.generators/dual_barabasi_albert_graph,
   :rdf/type :py/Function})

(def duplication
  "Functions for generating graphs based on the \"duplication\" method."
  {:db/ident :networkx.generators/duplication,
   :rdf/type :py/Function})

(def duplication_divergence_graph
  "Returns an undirected graph using the duplication-divergence model."
  {:db/ident :networkx.generators/duplication_divergence_graph,
   :rdf/type :py/Function})

(def ego
  "Ego graph."
  {:db/ident :networkx.generators/ego,
   :rdf/type :py/Function})

(def ego_graph
  "Returns induced subgraph of neighbors centered at node n within a given radius."
  {:db/ident :networkx.generators/ego_graph,
   :rdf/type :py/Function})

(def empty_graph
  "Returns the empty graph with n nodes and zero edges."
  {:db/ident :networkx.generators/empty_graph,
   :rdf/type :py/Function})

(def erdos_renyi_graph
  "Returns a $G_{n,p}$ random graph, also known as an Erdős-Rényi graph or a binomial graph."
  {:db/ident :networkx.generators/erdos_renyi_graph,
   :rdf/type :py/Function})

(def expanders
  "Provides explicit constructions of expander graphs."
  {:db/ident :networkx.generators/expanders,
   :rdf/type :py/Function})

(def expected_degree_graph
  "Returns a random graph with given expected degrees."
  {:db/ident :networkx.generators/expected_degree_graph,
   :rdf/type :py/Function})

(def extended_barabasi_albert_graph
  "Returns an extended Barabási–Albert model graph."
  {:db/ident :networkx.generators/extended_barabasi_albert_graph,
   :rdf/type :py/Function})

(def fast_gnp_random_graph
  "Returns a $G_{n,p}$ random graph, also known as an Erdős-Rényi graph or a binomial graph."
  {:db/ident :networkx.generators/fast_gnp_random_graph,
   :rdf/type :py/Function})

(def florentine_families_graph
  "Returns Florentine families graph."
  {:db/ident :networkx.generators/florentine_families_graph,
   :rdf/type :py/Function})

(def frucht_graph
  "Returns the Frucht Graph."
  {:db/ident :networkx.generators/frucht_graph,
   :rdf/type :py/Function})

(def full_rary_tree
  "Creates a full r-ary tree of `n` nodes."
  {:db/ident :networkx.generators/full_rary_tree,
   :rdf/type :py/Function})

(def gaussian_random_partition_graph
  "Generate a Gaussian random partition graph."
  {:db/ident :networkx.generators/gaussian_random_partition_graph,
   :rdf/type :py/Function})

(def general_random_intersection_graph
  "Returns a random intersection graph with independent probabilities for connections between node and attribute sets."
  {:db/ident :networkx.generators/general_random_intersection_graph,
   :rdf/type :py/Function})

(def geographical_threshold_graph
  "Returns a geographical threshold graph."
  {:db/ident :networkx.generators/geographical_threshold_graph,
   :rdf/type :py/Function})

(def geometric
  "Generators for geometric graphs."
  {:db/ident :networkx.generators/geometric,
   :rdf/type :py/Function})

(def geometric_edges
  "Returns edge list of node pairs within `radius` of each other."
  {:db/ident :networkx.generators/geometric_edges,
   :rdf/type :py/Function})

(def gn_graph
  "Returns the growing network (GN) digraph with `n` nodes."
  {:db/ident :networkx.generators/gn_graph,
   :rdf/type :py/Function})

(def gnc_graph
  "Returns the growing network with copying (GNC) digraph with `n` nodes."
  {:db/ident :networkx.generators/gnc_graph,
   :rdf/type :py/Function})

(def gnm_random_graph
  "Returns a $G_{n,m}$ random graph."
  {:db/ident :networkx.generators/gnm_random_graph,
   :rdf/type :py/Function})

(def gnp_random_graph
  "Returns a $G_{n,p}$ random graph, also known as an Erdős-Rényi graph or a binomial graph."
  {:db/ident :networkx.generators/gnp_random_graph,
   :rdf/type :py/Function})

(def gnr_graph
  "Returns the growing network with redirection (GNR) digraph with `n` nodes and redirection probability `p`."
  {:db/ident :networkx.generators/gnr_graph,
   :rdf/type :py/Function})

(def graph_atlas
  "Returns graph number `i` from the Graph Atlas."
  {:db/ident :networkx.generators/graph_atlas,
   :rdf/type :py/Function})

(def graph_atlas_g
  "Returns the list of all graphs with up to seven nodes named in the Graph Atlas."
  {:db/ident :networkx.generators/graph_atlas_g,
   :rdf/type :py/Function})

(def grid_2d_graph
  "Returns the two-dimensional grid graph."
  {:db/ident :networkx.generators/grid_2d_graph,
   :rdf/type :py/Function})

(def grid_graph
  "Returns the *n*-dimensional grid graph."
  {:db/ident :networkx.generators/grid_graph,
   :rdf/type :py/Function})

(def havel_hakimi_graph
  "Returns a simple graph with given degree sequence constructed using the Havel-Hakimi algorithm."
  {:db/ident :networkx.generators/havel_hakimi_graph,
   :rdf/type :py/Function})

(def heawood_graph
  "Returns the Heawood Graph, a (3,6) cage."
  {:db/ident :networkx.generators/heawood_graph,
   :rdf/type :py/Function})

(def hexagonal_lattice_graph
  "Returns an `m` by `n` hexagonal lattice graph."
  {:db/ident :networkx.generators/hexagonal_lattice_graph,
   :rdf/type :py/Function})

(def hoffman_singleton_graph
  "Returns the Hoffman-Singleton Graph."
  {:db/ident :networkx.generators/hoffman_singleton_graph,
   :rdf/type :py/Function})

(def house_graph
  "Returns the House graph (square with triangle on top)"
  {:db/ident :networkx.generators/house_graph,
   :rdf/type :py/Function})

(def house_x_graph
  "Returns the House graph with a cross inside the house square."
  {:db/ident :networkx.generators/house_x_graph,
   :rdf/type :py/Function})

(def hypercube_graph
  "Returns the *n*-dimensional hypercube graph."
  {:db/ident :networkx.generators/hypercube_graph,
   :rdf/type :py/Function})

(def icosahedral_graph
  "Returns the Platonic Icosahedral graph."
  {:db/ident :networkx.generators/icosahedral_graph,
   :rdf/type :py/Function})

(def internet_as_graphs
  "Generates graphs resembling the Internet Autonomous System network"
  {:db/ident :networkx.generators/internet_as_graphs,
   :rdf/type :py/Function})

(def intersection
  "Generators for random intersection graphs."
  {:db/ident :networkx.generators/intersection,
   :rdf/type :py/Function})

(def interval_graph
  "Generates an interval graph for a list of intervals given."
  {:db/ident :networkx.generators/interval_graph,
   :rdf/type :py/Function})

(def inverse_line_graph
  "Returns the inverse line graph of graph G."
  {:db/ident :networkx.generators/inverse_line_graph,
   :rdf/type :py/Function})

(def is_valid_directed_joint_degree
  "Checks whether the given directed joint degree input is realizable"
  {:db/ident :networkx.generators/is_valid_directed_joint_degree,
   :rdf/type :py/Function})

(def is_valid_joint_degree
  "Checks whether the given joint degree dictionary is realizable."
  {:db/ident :networkx.generators/is_valid_joint_degree,
   :rdf/type :py/Function})

(def joint_degree_graph
  "Generates a random simple graph with the given joint degree dictionary."
  {:db/ident :networkx.generators/joint_degree_graph,
   :rdf/type :py/Function})

(def joint_degree_seq
  "Generate graphs with a given joint degree and directed joint degree"
  {:db/ident :networkx.generators/joint_degree_seq,
   :rdf/type :py/Function})

(def k_random_intersection_graph
  "Returns a intersection graph with randomly chosen attribute sets for each node that are of equal size (k)."
  {:db/ident :networkx.generators/k_random_intersection_graph,
   :rdf/type :py/Function})

(def karate_club_graph
  "Returns Zachary's Karate Club graph."
  {:db/ident :networkx.generators/karate_club_graph,
   :rdf/type :py/Function})

(def krackhardt_kite_graph
  "Returns the Krackhardt Kite Social Network."
  {:db/ident :networkx.generators/krackhardt_kite_graph,
   :rdf/type :py/Function})

(def ladder_graph
  "Returns the Ladder graph of length n."
  {:db/ident :networkx.generators/ladder_graph,
   :rdf/type :py/Function})

(def lattice
  "Functions for generating grid graphs and lattices"
  {:db/ident :networkx.generators/lattice,
   :rdf/type :py/Function})

(def les_miserables_graph
  "Returns coappearance network of characters in the novel Les Miserables."
  {:db/ident :networkx.generators/les_miserables_graph,
   :rdf/type :py/Function})

(def line
  "Functions for generating line graphs."
  {:db/ident :networkx.generators/line,
   :rdf/type :py/Function})

(def line_graph
  "Returns the line graph of the graph or digraph `G`."
  {:db/ident :networkx.generators/line_graph,
   :rdf/type :py/Function})

(def lollipop_graph
  "Returns the Lollipop Graph; `K_m` connected to `P_n`."
  {:db/ident :networkx.generators/lollipop_graph,
   :rdf/type :py/Function})

(def margulis_gabber_galil_graph
  "Returns the Margulis-Gabber-Galil undirected MultiGraph on `n^2` nodes."
  {:db/ident :networkx.generators/margulis_gabber_galil_graph,
   :rdf/type :py/Function})

(def moebius_kantor_graph
  "Returns the Moebius-Kantor graph."
  {:db/ident :networkx.generators/moebius_kantor_graph,
   :rdf/type :py/Function})

(def mycielski
  "Functions related to the Mycielski Operation and the Mycielskian family of graphs."
  {:db/ident :networkx.generators/mycielski,
   :rdf/type :py/Function})

(def mycielski_graph
  "Generator for the n_th Mycielski Graph."
  {:db/ident :networkx.generators/mycielski_graph,
   :rdf/type :py/Function})

(def mycielskian
  "Returns the Mycielskian of a simple, undirected graph G"
  {:db/ident :networkx.generators/mycielskian,
   :rdf/type :py/Function})

(def navigable_small_world_graph
  "Returns a navigable small-world graph."
  {:db/ident :networkx.generators/navigable_small_world_graph,
   :rdf/type :py/Function})

(def newman_watts_strogatz_graph
  "Returns a Newman–Watts–Strogatz small-world graph."
  {:db/ident :networkx.generators/newman_watts_strogatz_graph,
   :rdf/type :py/Function})

(def nonisomorphic_trees
  "Returns a list of nonisomporphic trees"
  {:db/ident :networkx.generators/nonisomorphic_trees,
   :rdf/type :py/Function})

(def null_graph
  "Returns the Null graph with no nodes or edges."
  {:db/ident :networkx.generators/null_graph,
   :rdf/type :py/Function})

(def number_of_nonisomorphic_trees
  "Returns the number of nonisomorphic trees"
  {:db/ident :networkx.generators/number_of_nonisomorphic_trees,
   :rdf/type :py/Function})

(def octahedral_graph
  "Returns the Platonic Octahedral graph."
  {:db/ident :networkx.generators/octahedral_graph,
   :rdf/type :py/Function})

(def paley_graph
  "Returns the Paley (p-1)/2-regular graph on p nodes."
  {:db/ident :networkx.generators/paley_graph,
   :rdf/type :py/Function})

(def pappus_graph
  "Returns the Pappus graph."
  {:db/ident :networkx.generators/pappus_graph,
   :rdf/type :py/Function})

(def partial_duplication_graph
  "Returns a random graph using the partial duplication model."
  {:db/ident :networkx.generators/partial_duplication_graph,
   :rdf/type :py/Function})

(def path_graph
  "Returns the Path graph `P_n` of linearly connected nodes."
  {:db/ident :networkx.generators/path_graph,
   :rdf/type :py/Function})

(def petersen_graph
  "Returns the Petersen graph."
  {:db/ident :networkx.generators/petersen_graph,
   :rdf/type :py/Function})

(def planted_partition_graph
  "Returns the planted l-partition graph."
  {:db/ident :networkx.generators/planted_partition_graph,
   :rdf/type :py/Function})

(def powerlaw_cluster_graph
  "Holme and Kim algorithm for growing graphs with powerlaw degree distribution and approximate average clustering."
  {:db/ident :networkx.generators/powerlaw_cluster_graph,
   :rdf/type :py/Function})

(def prefix_tree
  "Creates a directed prefix tree from a list of paths."
  {:db/ident :networkx.generators/prefix_tree,
   :rdf/type :py/Function})

(def prefix_tree_recursive
  "Recursively creates a directed prefix tree from a list of paths."
  {:db/ident :networkx.generators/prefix_tree_recursive,
   :rdf/type :py/Function})

(def random_clustered
  "Generate graphs with given degree and triangle sequence."
  {:db/ident :networkx.generators/random_clustered,
   :rdf/type :py/Function})

(def random_clustered_graph
  "Generate a random graph with the given joint independent edge degree and triangle degree sequence."
  {:db/ident :networkx.generators/random_clustered_graph,
   :rdf/type :py/Function})

(def random_cograph
  "Returns a random cograph with $2 ^ n$ nodes."
  {:db/ident :networkx.generators/random_cograph,
   :rdf/type :py/Function})

(def random_degree_sequence_graph
  "Returns a simple random graph with the given degree sequence."
  {:db/ident :networkx.generators/random_degree_sequence_graph,
   :rdf/type :py/Function})

(def random_geometric_graph
  "Returns a random geometric graph in the unit cube of dimensions `dim`."
  {:db/ident :networkx.generators/random_geometric_graph,
   :rdf/type :py/Function})

(def random_graphs
  "Generators for random graphs."
  {:db/ident :networkx.generators/random_graphs,
   :rdf/type :py/Function})

(def random_internet_as_graph
  "Generates a random undirected graph resembling the Internet AS network"
  {:db/ident :networkx.generators/random_internet_as_graph,
   :rdf/type :py/Function})

(def random_k_out_graph
  "Returns a random `k`-out graph with preferential attachment."
  {:db/ident :networkx.generators/random_k_out_graph,
   :rdf/type :py/Function})

(def random_kernel_graph
  "Returns an random graph based on the specified kernel."
  {:db/ident :networkx.generators/random_kernel_graph,
   :rdf/type :py/Function})

(def random_lobster
  "Returns a random lobster graph."
  {:db/ident :networkx.generators/random_lobster,
   :rdf/type :py/Function})

(def random_partition_graph
  "Returns the random partition graph with a partition of sizes."
  {:db/ident :networkx.generators/random_partition_graph,
   :rdf/type :py/Function})

(def random_powerlaw_tree
  "Returns a tree with a power law degree distribution."
  {:db/ident :networkx.generators/random_powerlaw_tree,
   :rdf/type :py/Function})

(def random_powerlaw_tree_sequence
  "Returns a degree sequence for a tree with a power law distribution."
  {:db/ident :networkx.generators/random_powerlaw_tree_sequence,
   :rdf/type :py/Function})

(def random_regular_graph
  "Returns a random $d$-regular graph on $n$ nodes."
  {:db/ident :networkx.generators/random_regular_graph,
   :rdf/type :py/Function})

(def random_shell_graph
  "Returns a random shell graph for the constructor given."
  {:db/ident :networkx.generators/random_shell_graph,
   :rdf/type :py/Function})

(def random_tree
  "Returns a uniformly random tree on `n` nodes."
  {:db/ident :networkx.generators/random_tree,
   :rdf/type :py/Function})

(def relaxed_caveman_graph
  "Returns a relaxed caveman graph."
  {:db/ident :networkx.generators/relaxed_caveman_graph,
   :rdf/type :py/Function})

(def ring_of_cliques
  "Defines a \"ring of cliques\" graph."
  {:db/ident :networkx.generators/ring_of_cliques,
   :rdf/type :py/Function})

(def scale_free_graph
  "Returns a scale-free directed graph."
  {:db/ident :networkx.generators/scale_free_graph,
   :rdf/type :py/Function})

(def sedgewick_maze_graph
  "Return a small maze with a cycle."
  {:db/ident :networkx.generators/sedgewick_maze_graph,
   :rdf/type :py/Function})

(def small
  "Various small and named graphs, together with some compact generators."
  {:db/ident :networkx.generators/small,
   :rdf/type :py/Function})

(def social
  "Famous social networks."
  {:db/ident :networkx.generators/social,
   :rdf/type :py/Function})

(def soft_random_geometric_graph
  "Returns a soft random geometric graph in the unit cube."
  {:db/ident :networkx.generators/soft_random_geometric_graph,
   :rdf/type :py/Function})

(def spectral_graph_forge
  "Returns a random simple graph with spectrum resembling that of `G`"
  {:db/ident :networkx.generators/spectral_graph_forge,
   :rdf/type :py/Function})

(def star_graph
  "Return the star graph"
  {:db/ident :networkx.generators/star_graph,
   :rdf/type :py/Function})

(def stochastic
  "Functions for generating stochastic graphs from a given weighted directed graph."
  {:db/ident :networkx.generators/stochastic,
   :rdf/type :py/Function})

(def stochastic_block_model
  "Returns a stochastic block model graph."
  {:db/ident :networkx.generators/stochastic_block_model,
   :rdf/type :py/Function})

(def stochastic_graph
  "Returns a right-stochastic representation of directed graph `G`."
  {:db/ident :networkx.generators/stochastic_graph,
   :rdf/type :py/Function})

(def sudoku
  "Generator for Sudoku graphs"
  {:db/ident :networkx.generators/sudoku,
   :rdf/type :py/Function})

(def sudoku_graph
  "Returns the n-Sudoku graph. The default value of n is 3."
  {:db/ident :networkx.generators/sudoku_graph,
   :rdf/type :py/Function})

(def tetrahedral_graph
  "Returns the 3-regular Platonic Tetrahedral graph."
  {:db/ident :networkx.generators/tetrahedral_graph,
   :rdf/type :py/Function})

(def thresholded_random_geometric_graph
  "Returns a thresholded random geometric graph in the unit cube."
  {:db/ident :networkx.generators/thresholded_random_geometric_graph,
   :rdf/type :py/Function})

(def trees
  "Functions for generating trees."
  {:db/ident :networkx.generators/trees,
   :rdf/type :py/Function})

(def triad_graph
  "Returns the triad graph with the given name."
  {:db/ident :networkx.generators/triad_graph,
   :rdf/type :py/Function})

(def triads
  "Functions that generate the triad graphs, that is, the possible digraphs on three nodes."
  {:db/ident :networkx.generators/triads,
   :rdf/type :py/Function})

(def triangular_lattice_graph
  "Returns the $m$ by $n$ triangular lattice graph."
  {:db/ident :networkx.generators/triangular_lattice_graph,
   :rdf/type :py/Function})

(def trivial_graph
  "Return the Trivial graph with one node (with label 0) and no edges."
  {:db/ident :networkx.generators/trivial_graph,
   :rdf/type :py/Function})

(def truncated_cube_graph
  "Returns the skeleton of the truncated cube."
  {:db/ident :networkx.generators/truncated_cube_graph,
   :rdf/type :py/Function})

(def truncated_tetrahedron_graph
  "Returns the skeleton of the truncated Platonic tetrahedron."
  {:db/ident :networkx.generators/truncated_tetrahedron_graph,
   :rdf/type :py/Function})

(def turan_graph
  "Return the Turan Graph"
  {:db/ident :networkx.generators/turan_graph,
   :rdf/type :py/Function})

(def tutte_graph
  "Returns the Tutte graph."
  {:db/ident :networkx.generators/tutte_graph,
   :rdf/type :py/Function})

(def uniform_random_intersection_graph
  "Returns a uniform random intersection graph."
  {:db/ident :networkx.generators/uniform_random_intersection_graph,
   :rdf/type :py/Function})

(def watts_strogatz_graph
  "Returns a Watts–Strogatz small-world graph."
  {:db/ident :networkx.generators/watts_strogatz_graph,
   :rdf/type :py/Function})

(def waxman_graph
  "Returns a Waxman random graph."
  {:db/ident :networkx.generators/waxman_graph,
   :rdf/type :py/Function})

(def wheel_graph
  "Return the wheel graph"
  {:db/ident :networkx.generators/wheel_graph,
   :rdf/type :py/Function})

(def windmill_graph
  "Generate a windmill graph. A windmill graph is a graph of `n` cliques each of size `k` that are all joined at one node. It can be thought of as taking a disjoint union of `n` cliques of size `k`, selecting one point from each, and contracting all of the selected points. Alternatively, one could generate `n` cliques of size `k-1` and one node that is connected to all other nodes in the graph."
  {:db/ident :networkx.generators/windmill_graph,
   :rdf/type :py/Function})
