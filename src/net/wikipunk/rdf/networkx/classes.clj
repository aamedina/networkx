(ns net.wikipunk.rdf.networkx.classes
  {:rdf/type :owl/Ontology}
  (:require   
   [net.wikipunk.rdf.py]))

(def add_cycle
  "Add a cycle to the Graph G_to_add_to."
  {:db/ident :networkx.classes/add_cycle,
   :rdf/type :py/Function})

(def add_path
  "Add a path to the Graph G_to_add_to."
  {:db/ident :networkx.classes/add_path,
   :rdf/type :py/Function})

(def add_star
  "Add a star to Graph G_to_add_to."
  {:db/ident :networkx.classes/add_star,
   :rdf/type :py/Function})

(def all_neighbors
  "Returns all of the neighbors of a node in the graph."
  {:db/ident :networkx.classes/all_neighbors,
   :rdf/type :py/Function})

(def backends
  "Code to support various backends in a plugin dispatch architecture."
  {:db/ident :networkx.classes/backends,
   :rdf/type :py/Function})

(def common_neighbors
  "Returns the common neighbors of two nodes in a graph."
  {:db/ident :networkx.classes/common_neighbors,
   :rdf/type :py/Function})

(def coreviews
  "Views of core data structures such as nested Mappings (e.g. dict-of-dicts). These ``Views`` often restrict element access, with either the entire view or layers of nested mappings being read-only."
  {:db/ident :networkx.classes/coreviews,
   :rdf/type :py/Function})

(def create_empty_copy
  "Returns a copy of the graph G with all of the edges removed."
  {:db/ident :networkx.classes/create_empty_copy,
   :rdf/type :py/Function})

(def degree
  "Returns a degree view of single node or of nbunch of nodes. If nbunch is omitted, then return degrees of *all* nodes."
  {:db/ident :networkx.classes/degree,
   :rdf/type :py/Function})

(def degree_histogram
  "Returns a list of the frequency of each degree value."
  {:db/ident :networkx.classes/degree_histogram,
   :rdf/type :py/Function})

(def density
  "Returns the density of a graph."
  {:db/ident :networkx.classes/density,
   :rdf/type :py/Function})

(def digraph
  "Base class for directed graphs."
  {:db/ident :networkx.classes/digraph,
   :rdf/type :py/Function})

(def edge_subgraph
  "Returns a view of the subgraph induced by the specified edges."
  {:db/ident :networkx.classes/edge_subgraph,
   :rdf/type :py/Function})

(def edges
  "Returns an edge view of edges incident to nodes in nbunch."
  {:db/ident :networkx.classes/edges,
   :rdf/type :py/Function})

(def filters
  "Filter factories to hide or show sets of nodes and edges."
  {:db/ident :networkx.classes/filters,
   :rdf/type :py/Function})

(def freeze
  "Modify graph to prevent further change by adding or removing nodes or edges."
  {:db/ident :networkx.classes/freeze,
   :rdf/type :py/Function})

(def function
  "Functional interface to graph methods and assorted utilities."
  {:db/ident :networkx.classes/function,
   :rdf/type :py/Function})

(def get_edge_attributes
  "Get edge attributes from graph"
  {:db/ident :networkx.classes/get_edge_attributes,
   :rdf/type :py/Function})

(def get_node_attributes
  "Get node attributes from graph"
  {:db/ident :networkx.classes/get_node_attributes,
   :rdf/type :py/Function})

(def graph
  "Base class for undirected graphs."
  {:db/ident :networkx.classes/graph,
   :rdf/type :py/Function})

(def graphviews
  "View of Graphs as SubGraph, Reverse, Directed, Undirected."
  {:db/ident :networkx.classes/graphviews,
   :rdf/type :py/Function})

(def induced_subgraph
  "Returns a SubGraph view of `G` showing only nodes in nbunch."
  {:db/ident :networkx.classes/induced_subgraph,
   :rdf/type :py/Function})

(def is_directed
  "Return True if graph is directed."
  {:db/ident :networkx.classes/is_directed,
   :rdf/type :py/Function})

(def is_empty
  "Returns True if `G` has no edges."
  {:db/ident :networkx.classes/is_empty,
   :rdf/type :py/Function})

(def is_frozen
  "Returns True if graph is frozen."
  {:db/ident :networkx.classes/is_frozen,
   :rdf/type :py/Function})

(def is_negatively_weighted
  "Returns True if `G` has negatively weighted edges."
  {:db/ident :networkx.classes/is_negatively_weighted,
   :rdf/type :py/Function})

(def is_path
  "Returns whether or not the specified path exists."
  {:db/ident :networkx.classes/is_path,
   :rdf/type :py/Function})

(def is_weighted
  "Returns True if `G` has weighted edges."
  {:db/ident :networkx.classes/is_weighted,
   :rdf/type :py/Function})

(def multidigraph
  "Base class for MultiDiGraph."
  {:db/ident :networkx.classes/multidigraph,
   :rdf/type :py/Function})

(def multigraph
  "Base class for MultiGraph."
  {:db/ident :networkx.classes/multigraph,
   :rdf/type :py/Function})

(def neighbors
  "Returns a list of nodes connected to node n."
  {:db/ident :networkx.classes/neighbors,
   :rdf/type :py/Function})

(def nodes
  "Returns an iterator over the graph nodes."
  {:db/ident :networkx.classes/nodes,
   :rdf/type :py/Function})

(def nodes_with_selfloops
  "Returns an iterator over nodes with self loops."
  {:db/ident :networkx.classes/nodes_with_selfloops,
   :rdf/type :py/Function})

(def non_edges
  "Returns the non-existent edges in the graph."
  {:db/ident :networkx.classes/non_edges,
   :rdf/type :py/Function})

(def non_neighbors
  "Returns the non-neighbors of the node in the graph."
  {:db/ident :networkx.classes/non_neighbors,
   :rdf/type :py/Function})

(def number_of_edges
  "Returns the number of edges in the graph."
  {:db/ident :networkx.classes/number_of_edges,
   :rdf/type :py/Function})

(def number_of_nodes
  "Returns the number of nodes in the graph."
  {:db/ident :networkx.classes/number_of_nodes,
   :rdf/type :py/Function})

(def number_of_selfloops
  "Returns the number of selfloop edges."
  {:db/ident :networkx.classes/number_of_selfloops,
   :rdf/type :py/Function})

(def path_weight
  "Returns total cost associated with specified path and weight"
  {:db/ident :networkx.classes/path_weight,
   :rdf/type :py/Function})

(def reportviews
  "View Classes provide node, edge and degree \"views\" of a graph."
  {:db/ident :networkx.classes/reportviews,
   :rdf/type :py/Function})

(def restricted_view
  "Returns a view of `G` with hidden nodes and edges."
  {:db/ident :networkx.classes/restricted_view,
   :rdf/type :py/Function})

(def reverse_view
  "View of `G` with edge directions reversed"
  {:db/ident :networkx.classes/reverse_view,
   :rdf/type :py/Function})

(def selfloop_edges
  "Returns an iterator over selfloop edges."
  {:db/ident :networkx.classes/selfloop_edges,
   :rdf/type :py/Function})

(def set_edge_attributes
  "Sets edge attributes from a given value or dictionary of values."
  {:db/ident :networkx.classes/set_edge_attributes,
   :rdf/type :py/Function})

(def set_node_attributes
  "Sets node attributes from a given value or dictionary of values."
  {:db/ident :networkx.classes/set_node_attributes,
   :rdf/type :py/Function})

(def subgraph
  "Returns the subgraph induced on nodes in nbunch."
  {:db/ident :networkx.classes/subgraph,
   :rdf/type :py/Function})

(def subgraph_view
  "View of `G` applying a filter on nodes and edges."
  {:db/ident :networkx.classes/subgraph_view,
   :rdf/type :py/Function})

(def to_directed
  "Returns a directed view of the graph `graph`."
  {:db/ident :networkx.classes/to_directed,
   :rdf/type :py/Function})

(def to_undirected
  "Returns an undirected view of the graph `graph`."
  {:db/ident :networkx.classes/to_undirected,
   :rdf/type :py/Function})
