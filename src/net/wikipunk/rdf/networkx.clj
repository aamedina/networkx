(ns net.wikipunk.rdf.networkx
  {:rdf/type :owl/Ontology}
  (:require
   [net.wikipunk.rdf.py]))

(def Graph
  "Base class for undirected graphs."
  {:db/ident        :networkx/Graph,
   :rdf/type        :owl/Class,
   :rdfs/subClassOf [:py/Object]})

(def DiGraph
  "Base class for directed graphs."
  {:db/ident        :networkx/DiGraph,
   :rdf/type        :owl/Class,
   :rdfs/subClassOf :networkx/Graph})

(def ArborescenceIterator
  "Iterate over all spanning arborescences of a graph in either increasing or decreasing cost."
  {:db/ident        :networkx/ArborescenceIterator,
   :rdf/type        :owl/Class,
   :rdfs/subClassOf :py/Object})

(def EdgePartition
  "An enum to store the state of an edge partition. The enum is written to the edges of a graph before being pasted to `kruskal_mst_edges`. Options are:"
  {:db/ident        :networkx/EdgePartition,
   :rdf/type        :owl/Class,
   :rdfs/subClassOf :py/Enum})

(def GraphMLReader
  "Read a GraphML document.  Produces NetworkX graph objects."
  {:db/ident        :networkx/GraphMLReader,
   :rdf/type        :owl/Class,
   :rdfs/subClassOf :networkx/GraphML})

(def GraphMLWriter
  "GraphMLWriter"
  {:db/ident        :networkx/GraphMLWriter,
   :rdf/type        :owl/Class,
   :rdfs/subClassOf :networkx/GraphML})

(def MultiGraph
  "An undirected graph class that can store multiedges."
  {:db/ident        :networkx/MultiGraph,
   :rdf/type        :owl/Class,
   :rdfs/subClassOf :networkx/Graph})

(def MultiDiGraph
  "A directed graph class that can store multiedges."
  {:db/ident        :networkx/MultiDiGraph,
   :rdf/type        :owl/Class,
   :rdfs/subClassOf [:networkx/MultiGraph
                     :networkx/DiGraph]})

(def PlanarEmbedding
  "Represents a planar graph with its planar embedding."
  {:db/ident        :networkx/PlanarEmbedding,
   :rdf/type        :owl/Class,
   :rdfs/subClassOf :networkx/DiGraph})

(def SpanningTreeIterator
  "Iterate over all spanning trees of a graph in either increasing or decreasing cost."
  {:db/ident        :networkx/SpanningTreeIterator,
   :rdf/type        :owl/Class,
   :rdfs/subClassOf :py/Object})

#_(def Collection
  "Collection"
  {:db/ident        :networkx/Collection,
   :rdf/type        :owl/Class,
   :rdfs/subClassOf [:collections/Sized
                     :collections/Iterable
                     :collections/Container]})

#_(def Generator
  "Generator"
  {:db/ident        :networkx/Generator,
   :rdf/type        :owl/Class,
   :rdfs/subClassOf [:collections/Iterator
                     :collections/Iterable]})

#_(def Iterator
  "Iterator"
  {:db/ident        :networkx/Iterator,
   :rdf/type        :owl/Class,
   :rdfs/subClassOf :collections/Iterable})

;; Exceptions

(def NetworkXException
  "Base class for exceptions in NetworkX."
  {:db/ident        :networkx/NetworkXException,
   :rdf/type        :owl/Class,
   :rdfs/subClassOf :py/Exception})

(def AmbiguousSolution
  "Raised if more than one valid solution exists for an intermediary step of an algorithm."
  {:db/ident        :networkx/AmbiguousSolution,
   :rdf/type        :owl/Class,
   :rdfs/subClassOf :networkx/NetworkXException})

(def ExceededMaxIterations
  "Raised if a loop iterates too many times without breaking."
  {:db/ident        :networkx/ExceededMaxIterations,
   :rdf/type        :owl/Class,
   :rdfs/subClassOf :networkx/NetworkXException})

(def HasACycle
  "Raised if a graph has a cycle when an algorithm expects that it will have no cycles."
  {:db/ident        :networkx/HasACycle,
   :rdf/type        :owl/Class,
   :rdfs/subClassOf :networkx/NetworkXException})

(def NetworkXAlgorithmError
  "Exception for unexpected termination of algorithms."
  {:db/ident        :networkx/NetworkXAlgorithmError,
   :rdf/type        :owl/Class,
   :rdfs/subClassOf :networkx/NetworkXException})

(def NetworkXError
  "Exception for a serious error in NetworkX"
  {:db/ident        :networkx/NetworkXError,
   :rdf/type        :owl/Class,
   :rdfs/subClassOf :networkx/NetworkXException})

(def NetworkXNoCycle
  "Exception for algorithms that should return a cycle when running on graphs where such a cycle does not exist."
  {:db/ident        :networkx/NetworkXNoCycle,
   :rdf/type        :owl/Class,
   :rdfs/subClassOf :networkx/NetworkXUnfeasible})

(def NetworkXNoPath
  "Exception for algorithms that should return a path when running on graphs where such a path does not exist."
  {:db/ident        :networkx/NetworkXNoPath,
   :rdf/type        :owl/Class,
   :rdfs/subClassOf :networkx/NetworkXUnfeasible})

(def NetworkXNotImplemented
  "Exception raised by algorithms not implemented for a type of graph."
  {:db/ident        :networkx/NetworkXNotImplemented,
   :rdf/type        :owl/Class,
   :rdfs/subClassOf :networkx/NetworkXException})

(def NetworkXPointlessConcept
  "Raised when a null graph is provided as input to an algorithm that cannot use it."
  {:db/ident        :networkx/NetworkXPointlessConcept,
   :rdf/type        :owl/Class,
   :rdfs/subClassOf :networkx/NetworkXException})

(def NetworkXTreewidthBoundExceeded
  "Exception raised when a treewidth bound has been provided and it has been exceeded"
  {:db/ident        :networkx/NetworkXTreewidthBoundExceeded,
   :rdf/type        :owl/Class,
   :rdfs/subClassOf :networkx/NetworkXException})

(def NetworkXUnbounded
  "Exception raised by algorithms trying to solve a maximization or a minimization problem instance that is unbounded."
  {:db/ident        :networkx/NetworkXUnbounded,
   :rdf/type        :owl/Class,
   :rdfs/subClassOf :networkx/NetworkXAlgorithmError})

(def NetworkXUnfeasible
  "Exception raised by algorithms trying to solve a problem instance that has no feasible solution."
  {:db/ident        :networkx/NetworkXUnfeasible,
   :rdf/type        :owl/Class,
   :rdfs/subClassOf :networkx/NetworkXAlgorithmError})

(def NodeNotFound
  "Exception raised if requested node is not present in the graph"
  {:db/ident        :networkx/NodeNotFound,
   :rdf/type        :owl/Class,
   :rdfs/subClassOf :networkx/NetworkXException})

(def NotATree
  "Raised when a function expects a tree (that is, a connected undirected graph with no cycles) but gets a non-tree graph as input instead."
  {:db/ident        :networkx/NotATree,
   :rdf/type        :owl/Class,
   :rdfs/subClassOf :networkx/NetworkXException})

(def PowerIterationFailedConvergence
  "Raised when the power iteration method fails to converge within a specified iteration limit."
  {:db/ident        :networkx/PowerIterationFailedConvergence,
   :rdf/type        :owl/Class,
   :rdfs/subClassOf :networkx/ExceededMaxIterations})
