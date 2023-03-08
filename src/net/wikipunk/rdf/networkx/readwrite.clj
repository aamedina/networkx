(ns net.wikipunk.rdf.networkx.readwrite
  {:rdf/type :owl/Ontology}
  (:require   
   [net.wikipunk.rdf.py]))

(def adjacency
  "adjacency"
  {:db/ident :networkx.readwrite/adjacency,
   :rdf/type :py/Function})

(def adjacency_data
  "Returns data in adjacency format that is suitable for JSON serialization and use in Javascript documents."
  {:db/ident :networkx.readwrite/adjacency_data,
   :rdf/type :py/Function})

(def adjacency_graph
  "Returns graph from adjacency data format."
  {:db/ident :networkx.readwrite/adjacency_graph,
   :rdf/type :py/Function})

(def adjlist
  "************** Adjacency List ************** Read and write NetworkX graphs as adjacency lists."
  {:db/ident :networkx.readwrite/adjlist,
   :rdf/type :py/Function})

(def cytoscape
  "cytoscape"
  {:db/ident :networkx.readwrite/cytoscape,
   :rdf/type :py/Function})

(def cytoscape_data
  "Returns data in Cytoscape JSON format (cyjs)."
  {:db/ident :networkx.readwrite/cytoscape_data,
   :rdf/type :py/Function})

(def cytoscape_graph
  "Create a NetworkX graph from a dictionary in cytoscape JSON format."
  {:db/ident :networkx.readwrite/cytoscape_graph,
   :rdf/type :py/Function})

(def edgelist
  "********** Edge Lists ********** Read and write NetworkX graphs as edge lists."
  {:db/ident :networkx.readwrite/edgelist,
   :rdf/type :py/Function})

(def forest_str
  "Creates a nice utf8 representation of a directed forest"
  {:db/ident :networkx.readwrite/forest_str,
   :rdf/type :py/Function})

(def from_graph6_bytes
  "Read a simple undirected graph in graph6 format from bytes."
  {:db/ident :networkx.readwrite/from_graph6_bytes,
   :rdf/type :py/Function})

(def from_sparse6_bytes
  "Read an undirected graph in sparse6 format from string."
  {:db/ident :networkx.readwrite/from_sparse6_bytes,
   :rdf/type :py/Function})

(def generate_adjlist
  "Generate a single line of the graph G in adjacency list format."
  {:db/ident :networkx.readwrite/generate_adjlist,
   :rdf/type :py/Function})

(def generate_edgelist
  "Generate a single line of the graph G in edge list format."
  {:db/ident :networkx.readwrite/generate_edgelist,
   :rdf/type :py/Function})

(def generate_gexf
  "Generate lines of GEXF format representation of G."
  {:db/ident :networkx.readwrite/generate_gexf,
   :rdf/type :py/Function})

(def generate_gml
  "Generate a single entry of the graph `G` in GML format."
  {:db/ident :networkx.readwrite/generate_gml,
   :rdf/type :py/Function})

(def generate_graphml
  "Generate GraphML lines for G"
  {:db/ident :networkx.readwrite/generate_graphml,
   :rdf/type :py/Function})

(def generate_multiline_adjlist
  "Generate a single line of the graph G in multiline adjacency list format."
  {:db/ident :networkx.readwrite/generate_multiline_adjlist,
   :rdf/type :py/Function})

(def generate_pajek
  "Generate lines in Pajek graph format."
  {:db/ident :networkx.readwrite/generate_pajek,
   :rdf/type :py/Function})

(def gexf
  "Read and write graphs in GEXF format."
  {:db/ident :networkx.readwrite/gexf,
   :rdf/type :py/Function})

(def gml
  "Read graphs in GML format."
  {:db/ident :networkx.readwrite/gml,
   :rdf/type :py/Function})

(def graph6
  "Functions for reading and writing graphs in the *graph6* format."
  {:db/ident :networkx.readwrite/graph6,
   :rdf/type :py/Function})

(def graphml
  "******* GraphML ******* Read and write graphs in GraphML format."
  {:db/ident :networkx.readwrite/graphml,
   :rdf/type :py/Function})

(def json_graph
  "********* JSON data ********* Generate and parse JSON serializable data for NetworkX graphs."
  {:db/ident :networkx.readwrite/json_graph,
   :rdf/type :py/Function})

(def leda
  "Read graphs in LEDA format."
  {:db/ident :networkx.readwrite/leda,
   :rdf/type :py/Function})

(def multiline_adjlist
  "************************* Multi-line Adjacency List ************************* Read and write NetworkX graphs as multi-line adjacency lists."
  {:db/ident :networkx.readwrite/multiline_adjlist,
   :rdf/type :py/Function})

(def node_link
  "node_link"
  {:db/ident :networkx.readwrite/node_link,
   :rdf/type :py/Function})

(def node_link_data
  "Returns data in node-link format that is suitable for JSON serialization and use in Javascript documents."
  {:db/ident :networkx.readwrite/node_link_data,
   :rdf/type :py/Function})

(def node_link_graph
  "Returns graph from node-link data format. Useful for de-serialization from JSON."
  {:db/ident :networkx.readwrite/node_link_graph,
   :rdf/type :py/Function})

(def pajek
  "***** Pajek ***** Read graphs in Pajek format."
  {:db/ident :networkx.readwrite/pajek,
   :rdf/type :py/Function})

(def parse_adjlist
  "Parse lines of a graph adjacency list representation."
  {:db/ident :networkx.readwrite/parse_adjlist,
   :rdf/type :py/Function})

(def parse_edgelist
  "Parse lines of an edge list representation of a graph."
  {:db/ident :networkx.readwrite/parse_edgelist,
   :rdf/type :py/Function})

(def parse_gml
  "Parse GML graph from a string or iterable."
  {:db/ident :networkx.readwrite/parse_gml,
   :rdf/type :py/Function})

(def parse_graphml
  "Read graph in GraphML format from string."
  {:db/ident :networkx.readwrite/parse_graphml,
   :rdf/type :py/Function})

(def parse_leda
  "Read graph in LEDA format from string or iterable."
  {:db/ident :networkx.readwrite/parse_leda,
   :rdf/type :py/Function})

(def parse_multiline_adjlist
  "Parse lines of a multiline adjacency list representation of a graph."
  {:db/ident :networkx.readwrite/parse_multiline_adjlist,
   :rdf/type :py/Function})

(def parse_pajek
  "Parse Pajek format graph from string or iterable."
  {:db/ident :networkx.readwrite/parse_pajek,
   :rdf/type :py/Function})

(def read_adjlist
  "Read graph in adjacency list format from path."
  {:db/ident :networkx.readwrite/read_adjlist,
   :rdf/type :py/Function})

(def read_edgelist
  "Read a graph from a list of edges."
  {:db/ident :networkx.readwrite/read_edgelist,
   :rdf/type :py/Function})

(def read_gexf
  "Read graph in GEXF format from path."
  {:db/ident :networkx.readwrite/read_gexf,
   :rdf/type :py/Function})

(def read_gml
  "Read graph in GML format from `path`."
  {:db/ident :networkx.readwrite/read_gml,
   :rdf/type :py/Function})

(def read_graph6
  "Read simple undirected graphs in graph6 format from path."
  {:db/ident :networkx.readwrite/read_graph6,
   :rdf/type :py/Function})

(def read_graphml
  "Read graph in GraphML format from path."
  {:db/ident :networkx.readwrite/read_graphml,
   :rdf/type :py/Function})

(def read_leda
  "Read graph in LEDA format from path."
  {:db/ident :networkx.readwrite/read_leda,
   :rdf/type :py/Function})

(def read_multiline_adjlist
  "Read graph in multi-line adjacency list format from path."
  {:db/ident :networkx.readwrite/read_multiline_adjlist,
   :rdf/type :py/Function})

(def read_pajek
  "Read graph in Pajek format from path."
  {:db/ident :networkx.readwrite/read_pajek,
   :rdf/type :py/Function})

(def read_sparse6
  "Read an undirected graph in sparse6 format from path."
  {:db/ident :networkx.readwrite/read_sparse6,
   :rdf/type :py/Function})

(def read_weighted_edgelist
  "Read a graph as list of edges with numeric weights."
  {:db/ident :networkx.readwrite/read_weighted_edgelist,
   :rdf/type :py/Function})

(def relabel_gexf_graph
  "Relabel graph using \"label\" node keyword for node label."
  {:db/ident :networkx.readwrite/relabel_gexf_graph,
   :rdf/type :py/Function})

(def sparse6
  "Functions for reading and writing graphs in the *sparse6* format."
  {:db/ident :networkx.readwrite/sparse6,
   :rdf/type :py/Function})

(def text
  "Text-based visual representations of graphs"
  {:db/ident :networkx.readwrite/text,
   :rdf/type :py/Function})

(def to_graph6_bytes
  "Convert a simple undirected graph to bytes in graph6 format."
  {:db/ident :networkx.readwrite/to_graph6_bytes,
   :rdf/type :py/Function})

(def to_sparse6_bytes
  "Convert an undirected graph to bytes in sparse6 format."
  {:db/ident :networkx.readwrite/to_sparse6_bytes,
   :rdf/type :py/Function})

(def tree
  "tree"
  {:db/ident :networkx.readwrite/tree,
   :rdf/type :py/Function})

(def tree_data
  "Returns data in tree format that is suitable for JSON serialization and use in Javascript documents."
  {:db/ident :networkx.readwrite/tree_data,
   :rdf/type :py/Function})

(def tree_graph
  "Returns graph from tree data format."
  {:db/ident :networkx.readwrite/tree_graph,
   :rdf/type :py/Function})

(def write_adjlist
  "Write graph G in single-line adjacency-list format to path."
  {:db/ident :networkx.readwrite/write_adjlist,
   :rdf/type :py/Function})

(def write_edgelist
  "Write graph as a list of edges."
  {:db/ident :networkx.readwrite/write_edgelist,
   :rdf/type :py/Function})

(def write_gexf
  "Write G in GEXF format to path."
  {:db/ident :networkx.readwrite/write_gexf,
   :rdf/type :py/Function})

(def write_gml
  "Write a graph `G` in GML format to the file or file handle `path`."
  {:db/ident :networkx.readwrite/write_gml,
   :rdf/type :py/Function})

(def write_graph6
  "Write a simple undirected graph to a path in graph6 format."
  {:db/ident :networkx.readwrite/write_graph6,
   :rdf/type :py/Function})

(def write_graphml
  "Write G in GraphML XML format to path"
  {:db/ident :networkx.readwrite/write_graphml,
   :rdf/type :py/Function})

(def write_graphml_lxml
  "Write G in GraphML XML format to path"
  {:db/ident :networkx.readwrite/write_graphml_lxml,
   :rdf/type :py/Function})

(def write_graphml_xml
  "Write G in GraphML XML format to path"
  {:db/ident :networkx.readwrite/write_graphml_xml,
   :rdf/type :py/Function})

(def write_multiline_adjlist
  "Write the graph G in multiline adjacency list format to path"
  {:db/ident :networkx.readwrite/write_multiline_adjlist,
   :rdf/type :py/Function})

(def write_pajek
  "Write graph in Pajek format to path."
  {:db/ident :networkx.readwrite/write_pajek,
   :rdf/type :py/Function})

(def write_sparse6
  "Write graph G to given path in sparse6 format."
  {:db/ident :networkx.readwrite/write_sparse6,
   :rdf/type :py/Function})

(def write_weighted_edgelist
  "Write graph G as a list of edges with numeric weights."
  {:db/ident :networkx.readwrite/write_weighted_edgelist,
   :rdf/type :py/Function})
