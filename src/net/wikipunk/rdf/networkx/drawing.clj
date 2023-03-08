(ns net.wikipunk.rdf.networkx.drawing
  {:rdf/type :owl/Ontology}
  (:require   
   [net.wikipunk.rdf.py]))

(def arf_layout
  "Arf layout for networkx"
  {:db/ident :networkx.drawing/arf_layout,
   :rdf/type :py/Function})

(def bipartite_layout
  "Position nodes in two straight lines."
  {:db/ident :networkx.drawing/bipartite_layout,
   :rdf/type :py/Function})

(def circular_layout
  "Position nodes on a circle."
  {:db/ident :networkx.drawing/circular_layout,
   :rdf/type :py/Function})

(def draw
  "Draw the graph G with Matplotlib."
  {:db/ident :networkx.drawing/draw,
   :rdf/type :py/Function})

(def draw_circular
  "Draw the graph `G` with a circular layout."
  {:db/ident :networkx.drawing/draw_circular,
   :rdf/type :py/Function})

(def draw_kamada_kawai
  "Draw the graph `G` with a Kamada-Kawai force-directed layout."
  {:db/ident :networkx.drawing/draw_kamada_kawai,
   :rdf/type :py/Function})

(def draw_networkx
  "Draw the graph G using Matplotlib."
  {:db/ident :networkx.drawing/draw_networkx,
   :rdf/type :py/Function})

(def draw_networkx_edge_labels
  "Draw edge labels."
  {:db/ident :networkx.drawing/draw_networkx_edge_labels,
   :rdf/type :py/Function})

(def draw_networkx_edges
  "Draw the edges of the graph G."
  {:db/ident :networkx.drawing/draw_networkx_edges,
   :rdf/type :py/Function})

(def draw_networkx_labels
  "Draw node labels on the graph G."
  {:db/ident :networkx.drawing/draw_networkx_labels,
   :rdf/type :py/Function})

(def draw_networkx_nodes
  "Draw the nodes of the graph G."
  {:db/ident :networkx.drawing/draw_networkx_nodes,
   :rdf/type :py/Function})

(def draw_planar
  "Draw a planar networkx graph `G` with planar layout."
  {:db/ident :networkx.drawing/draw_planar,
   :rdf/type :py/Function})

(def draw_random
  "Draw the graph `G` with a random layout."
  {:db/ident :networkx.drawing/draw_random,
   :rdf/type :py/Function})

(def draw_shell
  "Draw networkx graph `G` with shell layout."
  {:db/ident :networkx.drawing/draw_shell,
   :rdf/type :py/Function})

(def draw_spectral
  "Draw the graph `G` with a spectral 2D layout."
  {:db/ident :networkx.drawing/draw_spectral,
   :rdf/type :py/Function})

(def draw_spring
  "Draw the graph `G` with a spring layout."
  {:db/ident :networkx.drawing/draw_spring,
   :rdf/type :py/Function})

(def fruchterman_reingold_layout
  "Position nodes using Fruchterman-Reingold force-directed algorithm."
  {:db/ident :networkx.drawing/fruchterman_reingold_layout,
   :rdf/type :py/Function})

(def kamada_kawai_layout
  "Position nodes using Kamada-Kawai path-length cost-function."
  {:db/ident :networkx.drawing/kamada_kawai_layout,
   :rdf/type :py/Function})

(def layout
  "****** Layout ******"
  {:db/ident :networkx.drawing/layout,
   :rdf/type :py/Function})

(def multipartite_layout
  "Position nodes in layers of straight lines."
  {:db/ident :networkx.drawing/multipartite_layout,
   :rdf/type :py/Function})

(def nx_agraph
  "*************** Graphviz AGraph ***************"
  {:db/ident :networkx.drawing/nx_agraph,
   :rdf/type :py/Function})

(def nx_latex
  "***** LaTeX *****"
  {:db/ident :networkx.drawing/nx_latex,
   :rdf/type :py/Function})

(def nx_pydot
  "***** Pydot *****"
  {:db/ident :networkx.drawing/nx_pydot,
   :rdf/type :py/Function})

(def nx_pylab
  "********** Matplotlib **********"
  {:db/ident :networkx.drawing/nx_pylab,
   :rdf/type :py/Function})

(def planar_layout
  "Position nodes without edge intersections."
  {:db/ident :networkx.drawing/planar_layout,
   :rdf/type :py/Function})

(def random_layout
  "Position nodes uniformly at random in the unit square."
  {:db/ident :networkx.drawing/random_layout,
   :rdf/type :py/Function})

(def rescale_layout
  "Returns scaled position array to (-scale, scale) in all axes."
  {:db/ident :networkx.drawing/rescale_layout,
   :rdf/type :py/Function})

(def rescale_layout_dict
  "Return a dictionary of scaled positions keyed by node"
  {:db/ident :networkx.drawing/rescale_layout_dict,
   :rdf/type :py/Function})

(def shell_layout
  "Position nodes in concentric circles."
  {:db/ident :networkx.drawing/shell_layout,
   :rdf/type :py/Function})

(def spectral_layout
  "Position nodes using the eigenvectors of the graph Laplacian."
  {:db/ident :networkx.drawing/spectral_layout,
   :rdf/type :py/Function})

(def spiral_layout
  "Position nodes in a spiral layout."
  {:db/ident :networkx.drawing/spiral_layout,
   :rdf/type :py/Function})

(def spring_layout
  "Position nodes using Fruchterman-Reingold force-directed algorithm."
  {:db/ident :networkx.drawing/spring_layout,
   :rdf/type :py/Function})

(def to_latex
  "Return latex code to draw the graph(s) in `Gbunch`"
  {:db/ident :networkx.drawing/to_latex,
   :rdf/type :py/Function})

(def to_latex_raw
  "Return a string of the LaTeX/TikZ code to draw `G`"
  {:db/ident :networkx.drawing/to_latex_raw,
   :rdf/type :py/Function})

(def write_latex
  "Write the latex code to draw the graph(s) onto `path`."
  {:db/ident :networkx.drawing/write_latex,
   :rdf/type :py/Function})
