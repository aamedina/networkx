(ns net.wikipunk.rdf.networkx.linalg
  {:rdf/type :owl/Ontology}
  (:require   
   [net.wikipunk.rdf.py]))

(def adjacency_matrix
  "Returns adjacency matrix of G."
  {:db/ident :networkx.linalg/adjacency_matrix,
   :rdf/type :py/Function})

(def adjacency_spectrum
  "Returns eigenvalues of the adjacency matrix of G."
  {:db/ident :networkx.linalg/adjacency_spectrum,
   :rdf/type :py/Function})

(def algebraic_connectivity
  "Returns the algebraic connectivity of an undirected graph."
  {:db/ident :networkx.linalg/algebraic_connectivity,
   :rdf/type :py/Function})

(def algebraicconnectivity
  "Algebraic connectivity and Fiedler vectors of undirected graphs."
  {:db/ident :networkx.linalg/algebraicconnectivity,
   :rdf/type :py/Function})

(def attr_matrix
  "Returns the attribute matrix using attributes from `G` as a numpy array."
  {:db/ident :networkx.linalg/attr_matrix,
   :rdf/type :py/Function})

(def attr_sparse_matrix
  "Returns a SciPy sparse array using attributes from G."
  {:db/ident :networkx.linalg/attr_sparse_matrix,
   :rdf/type :py/Function})

(def attrmatrix
  "Functions for constructing matrix-like objects from graph attributes."
  {:db/ident :networkx.linalg/attrmatrix,
   :rdf/type :py/Function})

(def bethe_hessian_matrix
  "Returns the Bethe Hessian matrix of G."
  {:db/ident :networkx.linalg/bethe_hessian_matrix,
   :rdf/type :py/Function})

(def bethe_hessian_spectrum
  "Returns eigenvalues of the Bethe Hessian matrix of G."
  {:db/ident :networkx.linalg/bethe_hessian_spectrum,
   :rdf/type :py/Function})

(def bethehessianmatrix
  "Bethe Hessian or deformed Laplacian matrix of graphs."
  {:db/ident :networkx.linalg/bethehessianmatrix,
   :rdf/type :py/Function})

(def directed_combinatorial_laplacian_matrix
  "Return the directed combinatorial Laplacian matrix of G."
  {:db/ident :networkx.linalg/directed_combinatorial_laplacian_matrix,
   :rdf/type :py/Function})

(def directed_laplacian_matrix
  "Returns the directed Laplacian matrix of G."
  {:db/ident :networkx.linalg/directed_laplacian_matrix,
   :rdf/type :py/Function})

(def directed_modularity_matrix
  "Returns the directed modularity matrix of G."
  {:db/ident :networkx.linalg/directed_modularity_matrix,
   :rdf/type :py/Function})

(def fiedler_vector
  "Returns the Fiedler vector of a connected undirected graph."
  {:db/ident :networkx.linalg/fiedler_vector,
   :rdf/type :py/Function})

(def graphmatrix
  "Adjacency matrix and incidence matrix of graphs."
  {:db/ident :networkx.linalg/graphmatrix,
   :rdf/type :py/Function})

(def incidence_matrix
  "Returns incidence matrix of G."
  {:db/ident :networkx.linalg/incidence_matrix,
   :rdf/type :py/Function})

(def laplacian_matrix
  "Returns the Laplacian matrix of G."
  {:db/ident :networkx.linalg/laplacian_matrix,
   :rdf/type :py/Function})

(def laplacian_spectrum
  "Returns eigenvalues of the Laplacian of G"
  {:db/ident :networkx.linalg/laplacian_spectrum,
   :rdf/type :py/Function})

(def laplacianmatrix
  "Laplacian matrix of graphs."
  {:db/ident :networkx.linalg/laplacianmatrix,
   :rdf/type :py/Function})

(def modularity_matrix
  "Returns the modularity matrix of G."
  {:db/ident :networkx.linalg/modularity_matrix,
   :rdf/type :py/Function})

(def modularity_spectrum
  "Returns eigenvalues of the modularity matrix of G."
  {:db/ident :networkx.linalg/modularity_spectrum,
   :rdf/type :py/Function})

(def modularitymatrix
  "Modularity matrix of graphs."
  {:db/ident :networkx.linalg/modularitymatrix,
   :rdf/type :py/Function})

(def normalized_laplacian_matrix
  "Returns the normalized Laplacian matrix of G."
  {:db/ident :networkx.linalg/normalized_laplacian_matrix,
   :rdf/type :py/Function})

(def normalized_laplacian_spectrum
  "Return eigenvalues of the normalized Laplacian of G"
  {:db/ident :networkx.linalg/normalized_laplacian_spectrum,
   :rdf/type :py/Function})

(def spectral_ordering
  "Compute the spectral_ordering of a graph."
  {:db/ident :networkx.linalg/spectral_ordering,
   :rdf/type :py/Function})

(def spectrum
  "Eigenvalue spectrum of graphs."
  {:db/ident :networkx.linalg/spectrum,
   :rdf/type :py/Function})

(def total_spanning_tree_weight
  "Returns the total weight of all spanning trees of `G`."
  {:db/ident :networkx.linalg/total_spanning_tree_weight,
   :rdf/type :py/Function})
