(ns net.wikipunk.networkx.boot
  {:rdf/type :jsonld/Context})

(def networkx
  {:rdf/type    :rdfa/PrefixMapping
   :rdfa/uri    "https://wikipunk.net/networkx/"
   :rdfa/prefix "networkx"})

(def networkx.algorithms
  {:rdf/type    :rdfa/PrefixMapping
   :rdfa/uri    "https://wikipunk.net/networkx/algorithms/"
   :rdfa/prefix "networkx.algorithms"})

(def networkx.classes
  {:rdf/type    :rdfa/PrefixMapping
   :rdfa/uri    "https://wikipunk.net/networkx/classes/"
   :rdfa/prefix "networkx.classes"})

(def networkx.generators
  {:rdf/type    :rdfa/PrefixMapping
   :rdfa/uri    "https://wikipunk.net/networkx/generators/"
   :rdfa/prefix "networkx.generators"})

(def networkx.linalg
  {:rdf/type    :rdfa/PrefixMapping
   :rdfa/uri    "https://wikipunk.net/networkx/linalg/"
   :rdfa/prefix "networkx.linalg"})

(def networkx.drawing
  {:rdf/type    :rdfa/PrefixMapping
   :rdfa/uri    "https://wikipunk.net/networkx/drawing/"
   :rdfa/prefix "networkx.drawing"})

(def networkx.readwrite
  {:rdf/type    :rdfa/PrefixMapping
   :rdfa/uri    "https://wikipunk.net/networkx/readwrite/"
   :rdfa/prefix "networkx.readwrite"})
