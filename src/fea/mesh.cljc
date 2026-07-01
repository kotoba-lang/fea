(ns fea.mesh
  "FEA mesh generation: nodes, elements, node/element sets, mesh stats, box
  mesher. Restored from kami-cae's `mesh` module (kami-engine/kami-cae/
  src/lib.rs, deleted PR #82). Node/element IDs are plain integers (the
  original's `NodeId(u32)`/`ElementId(u32)` newtypes add no behavior).
  Elements are `[kind id node-ids]` vectors (e.g. `[:beam2 0 [0 1]]`) — the
  original's enum variants each carry a fixed-size node array.")

(def element-kinds #{:beam2 :tri3 :tri6 :quad4 :tet4 :tet10 :hex8})

(defn beam2-element [id node-ids] [:beam2 id node-ids])
(defn tri3-element [id node-ids] [:tri3 id node-ids])
(defn tri6-element [id node-ids] [:tri6 id node-ids])
(defn quad4-element [id node-ids] [:quad4 id node-ids])
(defn tet4-element [id node-ids] [:tet4 id node-ids])
(defn tet10-element [id node-ids] [:tet10 id node-ids])
(defn hex8-element [id node-ids] [:hex8 id node-ids])

(defn element-id [element] (nth element 1))
(defn element-kind [element] (nth element 0))
(defn element-node-ids [element] (nth element 2))

(def element-orders #{:linear :quadratic})

(defn mesh-config
  ([] (mesh-config {}))
  ([{:keys [element-size min-size max-size curvature-refinement quality-threshold order]
     :or {element-size 1.0 min-size 0.1 max-size 10.0
          curvature-refinement true quality-threshold 0.3 order :linear}}]
   {:element-size element-size :min-size min-size :max-size max-size
    :curvature-refinement curvature-refinement :quality-threshold quality-threshold :order order}))

(defn fea-mesh
  "A fresh, empty FEA mesh."
  []
  {:nodes [] :elements [] :node-sets {} :element-sets {}})

(defn add-node
  "Append a node at `position` (`[x y z]`). Returns `[id mesh']`."
  [mesh position]
  (let [id (count (:nodes mesh))]
    [id (update mesh :nodes conj {:id id :position position})]))

(defn add-element [mesh element] (update mesh :elements conj element))

(defn create-node-set [mesh name ids] (assoc-in mesh [:node-sets name] ids))

(defn mesh-stats
  "Basic quality statistics. Element quality is a placeholder of 1.0 for
  all well-formed elements (matches the original — a real mesher would
  compute aspect-ratio/Jacobian metrics)."
  [mesh]
  (let [n-elements (count (:elements mesh))
        quality 1.0]
    {:node-count (count (:nodes mesh))
     :element-count n-elements
     :min-quality (if (zero? n-elements) 0.0 quality)
     :avg-quality (if (zero? n-elements) 0.0 quality)}))

(defn generate-box-mesh
  "Generate a regular hexahedral mesh of a `width` x `height` x `depth` box
  with `divisions` cells along each axis."
  [width height depth divisions]
  (let [n (inc divisions)
        node-idx (fn [iz iy ix] (+ (* iz n n) (* iy n) ix))
        positions (for [iz (range n) iy (range n) ix (range n)]
                    [(* (/ (double ix) divisions) width)
                     (* (/ (double iy) divisions) height)
                     (* (/ (double iz) divisions) depth)])
        nodes (mapv (fn [id p] {:id id :position p}) (range) positions)
        elements (vec
                  (for [iz (range divisions) iy (range divisions) ix (range divisions)
                        :let [eid (+ (* iz divisions divisions) (* iy divisions) ix)
                              base (fn [dz dy dx] (node-idx (+ iz dz) (+ iy dy) (+ ix dx)))
                              node-ids [(base 0 0 0) (base 0 0 1) (base 0 1 1) (base 0 1 0)
                                        (base 1 0 0) (base 1 0 1) (base 1 1 1) (base 1 1 0)]]]
                    (hex8-element eid node-ids)))]
    {:nodes nodes :elements elements :node-sets {} :element-sets {}}))
