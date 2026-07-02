(ns kotoba.fea.solver
  "Linear-static FEA solver — cljc port of kami-cae's `solver` module
  (kami-engine, retired per ADR-2607010000).

  Currently supports meshes composed entirely of `:beam2` (1-D bar)
  elements, matching kami-cae's scope exactly — it never implemented
  assembly for the other element topologies (`AnalysisType` and
  `SolverMethod` beyond `DirectCholesky` were declared as enums but had
  no dispatch/implementation in kami-cae either; ported here as data only,
  see below). Each bar element has a unit cross-section area (A = 1 m^2)
  — callers scale Young's modulus accordingly.

  Assembly: for each bar element connecting nodes i and j, local stiffness
  k = A*E / L is assembled into the global matrix. Only translational DOFs
  (x, y, z per node -> 3*N total DOFs) are considered.

  Errors are reported via `ex-info` (idiomatic Clojure; kami-cae used a
  `thiserror` `SolverError` enum) with `:type` one of `:singular-matrix`
  `:unsupported-element` `:no-loads` `:node-set-not-found`. No network,
  no I/O."
  (:require [kotoba.fea.vec3 :as v3]))

;; ---------------------------------------------------------------------------
;; analysis / solver method selectors (data port — see namespace docstring)
;; ---------------------------------------------------------------------------

(def analysis-types
  #{:linear-static :nonlinear-static :modal :thermal-steady :thermal-transient :buckling})

(def default-solver-method {:type :direct-cholesky})

(defn conjugate-gradient-method [max-iter tolerance]
  {:type :conjugate-gradient :max-iter max-iter :tolerance tolerance})

(def gmres-method {:type :gmres})

;; ---------------------------------------------------------------------------
;; dense linear algebra helpers (educational, small problems)
;; ---------------------------------------------------------------------------

(defn- at [m n i j] (nth m (+ (* i n) j)))

(defn- cholesky-decompose
  "Decompose symmetric positive-definite `a` (row-major, n x n) as
  A = L*L^T. Throws `ex-info` `{:type :singular-matrix}` if not SPD."
  [a n]
  (reduce
   (fn [l i]
     (reduce
      (fn [l j]
        (let [sum (reduce + (map (fn [k] (* (at l n i k) (at l n j k))) (range j)))]
          (if (= i j)
            (let [diag (- (at a n i i) sum)]
              (when (<= diag 0.0)
                (throw (ex-info "singular stiffness matrix — check constraints"
                                 {:type :singular-matrix})))
              (assoc l (+ (* i n) j) (v3/sqrt* diag)))
            (let [denom (at l n j j)]
              (when (< (v3/abs* denom) 1e-30)
                (throw (ex-info "singular stiffness matrix — check constraints"
                                 {:type :singular-matrix})))
              (assoc l (+ (* i n) j) (/ (- (at a n i j) sum) denom))))))
      l
      (range (inc i))))
   (vec (repeat (* n n) 0.0))
   (range n)))

(defn- forward-substitute
  "Solve L*y = b."
  [l b n]
  (reduce
   (fn [y i]
     (let [sum (reduce + (map (fn [k] (* (at l n i k) (nth y k))) (range i)))]
       (conj y (/ (- (nth b i) sum) (at l n i i)))))
   []
   (range n)))

(defn- backward-substitute
  "Solve L^T*x = y."
  [l y n]
  (reduce
   (fn [x i]
     (let [sum (reduce + (map (fn [k] (* (at l n k i) (nth x k))) (range (inc i) n)))]
       (assoc x i (/ (- (nth y i) sum) (at l n i i)))))
   (vec (repeat n 0.0))
   (range (dec n) -1 -1)))

(defn cholesky-solve
  "Solve a symmetric positive-definite system Ax = b via Cholesky
  decomposition. `a` is stored row-major, dimension `n x n`."
  [a b n]
  (let [l (cholesky-decompose a n)
        y (forward-substitute l b n)]
    (backward-substitute l y n)))

;; ---------------------------------------------------------------------------
;; assembly
;; ---------------------------------------------------------------------------

(defn- add-at [v idx x] (update v idx + x))

(defn- assemble-stiffness
  "Assemble the global stiffness matrix (dense, row-major, `ndof x ndof`)
  from `mesh`'s `:beam2` elements."
  [mesh youngs-modulus ndof]
  (let [nodes (:nodes mesh)
        cross-section-area 1.0]
    (reduce
     (fn [k-global elem]
       (when (not= (:type elem) :beam2)
         (throw (ex-info "unsupported element type for this solver"
                          {:type :unsupported-element})))
       (let [[ni nj] (:nodes elem)
             pi (:position (nth nodes ni))
             pj (:position (nth nodes nj))
             delta (v3/sub pj pi)
             length (v3/length delta)]
         (if (< length 1e-15)
           k-global
           (let [dir (v3/scale delta (/ 1.0 length))
                 ke (/ (* cross-section-area youngs-modulus) length)
                 dofs-i [(* ni 3) (inc (* ni 3)) (+ (* ni 3) 2)]
                 dofs-j [(* nj 3) (inc (* nj 3)) (+ (* nj 3) 2)]]
             (reduce
              (fn [k-global a]
                (reduce
                 (fn [k-global b]
                   (let [val (* ke (nth dir a) (nth dir b))
                         ii (nth dofs-i a) ib (nth dofs-i b)
                         ji (nth dofs-j a) jb (nth dofs-j b)]
                     (-> k-global
                         (add-at (+ (* ii ndof) ib) val)
                         (add-at (+ (* ji ndof) jb) val)
                         (add-at (+ (* ii ndof) jb) (- val))
                         (add-at (+ (* ji ndof) ib) (- val)))))
                 k-global
                 (range 3)))
              k-global
              (range 3))))))
     (vec (repeat (* ndof ndof) 0.0))
     (:elements mesh))))

(defn- apply-force-bcs [f-global mesh bcs]
  (let [force-bcs (filter #(= (:type %) :force) bcs)]
    (when (empty? force-bcs)
      (throw (ex-info "no force boundary conditions specified" {:type :no-loads})))
    (reduce
     (fn [f {:keys [node-set value]}]
       (let [ids (or (get (:node-sets mesh) node-set)
                      (throw (ex-info (str "node set '" node-set "' not found in mesh")
                                       {:type :node-set-not-found :node-set node-set})))
             [vx vy vz] value]
         (reduce (fn [f nid]
                    (let [base (* nid 3)]
                      (-> f (add-at base vx) (add-at (inc base) vy) (add-at (+ base 2) vz))))
                  f ids)))
     f-global
     force-bcs)))

(defn- apply-displacement-bcs [k-global f-global mesh ndof bcs]
  (let [masks [:x :y :z]]
    (reduce
     (fn [[k f] bc]
       (if (not= (:type bc) :displacement)
         [k f]
         (let [{:keys [node-set dof-mask value]} bc
               ids (or (get (:node-sets mesh) node-set)
                        (throw (ex-info (str "node set '" node-set "' not found in mesh")
                                         {:type :node-set-not-found :node-set node-set})))
               vals value]
           (reduce
            (fn [[k f] nid]
              (let [base (* nid 3)]
                (reduce
                 (fn [[k f] c]
                   (if (not (contains? dof-mask (nth masks c)))
                     [k f]
                     (let [d (+ base c)
                           v (nth vals c)
                           f' (reduce (fn [f r]
                                        (if (= r d) f (update f r - (* (at k ndof r d) v))))
                                      f (range ndof))
                           k' (reduce (fn [k j] (-> k (assoc (+ (* d ndof) j) 0.0)
                                                       (assoc (+ (* j ndof) d) 0.0)))
                                      k (range ndof))
                           k'' (assoc k' (+ (* d ndof) d) 1.0)
                           f'' (assoc f' d v)]
                       [k'' f''])))
                 [k f]
                 (range 3))))
            [k f]
            ids))))
     [k-global f-global]
     bcs)))

(defn- stabilize
  "Stabilise unconstrained zero-stiffness DOFs (e.g. transverse DOFs of
  bar elements with no stiffness contribution). Without this the matrix
  is singular; a tiny diagonal value effectively fixes these DOFs at
  zero without affecting other results."
  [k-global ndof]
  (let [stab 1.0]
    (reduce
     (fn [k d]
       (if (< (v3/abs* (at k ndof d d)) 1e-30)
         (assoc k (+ (* d ndof) d) stab)
         k))
     k-global
     (range ndof))))

(defn solve-linear-static
  "Solve a linear-static FEA problem. `material` must have a
  `:linear-elastic` model. Returns
  `{:analysis-id :displacement :stress :strain :max-displacement
  :max-stress}` — `:displacement` is a vector of `[x y z]` per node,
  `:stress`/`:strain` are per-element (Von Mises stress approximated as
  |sigma| for bar elements, matching kami-cae)."
  [mesh material bcs]
  (when (not= (:type (:model material)) :linear-elastic)
    (throw (ex-info "unsupported element type for this solver" {:type :unsupported-element})))
  (let [youngs-modulus (:youngs-modulus (:model material))
        nodes (:nodes mesh)
        n-nodes (count nodes)
        ndof (* n-nodes 3)
        k0 (assemble-stiffness mesh youngs-modulus ndof)
        f0 (vec (repeat ndof 0.0))
        f1 (apply-force-bcs f0 mesh bcs)
        [k2 f2] (apply-displacement-bcs k0 f1 mesh ndof bcs)
        k3 (stabilize k2 ndof)
        u (cholesky-solve k3 f2 ndof)
        displacement (mapv (fn [i] [(nth u (* i 3)) (nth u (inc (* i 3))) (nth u (+ (* i 3) 2))])
                            (range n-nodes))
        max-displacement (reduce max 0.0 (map v3/length displacement))
        stress-strain
        (mapv (fn [elem]
                (let [[ni nj] (:nodes elem)
                      pi (:position (nth nodes ni))
                      pj (:position (nth nodes nj))
                      delta (v3/sub pj pi)
                      length (v3/length delta)
                      dir (v3/scale delta (/ 1.0 length))
                      ui (nth displacement ni)
                      uj (nth displacement nj)
                      eps (/ (v3/dot (v3/sub uj ui) dir) length)
                      sig (* youngs-modulus eps)]
                  {:strain eps :stress (v3/abs* sig)}))
              (:elements mesh))
        stress (mapv :stress stress-strain)
        strain (mapv :strain stress-strain)
        max-stress (reduce max 0.0 stress)]
    {:analysis-id "linear-static-0"
     :displacement displacement
     :stress stress
     :strain strain
     :max-displacement max-displacement
     :max-stress max-stress}))
