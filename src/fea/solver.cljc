(ns fea.solver
  "FEA solvers: dense Cholesky linear-system solve + a linear-static bar
  (Beam2) assembly/BC-application/solve pipeline. Restored from kami-cae's
  `solver` module (deleted PR #82). Matrices are flat row-major vectors
  (`n*n` length) — coarse-grained per-analysis assembly/solve, not a
  per-frame hot loop, matching ADR-2607010930's CLJC scope."
  (:require [fea.boundary :as boundary]))

(def analysis-types #{:linear-static :nonlinear-static :modal :thermal-steady :thermal-transient :buckling})
(def solver-methods #{:direct-cholesky :conjugate-gradient :gmres})

;; ---- vector-3 helpers ----
(defn- v- [[ax ay az] [bx by bz]] [(- ax bx) (- ay by) (- az bz)])
(defn- v-dot [[ax ay az] [bx by bz]] (+ (* ax bx) (* ay by) (* az bz)))
(defn- v-length [v] (Math/sqrt (v-dot v v)))

;; ---- dense matrix helpers (flat row-major, length n*n) ----
(defn- mat-get [m n r c] (nth m (+ (* r n) c)))
(defn- mat-set [m n r c v] (assoc m (+ (* r n) c) v))
(defn- mat-add [m n r c v] (update m (+ (* r n) c) + v))

(defn cholesky-solve
  "Solve a symmetric positive-definite system Ax = b via Cholesky
  decomposition. `a` is flat row-major `n x n`. Returns `[:ok x]` or
  `[:error :singular-matrix]`."
  [a b n]
  (let [l0 (vec (repeat (* n n) 0.0))
        l (reduce
           (fn [l i]
             (reduce
              (fn [l j]
                (let [sum (reduce (fn [s k] (+ s (* (mat-get l n i k) (mat-get l n j k)))) 0.0 (range j))]
                  (if (= i j)
                    (let [diag (- (mat-get a n i i) sum)]
                      (if (<= diag 0.0)
                        (reduced ::singular)
                        (mat-set l n i j (Math/sqrt diag))))
                    (let [denom (mat-get l n j j)]
                      (if (< (Math/abs denom) 1e-30)
                        (reduced ::singular)
                        (mat-set l n i j (/ (- (mat-get a n i j) sum) denom)))))))
              l (range (inc i))))
           l0 (range n))]
    (if (= l ::singular)
      [:error :singular-matrix]
      (let [y (reduce
               (fn [y i]
                 (let [sum (reduce (fn [s k] (+ s (* (mat-get l n i k) (nth y k)))) 0.0 (range i))]
                   (assoc y i (/ (- (nth b i) sum) (mat-get l n i i)))))
               (vec (repeat n 0.0)) (range n))
            x (reduce
               (fn [x i]
                 (let [sum (reduce (fn [s k] (+ s (* (mat-get l n k i) (nth x k)))) 0.0 (range (inc i) n))]
                   (assoc x i (/ (- (nth y i) sum) (mat-get l n i i)))))
               (vec (repeat n 0.0)) (range (dec n) -1 -1))]
        [:ok x]))))

(defn- assemble-bar-stiffness
  "Assemble global stiffness for Beam2-only mesh. Returns `[:ok k-global]`
  or `[:error :unsupported-element]`."
  [mesh youngs-modulus ndof]
  (let [cross-section-area 1.0
        nodes (:nodes mesh)]
    (reduce
     (fn [k-global elem]
       (if (not= (first elem) :beam2)
         (reduced [:error :unsupported-element])
         (let [[n-i n-j] (nth elem 2)
               pi (:position (nth nodes n-i))
               pj (:position (nth nodes n-j))
               delta (v- pj pi)
               length (v-length delta)]
           (if (< length 1e-15)
             k-global
             (let [dir (mapv #(/ % length) delta)
                   ke (/ (* cross-section-area youngs-modulus) length)
                   dofs-i [(* n-i 3) (inc (* n-i 3)) (+ 2 (* n-i 3))]
                   dofs-j [(* n-j 3) (inc (* n-j 3)) (+ 2 (* n-j 3))]]
               (reduce
                (fn [k-global [a b]]
                  (let [val (* ke (nth dir a) (nth dir b))]
                    (-> k-global
                        (mat-add ndof (nth dofs-i a) (nth dofs-i b) val)
                        (mat-add ndof (nth dofs-j a) (nth dofs-j b) val)
                        (mat-add ndof (nth dofs-i a) (nth dofs-j b) (- val))
                        (mat-add ndof (nth dofs-j a) (nth dofs-i b) (- val)))))
                k-global (for [a (range 3) b (range 3)] [a b]))))))
       )
     (vec (repeat (* ndof ndof) 0.0)) (:elements mesh))))

(defn- apply-force-bcs [f-global bcs node-sets]
  (reduce
   (fn [f-global bc]
     (if (not= (:kind bc) :force)
       f-global
       (if-let [ids (get node-sets (:node-set bc))]
         (let [[vx vy vz] (:value bc)]
           (reduce (fn [f-global nid]
                     (let [base (* nid 3)]
                       (-> f-global (update base + vx) (update (inc base) + vy) (update (+ base 2) + vz))))
                   f-global ids))
         (reduced [:error (str "node set '" (:node-set bc) "' not found")]))))
   f-global bcs))

(defn- apply-displacement-bcs [k-global f-global bcs node-sets ndof]
  (reduce
   (fn [[k-global f-global] bc]
     (if (not= (:kind bc) :displacement)
       [k-global f-global]
       (if-let [ids (get node-sets (:node-set bc))]
         (let [vals (:value bc)
               masks [boundary/dof-x boundary/dof-y boundary/dof-z]]
           (reduce
            (fn [[k-global f-global] nid]
              (let [base (* nid 3)]
                (reduce
                 (fn [[k-global f-global] c]
                   (if-not (boundary/dof-contains? (:dof-mask bc) (nth masks c))
                     [k-global f-global]
                     (let [d (+ base c)
                           v (nth vals c)
                           f-global (reduce (fn [f-global r]
                                               (if (= r d) f-global
                                                   (update f-global r - (* (mat-get k-global ndof r d) v))))
                                             f-global (range ndof))
                           k-global (reduce (fn [k-global j]
                                               (-> k-global (mat-set ndof d j 0.0) (mat-set ndof j d 0.0)))
                                             k-global (range ndof))
                           k-global (mat-set k-global ndof d d 1.0)
                           f-global (assoc f-global d v)]
                       [k-global f-global])))
                 [k-global f-global] (range 3))))
            [k-global f-global] ids))
         (reduced [:error (str "node set '" (:node-set bc) "' not found")]))))
   [k-global f-global] bcs))

(defn solve-linear-static
  "Solve a linear-static FEA problem (Beam2/bar elements only, unit
  cross-section area). Returns `[:ok result]` or `[:error kind-or-msg]`."
  [mesh material bcs]
  (if (not= (:kind (:model material)) :linear-elastic)
    [:error :unsupported-element]
    (let [youngs-modulus (:youngs-modulus (:model material))
          n-nodes (count (:nodes mesh))
          ndof (* n-nodes 3)
          k-result (assemble-bar-stiffness mesh youngs-modulus ndof)]
      (if (and (vector? k-result) (= (first k-result) :error))
        k-result
        (let [k-global k-result
              f-global0 (vec (repeat ndof 0.0))
              f-result (apply-force-bcs f-global0 bcs (:node-sets mesh))]
          (if (and (vector? f-result) (keyword? (first f-result)) (= (first f-result) :error))
            f-result
            (if-not (some #(= (:kind %) :force) bcs)
              [:error :no-loads]
              (let [[k-global f-global] (apply-displacement-bcs k-global f-result bcs (:node-sets mesh) ndof)
                    stab 1.0
                    k-global (reduce (fn [k-global d]
                                        (if (< (Math/abs (mat-get k-global ndof d d)) 1e-30)
                                          (mat-set k-global ndof d d stab)
                                          k-global))
                                      k-global (range ndof))
                    [status u] (cholesky-solve k-global f-global ndof)]
                (if (= status :error)
                  [:error u]
                  (let [displacement (mapv (fn [i] [(nth u (* i 3)) (nth u (inc (* i 3))) (nth u (+ 2 (* i 3)))])
                                            (range n-nodes))
                        max-disp (reduce max 0.0 (map v-length displacement))
                        elem-results
                        (mapv (fn [elem]
                                (let [[n-i n-j] (nth elem 2)
                                      pi (:position (nth (:nodes mesh) n-i))
                                      pj (:position (nth (:nodes mesh) n-j))
                                      delta (v- pj pi)
                                      length (v-length delta)
                                      dir (mapv #(/ % length) delta)
                                      ui (nth displacement n-i)
                                      uj (nth displacement n-j)
                                      eps (/ (v-dot (v- uj ui) dir) length)
                                      sig (* youngs-modulus eps)]
                                  [eps (Math/abs sig)]))
                              (:elements mesh))
                        strain (mapv first elem-results)
                        stress (mapv second elem-results)
                        max-stress (reduce max 0.0 stress)]
                    [:ok {:analysis-id "linear-static-0" :displacement displacement
                          :stress stress :strain strain
                          :max-displacement max-disp :max-stress max-stress}]))))))))))
