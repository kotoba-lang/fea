(ns kotoba.fea.postprocess
  "Post-processing — cljc port of kami-cae's `postprocess` module
  (kami-engine, retired per ADR-2607010000). Result fields are keywords:
  `:displacement` `:von-mises-stress` `:principal-stress` `:strain`
  `:temperature` `:mode-shape` `:safety-factor`, matching kami-cae's
  `ResultField` enum. No network, no I/O."
  (:require [kotoba.fea.vec3 :as v3]))

(def result-fields
  #{:displacement :von-mises-stress :principal-stress :strain
    :temperature :mode-shape :safety-factor})

(defn field-range
  "Compute `{:min :max :avg}` from a collection of scalar `values`.
  `{:min 0.0 :max 0.0 :avg 0.0}` for an empty collection."
  [values]
  (if (empty? values)
    {:min 0.0 :max 0.0 :avg 0.0}
    {:min (reduce min values)
     :max (reduce max values)
     :avg (/ (reduce + values) (count values))}))

(defn probe-point
  "Interpolate a displacement result at an arbitrary point.

  kami-cae's own implementation only had access to `AnalysisResult`
  (displacements, no node positions) at this call site, so instead of the
  documented inverse-distance weighting it fell back to returning the
  average displacement across all nodes (`_point` unused) — ported
  verbatim, including that limitation. A production implementation needs
  the mesh's node positions alongside the result to do real IDW."
  [result _point]
  (if (empty? (:displacement result))
    v3/zero
    (let [ds (:displacement result)]
      (v3/scale (reduce v3/add v3/zero ds) (/ 1.0 (count ds))))))

(defn export-color-map-data
  "Build per-node/per-element scalar data suitable for rendering via a
  color-map pipeline (kami-cae targeted `kami-eng-render`; left as plain
  data here — rendering is a host-adapter concern, see README).

  For `:displacement` the scalar is the displacement magnitude at each
  node. For element-based fields (`:von-mises-stress` `:strain` etc.) the
  values are per element (caller may average to nodes)."
  [result field]
  (case field
    :displacement (mapv v3/length (:displacement result))
    (:von-mises-stress :principal-stress) (vec (:stress result))
    :strain (vec (:strain result))
    :safety-factor
    ;; Safety factor = yield / stress. Placeholder yield of 250 MPa (mild
    ;; steel), matching kami-cae.
    (let [yield-stress 250.0e6]
      (mapv (fn [s] (if (> s 0.0) (/ yield-stress s) ##Inf)) (:stress result)))
    (:temperature :mode-shape)
    ;; Not computed by the linear-static solver; return zeros.
    (vec (repeat (count (:displacement result)) 0.0))))
