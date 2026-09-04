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

(def provenance-required-keys
  "Keys that must be present (non-nil) in a provenance map handed to
  [[factor-of-safety]]: where the allowable came from and which failure
  mode it bounds. No default is ever substituted — an unprovenanced
  allowable is not an allowable (same rule as `kotoba.fea.convergence`,
  which requires a caller-owned `fs`)."
  [:source :basis])

(defn- require-allowable!
  [{:keys [allowable-stress provenance] :as _opts}]
  (when-not (number? allowable-stress)
    (throw (ex-info "factor-of-safety: :allowable-stress [Pa] is required and must be a number"
                    {:type :missing-allowable})))
  (when-not (pos? allowable-stress)
    (throw (ex-info "factor-of-safety: :allowable-stress must be positive"
                    {:type :invalid-allowable :allowable-stress allowable-stress})))
  (when-not (and (map? provenance)
                 (every? #(contains? provenance %) provenance-required-keys)
                 (string? (:source provenance))
                 (seq (:source provenance)))
    (throw (ex-info (str "factor-of-safety: :provenance must carry "
                         provenance-required-keys
                         " — an unprovenanced allowable is not an allowable")
                    {:type :missing-provenance
                     :required provenance-required-keys
                     :provenance provenance}))))

(defn factor-of-safety
  "Per-element factor of safety against a **caller-supplied, provenance-carrying**
  allowable stress. This is the executable counterpart of
  [[export-color-map-data]]'s `:safety-factor` path, which keeps a
  hardcoded 250 MPa mild-steel placeholder for kami-cae parity — that
  placeholder must never be used to grade a real design.

      stresses (per-element, e.g. solver's von Mises values)
      opts    {:allowable-stress 570.0e6        ; Pa — e.g. AZ31B yield, from a
                                                ;      dated source the caller owns
               :provenance {:source \"AMS 4375, 2026-09 retrieval\"
                            :basis :yield}}

  Returns {:factors [...]                       ; allowable / stress per element,
                                                ; ##Inf where stress <= 0 (no load
                                                ; seen — matches export-color-map-data)
           :min-factor {:value v :element idx}  ; governing element (nil if all Inf)
           :allowable-stress a :provenance p}

  Fails closed with `ex-info` `:type` `:missing-allowable`
  `:invalid-allowable` or `:missing-provenance`. No physical constant is
  invented here; the allowable and its provenance are the caller's."
  [stresses {:keys [allowable-stress provenance] :as opts}]
  (require-allowable! opts)
  (when-not (sequential? stresses)
    (throw (ex-info "factor-of-safety: stresses must be a sequential collection"
                    {:type :invalid-stresses})))
  (let [factors (mapv (fn [s] (if (pos? s) (/ allowable-stress s) ##Inf)) stresses)
        governing (first (sort-by :value (keep-indexed
                                          (fn [i f] (when (< f ##Inf) {:value f :element i}))
                                          factors)))]
    {:factors factors
     :min-factor governing
     :allowable-stress allowable-stress
     :provenance provenance}))

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
