(ns fea.postprocess
  "Post-processing: result field ranges, point probing, color-map data
  export. Restored from kami-cae's `postprocess` module (deleted PR #82).")

(def result-fields
  #{:displacement :von-mises-stress :principal-stress :strain :temperature
    :mode-shape :safety-factor})

(defn field-range-from-values
  "Min/max/avg of `values`, or all-zero if empty."
  [values]
  (if (empty? values)
    {:min 0.0 :max 0.0 :avg 0.0}
    {:min (reduce min values) :max (reduce max values) :avg (/ (reduce + 0.0 values) (count values))}))

(defn probe-point
  "Average displacement across all nodes (a safe default — matches the
  original: a real point probe needs mesh node positions alongside the
  result, which `AnalysisResult` doesn't carry; the original documents
  this as a simplification)."
  [result _point]
  (let [ds (:displacement result)]
    (if (empty? ds)
      [0.0 0.0 0.0]
      (let [n (count ds)
            sum (reduce (fn [[sx sy sz] [x y z]] [(+ sx x) (+ sy y) (+ sz z)]) [0.0 0.0 0.0] ds)]
        (mapv #(/ % n) sum)))))

(defn export-color-map-data
  "Per-node/per-element scalar data for `field`, suitable for rendering via
  kotoba-lang/engineer-render's color-map pipeline."
  [result field]
  (case field
    :displacement (mapv (fn [[x y z]] (Math/sqrt (+ (* x x) (* y y) (* z z)))) (:displacement result))
    (:von-mises-stress :principal-stress) (:stress result)
    :strain (:strain result)
    :safety-factor (let [yield-stress 250.0e6]
                      (mapv (fn [s] (if (pos? s) (/ yield-stress s) ##Inf)) (:stress result)))
    (:temperature :mode-shape) (vec (repeat (count (:displacement result)) 0.0))))
