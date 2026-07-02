(ns kotoba.fea.postprocess-test
  "Parity port of kami-cae's `tests::test_field_range_calculation`, plus
  coverage of `export-color-map-data` and `probe-point`."
  (:require [clojure.test :refer [deftest is]]
            [kotoba.fea.postprocess :as pp]))

(deftest field-range-calculation-test
  (let [values [1.0 5.0 3.0 7.0 2.0]
        range (pp/field-range values)]
    (is (< (Math/abs (- (:min range) 1.0)) 1e-12))
    (is (< (Math/abs (- (:max range) 7.0)) 1e-12))
    (is (< (Math/abs (- (:avg range) 3.6)) 1e-12))))

(deftest field-range-empty-test
  (is (= {:min 0.0 :max 0.0 :avg 0.0} (pp/field-range []))))

(def sample-result
  {:analysis-id "linear-static-0"
   :displacement [[0.0 0.0 0.0] [3.0 4.0 0.0]]
   :stress [1000.0 2000.0]
   :strain [1e-6 2e-6]
   :max-displacement 5.0
   :max-stress 2000.0})

(deftest export-color-map-displacement-test
  (is (= [0.0 5.0] (pp/export-color-map-data sample-result :displacement))))

(deftest export-color-map-stress-test
  (is (= [1000.0 2000.0] (pp/export-color-map-data sample-result :von-mises-stress)))
  (is (= [1000.0 2000.0] (pp/export-color-map-data sample-result :principal-stress))))

(deftest export-color-map-strain-test
  (is (= [1e-6 2e-6] (pp/export-color-map-data sample-result :strain))))

(deftest export-color-map-safety-factor-test
  (is (= [250000.0 125000.0] (pp/export-color-map-data sample-result :safety-factor))))

(deftest export-color-map-zeros-test
  (is (= [0.0 0.0] (pp/export-color-map-data sample-result :temperature)))
  (is (= [0.0 0.0] (pp/export-color-map-data sample-result :mode-shape))))

(deftest probe-point-averages-displacement-test
  (is (= [1.5 2.0 0.0] (pp/probe-point sample-result [0.0 0.0 0.0]))))

(deftest probe-point-empty-test
  (is (= [0.0 0.0 0.0] (pp/probe-point {:displacement []} [0.0 0.0 0.0]))))
