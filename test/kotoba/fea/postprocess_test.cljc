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

;; ---------------------------------------------------------------------------
;; factor-of-safety — caller-supplied provenance-carrying allowable
;; ---------------------------------------------------------------------------

(def fos-opts
  {:allowable-stress 200.0e6
   :provenance {:source "test fixture, dated 2026-09-04" :basis :yield}})

(deftest factor-of-safety-basic-test
  (let [res (pp/factor-of-safety [100.0e6 50.0e6 0.0] fos-opts)]
    ;; allowable/stress per element; non-positive stress -> ##Inf (no load seen)
    (is (= [2.0 4.0 ##Inf] (:factors res)))
    ;; governing element is the smallest finite factor
    (is (= {:value 2.0 :element 0} (:min-factor res)))
    ;; the provenance is carried through verbatim
    (is (= (:provenance fos-opts) (:provenance res)))))

(deftest factor-of-safety-all-unloaded-test
  (is (nil? (:min-factor (pp/factor-of-safety [0.0 -1.0] fos-opts)))))

(deftest factor-of-safety-missing-allowable-throws-test
  (is (thrown-with-msg? #?(:clj clojure.lang.ExceptionInfo :cljs js/Error)
                        #"allowable-stress"
                        (pp/factor-of-safety [1.0e6] {}))))

(deftest factor-of-safety-negative-allowable-throws-test
  (is (thrown-with-msg? #?(:clj clojure.lang.ExceptionInfo :cljs js/Error)
                        #"positive"
                        (pp/factor-of-safety [1.0e6]
                                             {:allowable-stress -5.0e6
                                              :provenance {:source "x" :basis :yield}}))))

(deftest factor-of-safety-missing-provenance-throws-test
  ;; an allowable without source/basis is not an allowable — fail closed
  (is (thrown-with-msg? #?(:clj clojure.lang.ExceptionInfo :cljs js/Error)
                        #"provenance"
                        (pp/factor-of-safety [1.0e6] {:allowable-stress 200.0e6}))))

(deftest factor-of-safety-empty-source-throws-test
  (is (thrown-with-msg? #?(:clj clojure.lang.ExceptionInfo :cljs js/Error)
                        #"provenance"
                        (pp/factor-of-safety [1.0e6]
                                             {:allowable-stress 200.0e6
                                              :provenance {:source "" :basis :yield}}))))

(deftest factor-of-safety-non-sequential-throws-test
  (is (thrown-with-msg? #?(:clj clojure.lang.ExceptionInfo :cljs js/Error)
                        #"sequential"
                        (pp/factor-of-safety 1.0e6 fos-opts))))

(deftest factor-of-safety-accepts-lazy-seq-test
  ;; solver returns vectors, but lazy seqs from map must work too
  (let [res (pp/factor-of-safety (map #(* 1.0e6 %) [100 100]) fos-opts)]
    (is (= [2.0 2.0] (:factors res)))))
