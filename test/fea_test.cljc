(ns fea-test
  "Restoration-fidelity tests — one per original kami-cae Rust test
  (kami-engine/kami-cae/src/lib.rs `mod tests`, deleted PR #82)."
  (:require [clojure.test :refer [deftest is testing]]
            [fea]
            [fea.mesh :as mesh]
            [fea.material :as material]
            [fea.boundary :as boundary]
            [fea.solver :as solver]
            [fea.postprocess :as postprocess]))

(deftest namespace-loads
  (testing "the restored CLJC namespace loads"
    (is (some? (the-ns 'fea)))))

;; mirrors `test_material_library_presets`
(deftest material-library-presets
  (let [lib (material/with-presets)]
    (is (= 4 (count (:materials lib))))
    (let [steel (material/get-material lib "Steel-Structural")]
      (is (< (Math/abs (- (:youngs-modulus (:model steel)) 200.0e9)) 1.0))
      (is (< (Math/abs (- (:poissons-ratio (:model steel)) 0.3)) 1e-6))
      (is (< (Math/abs (- (:density (:model steel)) 7850.0)) 1e-6)))
    (is (some? (material/get-material lib "Aluminum-6061")))
    (is (some? (material/get-material lib "Titanium-6Al4V")))
    (is (some? (material/get-material lib "Concrete")))))

;; mirrors `test_box_mesh_generation`
(deftest box-mesh-generation
  (let [m (mesh/generate-box-mesh 2.0 3.0 4.0 2)
        stats (mesh/mesh-stats m)]
    (is (= 27 (:node-count stats)))
    (is (= 8 (:element-count stats)))
    (is (> (:min-quality stats) 0.0))))

;; mirrors `test_boundary_condition_creation`
(deftest boundary-condition-creation
  (let [fix (boundary/displacement-bc "base" boundary/dof-all [0.0 0.0 0.0])]
    (is (boundary/dof-contains? (:dof-mask fix) boundary/dof-x))
    (is (boundary/dof-contains? (:dof-mask fix) boundary/dof-rz)))
  (let [conv (boundary/convection-bc "outer" 25.0 293.15)]
    (is (= :convection (:kind conv)))))

;; mirrors `test_1d_bar_fea_solve`
(deftest bar-fea-solve
  (let [m0 (mesh/fea-mesh)
        [n0 m1] (mesh/add-node m0 [0.0 0.0 0.0])
        [n1 m2] (mesh/add-node m1 [1.0 0.0 0.0])
        m3 (mesh/add-element m2 (mesh/beam2-element 0 [n0 n1]))
        m4 (-> m3 (mesh/create-node-set "fixed" [n0]) (mesh/create-node-set "load" [n1]))
        mat (material/steel-structural)
        bcs [(boundary/displacement-bc "fixed" boundary/dof-all [0.0 0.0 0.0])
             (boundary/force-bc "load" [1000.0 0.0 0.0])]
        [status result] (solver/solve-linear-static m4 mat bcs)]
    (is (= :ok status))
    (let [expected (/ 1000.0 200.0e9)
          computed (first (nth (:displacement result) 1))
          err (/ (Math/abs (- computed expected)) expected)]
      (is (< err 1e-6)))
    (is (< (Math/abs (- (first (:stress result)) 1000.0)) 1.0))))

;; mirrors `test_field_range_calculation`
(deftest field-range-calculation
  (let [values [1.0 5.0 3.0 7.0 2.0]
        range (postprocess/field-range-from-values values)]
    (is (< (Math/abs (- (:min range) 1.0)) 1e-12))
    (is (< (Math/abs (- (:max range) 7.0)) 1e-12))
    (is (< (Math/abs (- (:avg range) 3.6)) 1e-12))))
