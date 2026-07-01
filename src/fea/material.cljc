(ns fea.material
  "Finite-element material definitions + a library with built-in
  engineering presets. Restored from kami-cae's `material` module
  (deleted PR #82).")

(defn linear-elastic
  [{:keys [youngs-modulus poissons-ratio density thermal-expansion
           thermal-conductivity specific-heat]}]
  {:kind :linear-elastic :youngs-modulus youngs-modulus :poissons-ratio poissons-ratio
   :density density :thermal-expansion thermal-expansion
   :thermal-conductivity thermal-conductivity :specific-heat specific-heat})

(defn hyperelastic [] {:kind :hyperelastic})
(defn elasto-plastic [] {:kind :elasto-plastic})

(defn fe-material [name model] {:name name :model model})

(defn library
  "A fresh, empty material library."
  []
  {:materials []})

(defn add
  "Add a material. Returns `[index library']`."
  [lib mat]
  [(count (:materials lib)) (update lib :materials conj mat)])

(defn get-material [lib name] (some #(when (= (:name %) name) %) (:materials lib)))

;; ---- built-in presets ------------------------------------------------

(defn steel-structural []
  (fe-material "Steel-Structural"
               (linear-elastic {:youngs-modulus 200.0e9 :poissons-ratio 0.3 :density 7850.0
                                 :thermal-expansion 12.0e-6 :thermal-conductivity 50.0 :specific-heat 490.0})))

(defn aluminum-6061 []
  (fe-material "Aluminum-6061"
               (linear-elastic {:youngs-modulus 68.9e9 :poissons-ratio 0.33 :density 2700.0
                                 :thermal-expansion 23.6e-6 :thermal-conductivity 167.0 :specific-heat 896.0})))

(defn titanium-6al4v []
  (fe-material "Titanium-6Al4V"
               (linear-elastic {:youngs-modulus 113.8e9 :poissons-ratio 0.342 :density 4430.0
                                 :thermal-expansion 8.6e-6 :thermal-conductivity 6.7 :specific-heat 526.0})))

(defn concrete []
  (fe-material "Concrete"
               (linear-elastic {:youngs-modulus 30.0e9 :poissons-ratio 0.2 :density 2400.0
                                 :thermal-expansion 10.0e-6 :thermal-conductivity 1.7 :specific-heat 880.0})))

(defn with-presets
  "A library pre-populated with all built-in presets."
  []
  {:materials [(steel-structural) (aluminum-6061) (titanium-6al4v) (concrete)]})
