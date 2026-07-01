(ns fea.boundary
  "Boundary conditions: DOF masks + prescribed displacement/force/pressure/
  temperature/convection. Restored from kami-cae's `boundary` module
  (deleted PR #82). A DOF mask is a plain integer bitmask (matches the
  original's `DofMask(u8)` newtype).")

(def dof-x 2r000001)
(def dof-y 2r000010)
(def dof-z 2r000100)
(def dof-rx 2r001000)
(def dof-ry 2r010000)
(def dof-rz 2r100000)
(def dof-all 2r111111)

(defn dof-union [a b] (bit-or a b))
(defn dof-contains? [mask other] (= (bit-and mask other) other))

;; BoundaryCondition variants
(defn displacement-bc [node-set dof-mask value]
  {:kind :displacement :node-set node-set :dof-mask dof-mask :value value})
(defn force-bc [node-set value] {:kind :force :node-set node-set :value value})
(defn pressure-bc [face-set value] {:kind :pressure :face-set face-set :value value})
(defn temperature-bc [node-set value] {:kind :temperature :node-set node-set :value value})
(defn convection-bc [face-set coefficient ambient-temp]
  {:kind :convection :face-set face-set :coefficient coefficient :ambient-temp ambient-temp})
