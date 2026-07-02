(ns kotoba.fea.material-loader
  "JVM classpath resource loading for `kotoba.fea.material`'s built-in
  presets (`resources/kami/fea/materials.edn`).

  Deliberately a plain `.clj` file, not `.cljc`: this is the one place in
  this repo that does I/O (a classpath resource read), and it is
  inherently JVM-only, so it's kept out of `kotoba.fea.material` (which
  stays 100% pure/portable) rather than reader-conditional-guarded inside
  it. Callers on other hosts (cljs, SCI, GraalVM) read
  `kotoba.fea.material/presets-resource` via their own host's I/O and pass
  the parsed EDN to `kotoba.fea.material/find-material`/`add-material`,
  which are pure. Mirrors `kami-scene-contracts/src/kami/scene/contracts.cljc`'s
  `load-edn-resource` EDN-authority pattern."
  (:require [clojure.edn :as edn]
            [clojure.java.io :as io]
            [kotoba.fea.material :as material]))

(defn load-edn-resource [path]
  (let [resource (io/resource path)]
    (when-not resource
      (throw (ex-info "missing fea resource" {:path path})))
    (edn/read-string (slurp resource))))

(defn presets
  "Built-in material presets, loaded from the EDN resource."
  []
  (load-edn-resource material/presets-resource))
