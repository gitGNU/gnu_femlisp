(in-package :cl-user)

(defpackage "DDO-FEMLISP"
  (:use "COMMON-LISP"
        "FL.MACROS" "FL.UTILITIES" "FL.DEBUG" "FL.AMOP"
        "FL.MATLISP" "FL.MESH"
        "FL.DISCRETIZATION"
        "FL.ELLSYS"
        "FL.ITERATION" "FL.MULTIGRID" "FL.GEOMG"
        "FL.STRATEGY"
        "FL.PLOT"
        "FL.PARALLEL"
        "LFARM" "MPI"
        "RELATIONS" "NET.SCIPOLIS.GRAPHS" "DDO")
  (:import-from "TREES"
                "MAKE-BINARY-TREE" "INSERT" "DOTREE" "PPRINT-TREE")
  (:export
   ;; ddo-mesh
   "MESH-GRAPH" "DISTRIBUTE-MESH"
   ;; ddo-solve
   "ASO-MAKE-DISTRIBUTED" "<DISTRIBUTED-JACOBI>" "MAKE-CONSISTENT"
   "PARALLEL-NORM"
   ;; ddo-sparseas
   "SPARSE-VECTOR-TO-LIST" "COMPARE-PV-LISTS")
  (:documentation
   "This package contains Femlisp parallelization stuff."))

(defpackage "DDO-FEMLISP-TEST"
  (:use "COMMON-LISP"
        "FL.MACROS" "FL.UTILITIES" "FL.DEBUG" "FL.AMOP" "FL.PARALLEL"
        "FL.MATLISP" "FL.MESH"
        "FL.DISCRETIZATION"
        "FL.ELLSYS"
        "FL.ITERATION" "FL.MULTIGRID" "FL.GEOMG"
        "FL.STRATEGY" "FL.PLOT"
        "FL.DOMAINS"
        "LFARM" "MPI"
        "RELATIONS" "NET.SCIPOLIS.GRAPHS" "DDO" "DDO-FEMLISP"
        "FL.DOMAINS" "FL.STRATEGY")
  (:import-from "FL.APPLICATION" "*RESULT*" "STORING")
  (:export)
  (:documentation
   "This package uses what the others provide."))


