(in-package :cl-user)

(defpackage "MPI-WORKER"
  (:use "CL" "MPI"))

(defpackage "NET.SCIPOLIS.RELATIONS"
  (:nicknames "RELATIONS")
  (:use "COMMON-LISP"
        "FL.MACROS" "FL.UTILITIES" "FL.DEBUG" "FL.AMOP")
  (:import-from "TREES"
                "MAKE-BINARY-TREE" "RIGHT" "LEFT" "ROOT" "DATUM" "TEST" "PRED"
                "INSERT"
                "LOWER-BOUND-NODE-WITH-PATH" "LOWER-BOUND-NODE"
                "UPPER-BOUND-NODE-WITH-PATH" "UPPER-BOUND-NODE"
                "EXTREME-NODE-WITH-PATH"
                "MAKE-ITERATOR"
                "DOTREE" "PPRINT-TREE")
  (:export "LIST-COMPARISON" "MAKE-NUMBER-RELATION" "R-INSERT" "R-REMOVE" "R-SELECT" "R-SOME")
  (:documentation
   "This package provides relations built upon binary trees."))

(defpackage "DDO"
  (:use "COMMON-LISP"
        "FL.MACROS" "FL.UTILITIES" "FL.DEBUG" "FL.AMOP" "FL.PARALLEL"
        "NET.SCIPOLIS.RELATIONS"
        "CL-MPI-EXTENSIONS")
  (:import-from "TREES"
                "MAKE-BINARY-TREE" "RIGHT" "LEFT" "ROOT" "DATUM" "TEST" "PRED"
                "INSERT"
                "LOWER-BOUND-NODE-WITH-PATH" "LOWER-BOUND-NODE"
                "UPPER-BOUND-NODE-WITH-PATH" "UPPER-BOUND-NODE"
                "EXTREME-NODE-WITH-PATH"
                "MAKE-ITERATOR"
                "DOTREE" "PPRINT-TREE")
  (:import-from "CL-MPI"
                "MPI-INIT" "MPI-INITIALIZED" "MPI-COMM-SIZE" "MPI-COMM-RANK")
  (:export
   "*DEBUG-SHOW-DATA*"
   "*SYNCHRONIZATION-REAL-TIME*"
   "*COMMUNICATION-REAL-TIME*"
   "*COMMUNICATION-SIZE*"
   "NEW-LOCAL-ID" "LOCAL-ID" "DISTRIBUTED-P"
   "DISTRIBUTED-DATA" "RESET-DISTRIBUTED-OBJECTS"
   "DDO-MIXIN" "DDO-CONTAINER-MIXIN"
   "DISTRIBUTED-P" "DISTRIBUTED-CONTAINER-P"
   "DISTRIBUTED-SLOTS" "DISTRIBUTED-SLOT-NAMES" "DISTRIBUTED-SLOT-VALUES"
   "ALL-PROCESSORS" "OWNERS" "NEIGHBORS-FOR" "MASTERP"
   "ENSURE-DISTRIBUTED-CLASS" "MAKE-DISTRIBUTED-OBJECT" "MAKE-DISTRIBUTED-CONTAINER"
   "MINIMUM-ID-MERGER" "OP-MERGER"
   "INSERT-INTO-CHANGED"
   "SYNCHRONIZE" "*SYNCHRONIZATION-MERGER*"
   "DDO" "DDOX" "DDO-")
  (:DOCUMENTATION
   "This package provides distributed objects."))

(defpackage "DDO-TEST"
  (:use "COMMON-LISP"
        "FL.MACROS" "FL.UTILITIES" "FL.DEBUG" "FL.AMOP"
        "LFARM" "MPI"
        "RELATIONS" "DDO")
  (:export )
  (:documentation
   "This package uses what the others provide."))
