;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; strategy-defp.lisp
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;
;;; Copyright (C) 2003 Nicolas Neuss, University of Heidelberg.
;;; All rights reserved.
;;; 
;;; Redistribution and use in source and binary forms, with or without
;;; modification, are permitted provided that the following conditions are
;;; met:
;;; 
;;; 1. Redistributions of source code must retain the above copyright
;;; notice, this list of conditions and the following disclaimer.
;;; 
;;; 2. Redistributions in binary form must reproduce the above copyright
;;; notice, this list of conditions and the following disclaimer in the
;;; documentation and/or other materials provided with the distribution.
;;; 
;;; THIS SOFTWARE IS PROVIDED ``AS IS'' AND ANY EXPRESS OR IMPLIED
;;; WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
;;; MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN
;;; NO EVENT SHALL THE AUTHOR, THE UNIVERSITY OF HEIDELBERG OR OTHER
;;; CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
;;; EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
;;; PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
;;; PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
;;; LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
;;; NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
;;; SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(in-package "COMMON-LISP-USER")

(defpackage "STRATEGY"
  (:use "COMMON-LISP" "MATLISP" "MACROS" "UTILITIES" "ALGEBRA" "MESH"
	"PROBLEM" "DISCRETIZATION" "CDR" "CDR-FE" "ITERATIONS" "TESTS" "GEOMG"
	"PLOT")
  (:export
   "<REFINEMENT-INDICATOR>" "INDICATE" "<UNIFORM-REFINEMENT-INDICATOR>"
   "<LARGEST-ETA-INDICATOR>" "<ABOVE-THRESHOLD-INDICATOR>"
   "ESTIMATE" "<PROJECTION-ERROR-ESTIMATOR>" "<DUALITY-ERROR-ESTIMATOR>" 
   "SOLVE-WITH" "<STRATEGY>" "INITIAL-GUESS" "SUFFICIENT-P"
   "IMPROVE-GUESS" "<FE-STRATEGY>" "LINEAR-PROBLEM-FE-STRATEGY"
   "*STRATEGY-OUTPUT*" "STRATEGY-OBSERVE"
   "*STANDARD-STRATEGY-OBSERVE-QUANTITIES*"
   "STOP-IF" "GLOBAL-ESTIMATE-SMALLER-THAN" "NR-OF-LEVELS>="))
