;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; multigrid-defp.lisp
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

(defpackage "MULTIGRID"
  (:use "COMMON-LISP" "GENERAL" "MACROS" "UTILITIES" "ALGEBRA" "MATLISP"
	"TESTS" "ITERATIONS" "FEMLISP-DEBUG")
  (:export				; multigrid.lisp
   "<MG-ITERATION>" "PRE-SMOOTH" "PRE-STEPS" "POST-SMOOTH" "POST-STEPS"
   "BASE-LEVEL" "COARSE-GRID-ITERATION" "FMG"
   "<CORRECTION-SCHEME>" "<FAS>"
   "*DEFAULT-SMOOTHER*" "*DEFAULT-COARSE-GRID-ITERATION*"
   "MULTILEVEL-DECOMPOSITION")
  (:export				; amg.lisp
   "SLAVE-DOF-P" "DIRICHLET-DOF-P" "SLAVE-OR-DIRICHLET-DOF-P" 
   "COARSEN" "PREPROCESS-MATRIX" "PROLONGATION" "RESTRICTION"
   "COARSE-GRID-MATRIX")
  (:export				; selection-amg.lisp
   "<SELECTION-AMG>")
  (:export				; cluster-amg.lisp
   "<CLUSTER-AMG>")
  (:export				; stueben.lisp
   "<STUEBEN>" "*STANDARD-STUEBEN*")
  )

