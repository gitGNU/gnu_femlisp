;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; iteration-defp.lisp
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

(defpackage "ITERATION"
  (:use "COMMON-LISP" "FL.MACROS" "FL.UTILITIES" "FL.MATLISP"
	"FL.TESTS" "FL.DEBUG"
	"ALGEBRA" "FL.FUNCTION")
  (:export
   ;; function.lisp
   "<FUNCTION>" "EVALUATE" "FUNC-DOMAIN-DIMENSION" "FUNC-IMAGE-DIMENSION"
   "<DIFFERENTIABLE-FUNCTION>" "EVALUATE-GRADIENT" "INTERVAL-METHOD"
   ;; solve.lisp
   ;; newton.lisp
   "NEWTON"
   ;; iterate.lisp
   "<ITERATION>" "*TIME-OBSERVE*"
   "INITIALLY" "INTERMEDIATE" "TERMINATE-P" "FINALLY" "OBSERVE"
   "SUCCESS-IF" "FAILURE-IF"
   "NEXT-STEP" "ITERATE" "OUTPUT"
   ;; solve.lisp
   "SOLVE" "<SOLVER>" "*ITERATIVE-SOLVER-OBSERVE*" "<ITERATIVE-SOLVER>"
   "INITIAL-GUESS" "ENSURE-RESIDUAL" "<SPECIAL-SOLVER>"
   ;; linit.lisp
   "<ITERATION>" "<LINEAR-ITERATION>" "DAMPING-FACTOR" "COMPUTE-RESIDUAL"
   "<ITERATOR>" "MAKE-ITERATOR"
   "<MULTI-ITERATION>" "PRODUCT-ITERATOR"
   "<LU>" "*LU-ITERATION*" "<ILU>" "*STANDARD-ILU*"
   "<JACOBI>" "<BLOCK-JACOBI>" "*UNDAMPED-JACOBI*"
   "<SOR>" "<GAUSS-SEIDEL>" "*GAUSS-SEIDEL*"
   "<SOLVER-ITERATION>"
   ;; linsolve.lisp
   "LINSOLVE" "<GPS-LINEAR-SOLVER>" "SELECT-LINEAR-SOLVER"
   ;; krylow.lisp
   "<GRADIENT-METHOD>" "<CG>" "<PCG>" "*STANDARD-CG*" "*STANDARD-PCG*"
   ;; blockit.lisp
   "<BLOCK-ITERATION>" "SETUP-BLOCKS" "MAKE-BLOCK" "BLOCKS"
   "<PSC>" "<SSC>" "<CUSTOM-PSC>" "<CUSTOM-SSC>"
   "<BLOCK-JACOBI>"
   "<BLOCK-SOR>" "<BLOCK-GAUSS-SEIDEL>" "<CUSTOM-BLOCK-GAUSS-SEIDEL>"
   ;; linsolve.lisp
   "<LINEAR-SOLVER>" "*LU-SOLVER*"
   ))

