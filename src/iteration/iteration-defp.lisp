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

(defpackage "FL.ITERATION"
  (:use "COMMON-LISP" "FL.MACROS" "FL.UTILITIES" "FL.MATLISP"
	"FL.TESTS" "FL.DEBUG" "FL.ALGEBRA" "FL.FUNCTION" "FL.PROBLEM")
  (:export  ; function.lisp
   "<FUNCTION>" "EVALUATE" "FUNC-DOMAIN-DIMENSION" "FUNC-IMAGE-DIMENSION"
   "<DIFFERENTIABLE-FUNCTION>" "EVALUATE-GRADIENT" "INTERVAL-METHOD")
  (:export  ; iterate.lisp
   "<ITERATION>" "*TIME-OBSERVE*"
   "INITIALLY" "INTERMEDIATE" "TERMINATE-P" "FINALLY" "OBSERVE"
   "SUCCESS-IF" "FAILURE-IF"
   "NEXT-STEP" "ITERATE" "OUTPUT")
  (:export  ; solve.lisp
   "*ITERATIVE-SOLVER-OBSERVE*" "<ITERATIVE-SOLVER>"
   "<SPECIAL-SOLVER>")
  (:export  ; linit.lisp
   "<LINEAR-ITERATION>" "DAMPING-FACTOR" "COMPUTE-RESIDUAL"
   "<ITERATOR>" "MAKE-ITERATOR"
   "<MULTI-ITERATION>" "PRODUCT-ITERATOR"
   "<LU>" "*LU-ITERATION*" "<ILU>" "*STANDARD-ILU*"
   "<JACOBI>" "<BLOCK-JACOBI>" "*UNDAMPED-JACOBI*"
   "<SOR>" "<GAUSS-SEIDEL>" "*GAUSS-SEIDEL*"
   "<SOLVER-ITERATION>")
  (:export  ; linsolve.lisp
   "<LINEAR-SOLVER>" "*LU-SOLVER*"
   "LINSOLVE")
  (:export  ; nlsolve.lisp
   "<NEWTON>")
  (:export  ; krylow.lisp
   "<GRADIENT-METHOD>" "<CG>" "<PCG>" "*STANDARD-CG*" "*STANDARD-PCG*")
  (:export  ; blockit.lisp
   "<BLOCK-ITERATION>" "SETUP-BLOCKS" "MAKE-BLOCK" "BLOCKS"
   "<PSC>" "<SSC>" "<CUSTOM-PSC>" "<CUSTOM-SSC>"
   "<BLOCK-JACOBI>"
   "<BLOCK-SOR>" "<BLOCK-GAUSS-SEIDEL>" "<CUSTOM-BLOCK-GAUSS-SEIDEL>")
  (:export  ; newton.lisp
   "NEWTON")
  (:documentation "The @package{FL.ITERATION} package includes the
definition for the abstract classes @class{<solver>},
@class{<iterative-solver>}, as well as the generic functions
@function{iterate} and @function{solve} which constitutes the interface for
linear and non-linear solving.  Both functions work on a blackboard which
is passed together with the iteration used as argument.

Several instances of iterative solvers are implemented, e.g. Gauss-Seidel,
SOR, ILU (in @file{linit.lisp}) and CG (in @file{krylow.lisp}).  A larger
block of code is contained in a separate package @package{FL.MULTIGRID} and
contains the multigrid iteration.  From this class, an algebraic multigrid
iteration is derived in @file{amg.lisp} and a geometric multigrid iteration
in a separate file @file{geomg.lisp} and package @package{FL.GEOMG}."))



