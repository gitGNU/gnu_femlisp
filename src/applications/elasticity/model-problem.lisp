;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; model-problem.lisp - Model problems for linear elasticity
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

(in-package :application)

(defun test-elasticity-model-problems ()

;;; Linear elasticity problem
(defparameter *result*
  (time
   (let* ((dim 2) (level 2) (order 3)
	  (problem
	   (standard-elasticity-problem
	    (n-cube-domain dim) :lambda 1.0 :mu 1.0
	    :force (constantly (unit-vector dim 0))))
	  (mm (uniformly-refined-hierarchical-mesh (domain problem) level))
	  (fe-class (lagrange-fe order :nr-comps dim)))
     (multiple-value-bind (mat rhs)
	 (discretize-globally problem mm fe-class)
       #+(or)(show mat)
       #+(or)(show rhs)
       #+(or)(m* (sparse-ldu mat) rhs)
       #-(or)(linsolve mat rhs :output t :iteration (geometric-cs :fmg t) :maxsteps 10)
       ))))
(plot *result* :component 1)
;;; Test for diffusion problem

(defparameter *result*
  (time
   (let* ((dim 2) (nr-comps 1) (level 2) (order 1)
	  (problem
	   (system-diffusion-problem
	    (n-cube-domain dim) :D 1.0 :nr-comps nr-comps
	    :force (constantly (make-array nr-comps :initial-element (ones 1)))))
	  (mm (uniformly-refined-hierarchical-mesh (domain problem) level))
	  (fe-class (lagrange-fe order :nr-comps nr-comps)))
     (multiple-value-bind (mat rhs)
	 (discretize-globally problem mm fe-class)
       ;;(show mat)
       #+(or)(m* (sparse-ldu mat) rhs)
       (linsolve mat rhs :output t :iteration (geometric-cs :fmg t) :maxsteps 2)
       ))))
(plot *result* :component 0)

;;; comparison with scalar case
(defparameter *result*
  (time
   (let* ((dim 2) (level 2) (order 1)
	  (problem (cdr-model-problem dim))
	  (mm (uniformly-refined-hierarchical-mesh (domain problem) level))
	  (fe-class (lagrange-fe order)))
     (multiple-value-bind (mat rhs)
	 (discretize-globally problem mm fe-class)
       #+(or)(m* (sparse-ldu mat) rhs)
       (linsolve mat rhs :output t :iteration (geometric-cs :fmg t) :maxsteps 2)
       ))))

)
