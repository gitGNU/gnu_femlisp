;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;  elasticity.lisp
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

(in-package :CL-USER)
(defpackage "FL.ELASTICITY-FE"
  (:use "COMMON-LISP" "FL.MATLISP" "FL.MACROS" "FL.UTILITIES"
	"FL.FUNCTION"
	"FL.MESH"
	"FL.PROBLEM" "FL.ELLSYS" "FL.ELASTICITY"
	"FL.DISCRETIZATION")
  (:export)
  (:documentation "This package specializes the discretization for systems
of elasticity.  Since elasticity is a special case of elliptic systems
which are handled in @path{ellsys-fe.lisp}, not much remains to do."))
(in-package "FL.ELASTICITY-FE")

(defun discretize-elasticity (dim order level)
  "Local routine for testing purposes."
  (let* ((problem (elasticity-model-problem
		   (n-cube-domain dim) :lambda 1.0 :mu 1.0
		   :force (constantly
			   (coerce (loop repeat dim collect (eye 1)) 'vector))))
	 (h-mesh (uniformly-refined-hierarchical-mesh (domain problem) level))
	 (fedisc (lagrange-fe order :nr-comps dim)))
    (discretize-globally problem h-mesh fedisc)))

;;; Testing
(defun elasticity-fe-tests ()
  
  (multiple-value-bind (mat rhs)
      (discretize-elasticity 1 1 1)
    (getrs (sparse-ldu mat) rhs))

  #+(or)
  (multiple-value-bind (matrix rhs)
      (time (discretize-elasticity 2 1 5))
    (fl.plot:plot (gesv matrix rhs) :component 0))
  
  )

;;; (elasticity-fe-tests)
(fl.tests:adjoin-test 'elasticity-fe-tests)
