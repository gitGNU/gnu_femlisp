;;; -*- mode: lisp; fill-column: 64; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; discretization-demos.lisp
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

(in-package :fl.application)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Lagrange basis plots
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun 1d-graph (func &key (left 0.0) (right 1.0) (step 0.01))
  (loop for x from left upto right by step
	collect (vector x (evaluate func x))))

(defun plot-lagrange-basis (order type)
  "1D-~type~ - Plots 1D Lagrange basis (points=~type~)"
  (plot
   (loop for phi in (fe-basis (get-fe (lagrange-fe order :type type) *unit-interval*))
	 and i from 1
	 collect
	 (cons (format nil "phi-~D" i)
	       (1d-graph phi :left 0.0 :right 1.0 :step 0.01)))))

(defun make-lagrange-basis-demo (type)
  (multiple-value-bind (title short long)
      (extract-demo-strings
       (documentation 'plot-lagrange-basis 'function)
       `(("~type~" . ,(symbol-name type))))
    (let ((demo (make-demo
		 :name title :short short :long long
		 :execute
		 (lambda ()
		   (plot-lagrange-basis 
		    (user-input "Order (1-9): "
				#'(lambda (x) (and (integerp x) (plusp x))))
		    type))
		 :test-input (format nil "3~%"))))
      (adjoin-demo demo *discretization-demo*))))

(make-lagrange-basis-demo :uniform)
(make-lagrange-basis-demo :gauss-lobatto)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Stiffness matrices
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun model-problem-stiffness-matrix (dim level order)
  "Returns the stiffness matrix for the @arg{dim}-dimensional
model problem discretized on level @arg{level} with finite
elements of order @arg{order}."
  (let* ((problem (cdr-model-problem dim))
	 (h-mesh (uniformly-refined-hierarchical-mesh (domain problem) level))
	 (fedisc (lagrange-fe order)))
    (discretize-globally problem h-mesh fedisc)))

(defun display-stiffness-matrix (&rest args)
  "Shows the stiffness matrix of the model problem.  Maybe you
will have to use hscroll-mode in your Emacs buffer for
comfortably looking at the matrix."
  (let ((mat (apply 'model-problem-stiffness-matrix args)))
    (display mat :order (sort-lexicographically (row-keys mat)))))

(defun make-stiffness-matrix-demo (dim max-level order)
  (let ((demo
	 (make-demo
	  :name (format nil "~DD-stiffness-mat" dim)
	  :short (format nil "Stiffness matrix for Delta u = 0 in ~Dd" dim)
	  :long (format nil "~A~%" (documentation 'display-stiffness-matrix 'function))
	  :execute
	  (lambda ()
	    (dotimes (i max-level)
	      (format t "~&~%Level = ~D~%~%" i)
	      (display-stiffness-matrix dim i order)
	      (sleep 1.5))))))
    (adjoin-demo demo *discretization-demo*)))

(make-stiffness-matrix-demo 1 4 1)
(make-stiffness-matrix-demo 2 3 1)

