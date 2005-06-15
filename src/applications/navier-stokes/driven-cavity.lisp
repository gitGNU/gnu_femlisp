;;; -*- mode: lisp; fill-column: 64; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; driven-cavity.lisp - Driven cavity computations
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

(defun watch-dc-center-velocity (dim)
  "Returns a observe list for watching the velocity at the center of the
driven cavity."
  (list (format nil "~{                 u~1D~}" (range<= 1 dim))
	"~{~19,10,2E~}"
	#'(lambda (blackboard)
	    (with-items (&key solution) blackboard
	      (let ((val (fe-value solution (make-double-vec dim 0.5))))
		(loop for i below dim collect (vref (aref val i) 0)))))))

(defun ns-driven-cavity-demo (dim order levels &key plot output (reynolds 0.0) smooth-p)
  "Performs the driven cavity demo."
  (defparameter *result*
    (solve 
     (blackboard
      :fe-class (navier-stokes-lagrange-fe order dim 1)
      :problem (driven-cavity dim :reynolds reynolds :smooth-p smooth-p)
      :base-level (if (> order 1) 0 1)
      :success-if (if levels
		      `(= :nr-levels ,levels)
		      `(> :time ,*demo-time*))
      :output output :observe
      (append *stationary-fe-strategy-observe*
	      (list (watch-dc-center-velocity dim))))))
  (when plot
    (let ((solution (getbb *result* :solution)))
      ;; plot components of cell solution tensor
      (dotimes (i (1+ dim))
	(plot solution :component i :depth 2)
	(sleep 1.0)))))

;;; (ns-driven-cavity-demo 2 2 5 :output :all :plot nil :reynolds 100.0)
;;; (plot (getbb *result* :solution) :component 0 :depth 2)

(defun make-driven-cavity-demo (dim order reynolds)
  (let ((title (format nil "DC-~DD-~A" dim reynolds))
	(short (format nil "Solves the driven cavity problem (Re=~A)." reynolds))
	(long (format nil "Solve the ~DD driven cavity problem
for the Navier-Stokes equation using Taylor-Hood finite elements
[Q^~D]^~D-Q^~D." dim (1+ order) dim order)))
    (let ((demo
	   (make-demo
	    :name title :short short :long long :execute
	    (lambda ()
	      (ns-driven-cavity-demo dim order 3 :output 1 :plot t
				     :reynolds (float reynolds 1.0))))))
      (adjoin-demo demo *navier-stokes-demo*))))

;;(ns-driven-cavity-demo 2 4 3 :output :all :plot t :reynolds 0.0)
(make-driven-cavity-demo 2 2 0)
(make-driven-cavity-demo 2 2 100)

;;;; Testing:

(defun test-driven-cavity ()
  (describe (driven-cavity 2))
  (describe (domain (driven-cavity 2 :smooth-p nil)))
  (ns-driven-cavity-demo 2 1 4 :output :all :plot nil)
  (let ((sol (getbb *result* :solution)))
    (fe-value sol #d(0.5 0.5)))
  #+(or)
  (defparameter *result*
    (let* ((dim 2) (order 3) (delta 1)
	   (problem (driven-cavity dim :smooth-p nil))
	   (mesh
	    (change-class
	     (triangulate (domain problem) :meshsize 0.01 :indicator
			  #'(lambda (patch x h)
			      (declare (ignore patch))
			      (let ((d (min (norm (m- x #d(0.0 1.0)))
					    (norm (m- x #d(1.0 1.0))))))
				(cond
				  ((>= h 0.25) :yes)
				  ((<= h (* 0.5 d)) :no)))))
	     '<hierarchical-mesh>))
	   (as (make-fe-ansatz-space (navier-stokes-lagrange-fe order dim delta)
				     problem mesh)))
      (solve (blackboard :problem problem :mesh mesh :ansatz-space as
			 :output :all :success-if '(> :time 30.0) :observe
			 (append *stationary-fe-strategy-observe*
				 (list (watch-dc-center-velocity dim)))))))
  (fe-extreme-values (getbb *result* :solution))
  (time (plot (getbb *result* :solution) :component 1))
  (time (plot (component (getbb *result* :solution) 1)))
)

;;; (test-driven-cavity)
(fl.tests:adjoin-test 'test-driven-cavity)


