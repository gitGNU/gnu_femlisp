;;; -*- mode: lisp; fill-column: 64; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; effcoeff.lisp
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

;;; Demo for a talk about effective coefficients which was given
;;; at Oberwolfach in June 2003 and the presentation of my
;;; habilitation in July 2003.

(defvar *effcoeff-root*
  (make-demo
   :name "effective-coefficients-demo"
   :short "Calculation of effective coefficients"
   :long "Femlisp demos for the talk about computing effective
coefficients."))

(defun smoother-performance-test (&key dim order (level 2) smoother output)
  "Tests performance of smoother on a Laplace model problem.
See make-smoother-demo for more information."
  (let* ((problem (laplace-test-problem-on-domain (n-cube-domain dim)))
	 (mm (uniformly-refined-hierarchical-mesh
	      (domain problem) level))
	 (fe-class (lagrange-fe order)))
     (multiple-value-bind (mat rhs)
	 (discretize-globally problem mm fe-class)
       (getf (nth-value 1 (linsolve mat rhs :output output :iteration smoother
				    :maxsteps 10 :threshold 1.0e-10))
	     :last-step-reduction))))

#+(or)
(let ((dim 1) (order 6))
  (smoother-performance-test :dim dim :order order :smoother *gauss-seidel*))

(defun smoother-demo-execute (smoother)
  (lambda ()
    (loop for dim = (user-input "Dimension (1..4): "
				#'(lambda (x) (and (integerp x) (< 0 x 5))))
	  until (eq dim :up) do
	  (loop with max-order = (case dim (1 8) (2 7) (t 4))
		for order =
		(user-input (format nil "Order (1..~A): " max-order)
			    #'(lambda (x) (and (integerp x) (<= 1 x max-order))))
		until (eq order :up) do
		(smoother-performance-test
		 :dim dim :order order :smoother smoother :output t)))))

(defun make-smoother-demo (smoother smoother-name)
  "~name~-performance - Tests smoother ~name~

Tests the performance of ~name~ applied to discretizations of
different order of a Laplace model problem on cubes of different
dimensions."
  (multiple-value-bind (name short long)
      (extract-demo-strings
       (documentation 'make-smoother-demo 'function)
       (list (cons "~name~" smoother-name)))
    (let ((demo
	   (make-demo
	    :name name :short short :long long
	    :execute (smoother-demo-execute smoother))))
      (adjoin-demo demo *effcoeff-root*))))

#+(or)(make-smoother-demo *gauss-seidel* "GS")
#+(or)(make-smoother-demo (make-instance '<local-bgs> :type :vertex-centered) "VC-BGS")
#+(or)(make-smoother-demo (make-instance '<local-bgs> :type :cell-centered) "CC-BGS")

(defun smoother-graph-execute (smoother)
  (lambda ()
    (loop for dim = (user-input "Dimension (1..4): "
				#'(lambda (x) (and (integerp x) (< 0 x 5))))
	  until (eq dim :up) do
	  (plot
	   (list
	    (cons
	     "convergence-rate"
	     (loop for order from 1 upto (case dim (1 8) (2 4) (t 3))
		   collect
		   (vector order
			   (smoother-performance-test
			    :dim dim :order order :smoother smoother)))))
	   :left 0 :right 8 :top 1.0 :bottom 0.0
	   ))))

(defun make-smoother-performance-graph-demo (smoother smoother-name)
  "~name~-cr-graph - Plots graph 'order->CR(order)' for ~name~

Plots a graph of the convergence rate for ~name~ smoother
applied to discretizations of different order of a Laplace model
problem on cubes of different dimensions."
  (multiple-value-bind (name short long)
      (extract-demo-strings
       (documentation 'make-smoother-performance-graph-demo 'function)
       (list (cons "~name~" smoother-name)))
    (let ((demo
	   (make-demo
	    :name name :short short :long long
	    :execute (smoother-graph-execute smoother))))
      (adjoin-demo demo *effcoeff-root*))))

(make-smoother-performance-graph-demo *gauss-seidel* "GS")
(make-smoother-performance-graph-demo
 (make-instance '<local-bgs> :type :vertex-centered) "VC-BGS")

(let ((demo (find-demo "bl-diffusion-2d" *laplace-demo*)))
  (adjoin-demo demo *effcoeff-root*))

;;; finally also the standard demos can be accessed
(adjoin-demo *demo-root* *effcoeff-root*)


(loop for demo in (femlisp-demo::leaves *effcoeff-root*) do
      (remhash (femlisp-demo::name demo) *visited-demos*))

#+(or)(demo *effcoeff-root*)

