;;; -*- mode: lisp; fill-column: 64; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; cdr-bl.lisp
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

;;; In this file we compute effective boundary-layer constants for
;;; oscillating boundaries.  This situation was done theoretically by
;;; [Jaeger&Mikelic_1995] and by [Allaire&Amar_1999].  A further
;;; presentation was given in [Neuss_2002].

(defun boundary-layer-cell-problem (bl-cell-domain &key &allow-other-keys)
  "Returns a boundary layer problem for the given boundary layer
cell domain.  We assume that the oscillating boundary on which
Dirichlet nodes are prescribed is identified by $x_n<1$ for its
coordinates, and the patch for which the source is prescribed is
$x_n=0$."
  (let ((dim (dimension bl-cell-domain)))
    (make-instance
     '<cdr-problem> :domain bl-cell-domain
     :patch->coefficients
     #'(lambda (patch)
	 (cond ((and (= (dimension patch) (1- dim))
		     (bl-patch-on-artificial-boundary bl-cell-domain patch))
		(list 'CDR::SOURCE *cf-constantly-1.0d0*))
	       ((bl-patch-on-lower-boundary bl-cell-domain patch)
		(list 'CDR::DIRICHLET *cf-constantly-0.0d0*))
	       ((= (dimension patch) dim)
		(list 'CDR::DIFFUSION (identity-diffusion-tensor dim)))
	       (t nil))))))

(defun sinusoidal-boundary-layer-cell-problem (dim &rest key-args)
  "Returns a boundary layer problem for a sinusoidally
oscillating domain."
  (boundary-layer-cell-problem
   (apply #'sinusoidal-bl-cell dim key-args)))

(defun compute-cbl (result)
  (dot (getf result :solution) (getf result :rhs)))

(defparameter *result* nil)

(defun bl-computation (&key domain order max-levels output plot &allow-other-keys)
  "Computes a boundary layer and a boundary law coefficient
constant for a given boundary layer cell domain."
  (setf (getf *strategy-output* :observe)
	`((:Cbl "                Cbl"
	   ,#'(lambda (assembly-line)
		(with-items (&key solution rhs) assembly-line
		  (format t "~19,10,2E" (dot solution rhs)))))
	  (:eta "         ETA"
	   ,#'(lambda (assembly-line)
		(with-items (&key global-eta) assembly-line
		  (format t "~12,2,2E" global-eta))))))
  (setf (getf *strategy-output* :plot-mesh) t)
  (defparameter *result*
    (let* ((dim (dimension domain))
	   (problem (boundary-layer-cell-problem domain))
	   (solver
	    (?2 *lu-solver*
		(make-instance
		 '<linear-solver>
		 :iteration
		 (let ((smoother
			(if (> dim 2)
			    *gauss-seidel*
			    (make-instance '<local-bgs> :type :vertex-centered))))
		   (make-instance '<s1-reduction> :max-depth 2 :pre-steps 1 :pre-smooth smoother
				  :post-steps 1 :post-smooth smoother
				  :gamma 2 :coarse-grid-iteration
				  (?2 *lu-iteration*
				      (make-instance '<s1-coarse-grid-iterator>
						     :output (eq output :all)))))
		 :success-if `(:defnorm< 1.0e-10)
		 :failure-if `(and (:step> 2) (:step-reduction> 0.9))
		 :output (eq output :all)))))
      (solve-with
       (make-instance
	'<fe-strategy> :fe-class (lagrange-fe order)
	:solver solver
	:estimator (make-instance '<duality-error-estimator> :functional :load-functional)
	:indicator (make-instance '<largest-eta-indicator> :pivot-factor 0.01
				  :from-level 1 :block-p t)
	:appraise (stop-if :eta<= 1.0e-10 :nr-levels>= max-levels))
       problem
       :output t)))
  (when plot
    (plot (getf *result* :solution))))


(defun cdr-bl-computation (dim order max-levels &rest rest)
  "bl-diffusion-~Dd - Computes a boundary law coefficient

Computes a boundary law coefficient constant for an oscillating
function which has to be 1-periodic and negative.  This involves
solving a cell problem on a semi-infinite domain with a smooth
and exponentially decaying solution.  The constant we search for
is equal to the energy of the solution.  Our adaptive strategy
with high order approximations and local mesh refinement is
perfectly suited for computing this constant with high accuracy.
Error estimation is achieved by using a duality error estimator
for the load functional."
  (apply #'bl-computation
	 :domain (apply #'sinusoidal-bl-cell dim rest)
	 :order order :max-levels max-levels rest))

#+(or) (cdr-bl-computation
	2 4 2 :plot t :amplitude 0.15 :extensible-p nil :output :all)

#+(or) (show (getf *result* :matrix))
;; before change (+ 1 sec on another run)
;;    2        36       784    1.8   9.3593619483d-01   6.6484919470d-04
;;    8       172      4462   10.9   9.4073588298d-01   1.6845432818d-04
;;   20       444     12426   34.6   9.4064954875d-01   2.4130130651d-06


(defun make-cdr-bl-demo (dim order levels)
  (multiple-value-bind (title short long)
      (extract-demo-strings (documentation 'cdr-bl-computation 'function))
    (let ((demo
	   (make-demo
	    :name (format nil title dim) :short short
	    :long (format nil "~A~%~%Parameters: dim=~D, p=~D, ~D levels~%"
			  long dim order levels)
	    :execute (lambda () (cdr-bl-computation dim order levels :plot t)))))
      (adjoin-demo demo *laplace-demo*)
      (adjoin-demo demo *adaptivity-demo*)
      (adjoin-demo demo *boundary-coeffs-demo*))))

;;; 2D and 3D demo
(make-cdr-bl-demo 2 4 3)
(make-cdr-bl-demo 3 3 2)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Testing
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun test-cdr-bl ()
  
  (bl-computation :domain (spline-interpolated-bl-cell #(1.0d0 0.8d0))
		  :order 4 :max-levels 3 :plot t)
  
  ;;; testing if identification with only one cell width works
  (let* ((problem (sinusoidal-boundary-layer-cell-problem 2))
	 (mesh (make-hierarchical-mesh-from-domain (domain problem)))
	 (fe-class (lagrange-fe 1)))
    (refine mesh :test #'(lambda (cell) (<= (aref (midpoint cell) 1) 0.0)))
    (multiple-value-bind (mat rhs)
	(discretize-globally problem mesh fe-class)
      #+(or)(show mat)
      #-(or)
      (plot (linsolve mat rhs :iteration *lu-iteration*)
	    :depth 3)
      ))

  (plot (getf *result* :mesh))
  (plot (strategy::eta->p2-vec (getf *result* :eta) (getf *result* :problem)
			       (getf *result* :mesh)))
  (display-ht (getf *result* :eta))
  (apply #'strategy::compute-local-estimate (strategy::estimator (getf *result* :strategy)) *result*)

  (dohash (key (getf *result* :eta))
    (let ((mp (midpoint key)))
      (when (< (norm (vec- mp #(0.9375 -0.3125))) 1.0e-10)
	(display-ht (matrix-row (getf *result* :matrix) key)))))
 
  (plot (getf *result* :solution))

  ;; reiteration
  (setf (getf *result* :solution)
	(destructuring-bind (&key matrix solution rhs &allow-other-keys)
	    *result*
	  (solve (s1-reduction-amg-solver 4 :reduction 1.0e-3 :output t)
		 :matrix matrix :rhs rhs :solution solution)))
  (plot *result*)
  (plot (getf *result* :mesh))
  (getf *result* :global-eta)
  (print-svec (getf *result* :solution))
  (plot (getf *result* :solution) :depth 2 :plot :file :filename "bl-cell-sol" :format "eps")
  (plot (mesh (getf *result* :solution)) :plot :file :filename "bl-cell-mesh.ps" :format "eps")
  (compute-cbl *result*)
  (nr-of-levels (mesh (getf *result* :solution)))
  (nr-of-cells (mesh (getf *result* :solution)))
  (destructuring-bind (&key matrix solution rhs &allow-other-keys)
      *result*
    (plot (m- rhs (m* matrix solution))))

;O=1
;1.0
;0.9994139684729783d0
;0.9478047281289215d0  2.34
;0.9427405097092483d0 10.9
;0.9412018564593625d0 56
;0.940791738100939d0 650 (cells 11300)

;O=4:
;0.9341601270140.
;0.940503094101
;0.940644213164

;O=5
;0.9414553049022709d0 2.45
;0.9406756579556297d0 8.96
;0.9406523275905488d0 27.3
;0.9406521294202064d0 87
;0.9406521488254819d0 550

  (problem-info (sinusoidal-boundary-layer-cell-problem 2))

  )
