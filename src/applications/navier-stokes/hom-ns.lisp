;;; -*- mode: lisp; fill-column: 64; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; hom-ns.lisp - Computing effective coefficients for Navier-Stokes
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

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Cell problem for Navier-Stokes in a porous cell
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun ns-cell-problem-force (dim)
  (constant-coefficient
   (coerce (loop for i upto dim collect
		 (if (< i dim)
		     (ensure-matlisp (unit-vector dim i) :row)
		     (zeros 1 dim)))
	   'vector)))

(defun ns-hole-cell-problem (dim &key (viscosity 1.0) (reynolds 0.0))
  (let* ((domain (n-cell-with-n-ball-hole dim)))
    (make-instance
     '<navier-stokes-problem>
     :domain domain
     :patch->coefficients
     #'(lambda (patch)
	 (cond ((patch-on-inner-boundary-p patch)
		(list 'CONSTRAINT (no-slip-boundary dim)))
	       ((= dim (dimension patch))  ; inner coeffs
		(list 'NAVIER-STOKES::VISCOSITY (constant-coefficient viscosity)
		      'NAVIER-STOKES::REYNOLDS (constant-coefficient reynolds)
		      'NAVIER-STOKES::FORCE (ns-cell-problem-force dim)))
	       (t ())))
     :multiplicity dim)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Computation of the permeability tensor
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun permeability-tensor (blackboard)
  (with-items (&key solution rhs) blackboard
    (correction-tensor solution rhs)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Demos
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun stokes-darcy-demo (problem &key order levels plot output (store-p t) (delta 1))
  "Stokes-Darcy - Computes a Darcy permeability tensor

Computes the effective permeability for Stokes flow in a domain
with periodically distributed ball-shaped holes.  The
approximation is done with finite elements and blending.
Uniform refinement is used, the linear solver is a geometric
multigrid scheme used in a nested iteration fashion.  Smoothing
is done by a Vanka type smoother.

The solution to this cell problem in ~d~ dimensions consists of
~d~ velocity/pressure pairs, i.e. in total ~d~*(~d~+1)
components."
  (let* ((domain (domain problem))
	 (dim (dimension domain)))
    (defparameter *result*
      (solve
       (make-instance
	'<stationary-fe-strategy>
	:fe-class (navier-stokes-lagrange-fe order dim delta)
	:estimator (make-instance '<projection-error-estimator>)
	:indicator (make-instance '<largest-eta-indicator> :fraction 1.0)
	:success-if `(>= :nr-levels ,levels)
	:solver
	(make-instance
	 '<linear-solver> :iteration
	 (let ((smoother (make-instance '<vanka> :store-p store-p)))
	   (make-instance
	    '<geometric-cs>
	    :coarse-grid-iteration
	    (make-instance '<multi-iteration> :nr-steps 3 :base smoother)
	    :pre-steps 1 :pre-smooth smoother
	    :post-steps 1 :post-smooth smoother
	    :gamma 2))
	 :success-if `(and (> :step 2) (> :step-reduction 0.9) (< :defnorm 1.0e-9))
	 :failure-if `(and (> :step 2) (> :step-reduction 0.9) (>= :defnorm 1.0e-9))
	 :output (eq output :all))
	:output t :observe
	(append *stationary-fe-strategy-observe*
		(list
		 (list (format nil "~19@A~19@A~19@A" "K_00" "K_01" "K_11") "~57A"
		      #'(lambda (blackboard)
			  (let ((tensor (permeability-tensor blackboard)))
			    (format t "~19,10,2E~19,10,2E~19,10,2E"
				    (and tensor (matrix-ref tensor 0 0))
				    (and tensor (matrix-ref tensor 0 1))
				    (and tensor (matrix-ref tensor 1 1)))))))))
       (blackboard :problem problem :output t)))
    (when plot
      ;; plot components of cell solution tensor
      (dotimes (i dim)
	(dotimes (j (1+ dim))
	  (plot (getbb *result* :solution) :component j :index i)
	  (sleep 1.0))))
    ;; compute the homogenized coefficient
    (format t "The permeability tensor is:~%~A~%"
	    (permeability-tensor *result*))))

#+(or)
(stokes-darcy-demo (ns-hole-cell-problem 2)
		   :order 4 :levels 2 :output :all :plot t :delta 2)
#+(or)(plot (getbb *result* :solution) :component 0 :index 0)

#+(or)
(stokes-darcy-demo (ns-hole-cell-problem 3)
		   :order 1 :levels 3 :output :all :plot nil :store-p nil)

#+(or)
(let ((dim 2)
      (counter -1))
  (dotimes (i dim)
    (dotimes (j (1+ dim))
      (plot (getbb *result* :solution) :component j :index i :depth 3
	    :plot :file :format "tiff" :background :white
	    :filename (format nil "hom-ns-cell-sol-~D" (incf counter))))))
;; montage -geometry 480x480 -tile 2x3 hom-ns-cell-*.tiff hom-ns-cell-sols.eps
;; only velocities: mogrify -crop 1380x900+30+30 hom-ns-cell-sols.eps


(defparameter *effective-permeability-demo*
  (make-demo
   :name "effective-permeability"
   :short "Computes a permeability tensor"
   :long (documentation 'stokes-darcy-demo 'function)))

(adjoin-demo *effective-permeability-demo* *interior-coeffs-demo*)
(adjoin-demo *effective-permeability-demo* *navier-stokes-demo*)


(defun make-effective-permeability-demo (dim order levels)
  (let* ((short (format nil "Computes a permeability tensor in ~DD" dim))
	 (demo (make-demo
		:name (format nil "~DD-cell" dim)
		:short short
		:long (format nil "We compute the effective
permeability tensor for a medium with periodically distributed
holes by solving a cell problem on a periodicity cell.~%
Parameters: order=~D, max-levels=~D~%~%"
			      order levels)
		:execute
		(lambda ()
		  (stokes-darcy-demo (ns-hole-cell-problem dim)
				     :order order :levels levels :plot t)))))
    (adjoin-demo demo *effective-permeability-demo*)))

(make-effective-permeability-demo 2 4 2)
(make-effective-permeability-demo 3 2 1)

(defun normalize-pressure (svec &key (normalize t))
  "This stuff is apparently not needed for good convergence."
  (let* ((pressure (component svec 2))
	 (mult (multiplicity pressure))
	 (sum (make-double-vec mult))
	 (count (make-fixnum-vec mult)))
    ;; compute average
    (for-each-entry
     #'(lambda (entry)
	 (for-each-key-and-entry
	  #'(lambda (i j entry)
	      (declare (ignore i))
	      (incf (aref sum j) entry)
	      (incf (aref count j)))
	  entry))
     pressure)
    (if normalize
	(for-each-entry
	 #'(lambda (entry)
	     (for-each-key
	      #'(lambda (i j)
		  (decf (matrix-ref entry i j)
			(/ (aref sum j) (aref count j))))
	      entry))
	 pressure)
	(map 'vector #'/ sum count))))

;;; Testing: (hom-ns-tests)
(tests::adjoin-femlisp-test 'hom-ns-tests)

(defun hom-ns-tests ()
  (stokes-darcy-demo
   (ns-hole-cell-problem 2)
   :order 4 :levels 2 :output :all :plot nil :delta 2))