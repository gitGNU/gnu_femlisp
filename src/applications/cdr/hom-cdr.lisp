;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; hom-cdr.lisp - computing effective coefficients for cdr problems
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
;;; Utilities
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun cell-solve (problem &key (level 1) (order 1) parametric (solver (lu-solver)))
  "Solves the given problem and returns a result on a blackboard."
  (let* ((domain (domain problem))
	 (mesh (uniformly-refined-hierarchical-mesh domain level :parametric parametric))
	 (fe-class (lagrange-fe order))
	 (as (make-fe-ansatz-space fe-class problem mesh))
	 (blackboard (fe-discretize (blackboard :ansatz-space as))))
    (setf (getbb blackboard :problem) problem)
    (setf (getbb blackboard :mesh) mesh)
    (with-items (&key solution matrix rhs)
	blackboard
      (setq solution
	    (getbb (solve solver (blackboard :problem (lse :matrix matrix :rhs rhs)))
		   :solution))
      blackboard)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Cell problems
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; cell problems defined by the coefficient function on regular meshes

(defun cdr-cell-problem (dim/domain &key diffusion-function)
  "Returns the cell problem on the n-cell-domain of the given dimension.
Gamma yields dim right hand sides of the form A e_k, i.e. the columns of
the diffusion tensor."
  (multiple-value-bind (dim domain)
      (etypecase dim/domain
	(number (values dim/domain (n-cell-domain dim/domain)))
	(<domain> (values (dimension dim/domain) dim/domain)))
    (ensure diffusion-function (constantly (eye dim)))
    (make-instance
     '<cdr-problem> :domain domain
     :multiplicity dim :patch->coefficients
     #'(lambda (patch)
	 (cond
	   ((= (dimension patch) dim)
	    (list 'FL.CDR::DIFFUSION (ensure-coefficient diffusion-function)
		  'FL.CDR::GAMMA (constant-coefficient (eye dim))))
	   ((mzerop (midpoint patch))
	    (list 'FL.CDR::CONSTRAINT (constant-coefficient 0.0)))))
	 )))

(defun simple-square-inlay-cell-problem (dim)
  (cdr-cell-problem
   dim :diffusion-function
   #'(lambda (x) 
       (let ((result (eye dim)))
	 (when (every #'(lambda (coord) (<= 0.25 coord 0.75)) x)
	    (scal! 0.1 result))
	  result))))

(defun smooth-coefficient-cell-problem (dim)
  (cdr-cell-problem
   dim :diffusion-function
   #'(lambda (x)
       (scal (reduce #'* (map 'vector #'(lambda (xc) #I(2.0+sin(2*pi*xc))) x))
	     (eye dim)))))
  
(defun simple-ball-inlay-cell-problem (dim eps)
  (cdr-cell-problem
   dim :diffusion-function
   #'(lambda (x)
       (let ((result (eye dim)))
	 (when (<= (norm (m- (make-double-vec dim 0.5) x)) 0.25)
	   (scal! eps result))
	 result))))

(defun chequerboard-problem (dim eps)
  "Returns cell problem on the n-cell-domain of the given dimension with
a chequerboard pattern."
  (cdr-cell-problem
   dim :diffusion-function
   #'(lambda (x)
       (scal! (if (evenp (loop for coord across x count (>= coord 0.5)))
		  1.0 eps)
	      (eye dim)))))

;;; cell problem defined with inlay-adapted subdomain

(defun inlay-cell-problem (dim eps)
  "Generates the inlay cell problem.  Parameters are the coefficient of the
inlay and the component of the cell vector."
  (make-instance
   '<cdr-problem> :domain (n-cell-with-ball-inlay dim)
   :multiplicity dim :patch->coefficients
   #'(lambda (patch)
       (when (= (dimension patch) dim)
	 (list 'FL.CDR::DIFFUSION
	       (constant-coefficient
		(scal! (if (patch-in-inlay-p patch) eps 1.0) (eye dim)))
	  'FL.CDR::GAMMA (constant-coefficient (eye dim)))))))

(defun porous-cell-problem (dim &key (radius 0.25) A)
  (cdr-cell-problem
   (if A
       (n-cell-with-ellipsoidal-hole dim :A A)
       (n-cell-with-ball-hole dim :radius radius))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Demos
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun cdr-interior-effective-coeff-demo (problem order levels
					  &key (output 1) plot)
  "Computes an effective diffusion tensor for different configurations.
The approximation is done with isoparametric finite elements.  Uniform
refinement is used, the linear solver is a geometric multigrid scheme used
in a nested iteration fashion.  The solution to the cell problem is a
vector field whose components are plotted one after the other.  The setting
has cubic symmetry, from which one can easily see that the effective tensor
must be a scalar multiple of the identity."
  (let ((*output-depth* output))
    (setq *result*
      (solve
       (make-instance
	'<stationary-fe-strategy> :fe-class (lagrange-fe order)
	:estimator (make-instance '<projection-error-estimator>)
	:indicator (make-instance '<largest-eta-indicator> :fraction 1.0)
	:success-if `(>= :nr-levels ,levels)
	:solver
	(make-instance
	 '<linear-solver> :iteration
	 (let ((smoother
		(if (>= (domain-dimension problem) 3)
		    *gauss-seidel*
		    (geometric-ssc :store-p nil))))
	   (geometric-cs
	    :gamma 2 :fmg nil :coarse-grid-iteration
	    (make-instance '<multi-iteration> :base smoother
			   :nr-steps (if (eq smoother *gauss-seidel*) 10 3))
	    :smoother smoother :pre-steps 2 :post-steps 2))
	 :success-if `(and (> :step 2) (> :step-reduction 0.9) (< :defnorm 1.0e-10))
	 :failure-if `(> :step 20))
	:observe
	(append *stationary-fe-strategy-observe*
		(list (list (format nil "~19@A" "Ahom") "~19,10,2E"
			    #'(lambda (blackboard)
				(mref (effective-tensor blackboard) 0 0))))))
       (blackboard :problem problem)))
    ;; plot cell solutions and compute the homogenized coefficient
    (when plot
      (let ((solution (getbb *result* :solution)))
	(plot solution :index 0)
	(sleep 1.0)
	(plot solution :index 1)))
    (format t "The effective tensor is:~%~A~%"
	    (effective-tensor *result*))))

;;; Testing:
#+(or)
(let ((*output-depth* 2))
	(cdr-interior-effective-coeff-demo (porous-cell-problem 2) 4 2 :plot t))
#+(or)(cdr-interior-effective-coeff-demo (inlay-cell-problem 2 0.1) 4 2)
#+(or)
(let ((A (FL.algebra::ellipse-matrix 0.25 0.3 0.7854)))
  (cdr-interior-effective-coeff-demo (porous-cell-problem 2 :A A) 4 2 :plot t))

(defparameter *effective-diffusion-demo*
  (make-demo
   :name "effective-diffusion"
   :short "Computes an effective diffusion tensor"
   :long (documentation 'cdr-interior-effective-coeff-demo 'function)))

(adjoin-demo *effective-diffusion-demo* *interior-coeffs-demo*)
(adjoin-demo *effective-diffusion-demo* *cdr-demo*)

(defun make-effective-diffusion-porous-domain-demo (dim order levels)
  (let* ((short (format nil "Computes effective diffusion for a ~DD cell with hole" dim))
	 (demo (make-demo
		:name (format nil "~DD-hole-cell" dim)
		:short short
		:long (format nil "Parameters: order=~D, max-levels=~D~%~%"
			      order levels)
		:execute (lambda () (cdr-interior-effective-coeff-demo
				     (porous-cell-problem dim) order levels :plot t)))))
    (adjoin-demo demo *effective-diffusion-demo*)))

(make-effective-diffusion-porous-domain-demo 2 5 2)
(make-effective-diffusion-porous-domain-demo 3 4 1)

;;; 3D Ordnung 4 und 5
;;;    6      1287     73393      30.8   8.3975311514d-01
;;;    6      2463    236081      59.7   8.3943428490d-01

(defun make-effective-diffusion-inlay-domain-demo (dim order levels)
  (let* ((short (format nil "Computes effective diffusion for a ~DD cell with a soft inlay." dim))
	 (demo (make-demo
		:name (format nil "~DD-inlay-cell" dim) :short short
		:long (format nil "~A~%The inlay has diffusivity 0.1 the rest 1.0.
~%Parameters: order=~D, levels=~D~%~%"
			      short order levels)
		:execute (lambda () (cdr-interior-effective-coeff-demo
				     (inlay-cell-problem dim 0.1) order levels :plot t)))))
    (adjoin-demo demo *effective-diffusion-demo*)))

(make-effective-diffusion-inlay-domain-demo 2 4 3)
(make-effective-diffusion-inlay-domain-demo 3 1 1)

#+(or)(demo *effective-diffusion-demo*)
     
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Testing: (cdr-hom-tests)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun cdr-hom-tests ()

(cdr-interior-effective-coeff-demo (porous-cell-problem 2) 4 2 :output :all)

(cdr-interior-effective-coeff-demo (inlay-cell-problem 2 0.1) 1 3 :output :all)

;;grid is not inlay-adapted
(plot (simple-ball-inlay-cell-problem 1 0.1) :refinements 3
	:coefficient 'FL.CDR::DIFFUSION :key (rcurry #'mref 0 0))

;; solve cell problem and compute homogenized coefficient
(setq *result*
  (cell-solve (simple-ball-inlay-cell-problem 2 0.1)
	      :level 3 :order 1))
(effective-tensor *result*)
(plot (getbb *result* :solution) :index 0)

;; generate 3D cell problem
(simple-ball-inlay-cell-problem 3 0.1)

;; inlay adapted grid
(plot (inlay-cell-problem 2 0.1) :refinements 0 :depth 2 :parametric (lagrange-mapping 3)
	:coefficient 'FL.CDR::DIFFUSION :key (rcurry #'mref 0 0))

;; first order with lu-solver
(setq *result*
      (cell-solve (inlay-cell-problem 2 0.1) :level 1
		  :order 1 :parametric (lagrange-mapping 2)))
(effective-tensor *result*)
(plot (getbb *result* :solution) :index 0)

;; first order with full multigrid
(setq *result*
  (time
   (let ((level 2) (order 1))
     (cell-solve
      (inlay-cell-problem 2 0.1) :level level
      :order order :parametric (lagrange-mapping order)
      :solver				; (lu-solver) #+(or)
      (make-instance
       '<linear-solver> :iteration
       (geometric-cs
	:fmg t :coarse-grid-iteration
	(make-instance '<multi-iteration> :base *gauss-seidel* :nr-steps 3))
       :output t
       :success-if `(< :reduction 1.0d-5)
       :failure-if `(> :step 20)
       )))))
(plot (getbb *result* :solution) :index 0)
(plot (getbb *result* :solution) :index 1)

(setq *result*
  (cell-solve (inlay-cell-problem 2 0.1) :level 2
		:order 1 :parametric (lagrange-mapping 2)))
(plot (getbb *result* :mesh)) ; :plot :file :format "tiff")

;; higher order
(time
 (setq *result*
   (let ((dim 2) (level 1) (order 4))
     (cell-solve
      (inlay-cell-problem dim 0.1) :level level
      :order order :parametric (lagrange-mapping (max order 2))
      :solver ;(lu-solver) #+(or)
      (make-instance
       '<linear-solver> :iteration
       (let ((smoother *gauss-seidel*))
	 (geometric-cs
	  :gamma 2 :fmg t :coarse-grid-iteration
	  (make-instance '<multi-iteration> :base *gauss-seidel* :nr-steps 10)
	  :smoother smoother :pre-steps 2 :post-steps 2))
       :success-if `(and (> :step 2) (> :step-reduction 0.9) (< :defnorm 1.0e-10))
       :failure-if `(> :step 20))))))

;; (l=2, o=4, cgs=10 ILU, 3 smooth, fmg)*2=146.7 URT
(effective-tensor *result*)
(plot (getbb *result* :solution) :index 1)

;;; an adaptive calculation (working?)
#+(or)
(setq *result*
  (time
   (let ((dim 2) (order 2))
     (solve
      (make-instance
       '<stationary-fe-strategy>
       :fe-class (lagrange-fe order)
       :solver (s1-reduction-amg-solver order) ; #-(or)(lu-solver)
       :estimator (make-instance '<projection-error-estimator>)
       :indicator (make-instance '<largest-eta-indicator> :fraction 0.5)
       :success-if '(and (>= :nr-levels 2)(< :global-eta 1.0d-3))
       :output t)
      (blackboard :problem (inlay-cell-problem dim 0.1))))))

(effective-tensor *result*)
(getbb *result* :global-eta)
(fe-extreme-values (getbb *result* :solution))

;; chequerboard cell
(let ((dim 2) (order 1))
  (effective-tensor
   (cell-solve (chequerboard-problem dim 0.1):level 3 :order order)))

;; examine properties
(describe (simple-ball-inlay-cell-problem 1 0.1))

)

;;; (cdr-hom-tests)
(fl.tests:adjoin-test 'cdr-hom-tests)

