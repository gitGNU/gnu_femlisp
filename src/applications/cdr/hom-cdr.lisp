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

(in-package application)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Utilities
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun cell-solve (problem &key (level 1) (order 1) parametric (solver *lu-solver*))
  "Solves the given problem and returns a result in form of an
assembly-line."
  (let* ((domain (domain problem))
	 (mesh (uniformly-refined-hierarchical-mesh domain level :parametric parametric))
	 (fe-class (lagrange-fe order))
	 (as (make-fe-ansatz-space fe-class problem mesh))
	 (assembly-line (fe-discretize (make-assembly-line :ansatz-space as))))
    (setf (get-al assembly-line :problem) problem)
    (setf (get-al assembly-line :mesh) mesh)
    (with-items (&key solution matrix rhs)
	assembly-line
      (setq solution (solve solver :matrix matrix :rhs rhs))
      assembly-line)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Cell problems
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; cell problems defined by the coefficient function on regular meshes

(defun cell-problem (dim &key diffusion-function)
  "Returns cell problem on the n-cell-domain of the given dimension.  Gamma
yields dim right hand sides of the form A e_k, i.e. the columns of the
diffusion tensor."
  (make-instance
   '<cdr-problem> :domain (n-cell-domain dim)
   :multiplicity dim :patch->coefficients
   #'(lambda (patch)
       (when (= (dimension patch) dim)
	 (list 'CDR::DIFFUSION (function->coefficient diffusion-function)
	       'CDR::GAMMA (constant-coefficient (eye dim)))))))

(defun simple-square-inlay-cell-problem (dim)
  (cell-problem
   dim :diffusion-function
   #'(lambda (x) 
       (let ((result (eye dim)))
	 (when (every #'(lambda (coord) (<= 0.25 coord 0.75)) x)
	    (scal! 0.1 result))
	  result))))

(defun smooth-coefficient-cell-problem (dim)
  (cell-problem
   dim :diffusion-function
   #'(lambda (x)
       (scal (reduce #'* (map 'vector #'(lambda (xc) #I(2.0+sin(2*pi*xc))) x))
	     (eye dim)))))
  
(defun simple-ball-inlay-cell-problem (dim eps)
  (cell-problem
   dim :diffusion-function
   #'(lambda (x)
       (let ((result (eye dim)))
	 (when (<= (norm (m- (make-double-vec dim 0.5d0) x)) 0.25)
	   (scal! eps result))
	 result))))

(defun chequerboard-problem (dim eps)
  "Returns cell problem on the n-cell-domain of the given dimension with
a chequerboard pattern."
  (cell-problem
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
   '<cdr-problem> :domain (n-cell-with-n-ball-inlay dim)
   :multiplicity dim :patch->coefficients
   #'(lambda (patch)
       (when (= (dimension patch) dim)
	 (list 'CDR::DIFFUSION
	       (constant-coefficient
		(scal! (if (patch-in-inlay-p patch) eps 1.0) (eye dim)))
	  'CDR::GAMMA (constant-coefficient (eye dim)))))))

(defun porous-cell-problem (dim)
  (make-instance
   '<cdr-problem> :domain (n-cell-with-n-ball-hole dim :radius 0.3)
   :multiplicity dim :patch->coefficients
   #'(lambda (patch)
       (when (= (dimension patch) dim)
	 (list 'CDR::DIFFUSION (constant-coefficient (eye dim))
	       'CDR::GAMMA (constant-coefficient (eye dim)))))))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Demos
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun cdr-interior-effective-coeff-demo (problem order levels &key output plot)
  "Computes an effective diffusion tensor for different configurations.
The approximation is done with isoparametric finite elements.  Uniform
refinement is used, the linear solver is a geometric multigrid scheme used
in a nested iteration fashion.  The solution to the cell problem is a
vector field whose components are plotted one after the other.  The setting
has cubic symmetry, from which one can easily see that the effective tensor
must be a scalar multiple of the identity."
  (setf (getf *strategy-output* :observe)
	`((:Ahom "               Ahom"
	   ,#'(lambda (assembly-line)
		(with-items (&key solution rhs) assembly-line
		  (format t "~19,10,2E"
			  (and solution rhs
			       (matrix-ref (apply #'effective-tensor assembly-line) 0 0))))))))
  (defparameter *result*
    (solve-with
     (make-instance
      '<fe-strategy> :fe-class (lagrange-fe order)
      :estimator (make-instance '<projection-error-estimator>)
      :indicator (make-instance '<largest-eta-indicator> :fraction 1.0)
      :appraise (stop-if :nr-levels>= levels)
      :solver (make-instance
	       '<linear-solver> :iteration
	       (let ((smoother
		      (if (>= (dimension (domain problem)) 3)
			  *gauss-seidel*
			  (make-instance '<local-bgs> :type :vertex-centered))))
		 (geometric-cs
		  :gamma 2 :fmg nil :coarse-grid-iteration
		  (make-instance '<multi-iteration> :base smoother
				 :nr-steps (if (eq smoother *gauss-seidel*) 10 3))
		  :pre-steps 2 :pre-smooth smoother :post-steps 2 :post-smooth smoother))
	       :success-if `(:defnorm< 1.0e-10) ; (:reduction< ,(* 0.1 (expt 0.5 (* 2 (+ order 1)))))
	       :failure-if `(:step-reduction> 0.9)
	       :output (eq output :all))
      :output t)
     problem
     ))
  ;; plot cell solutions and compute the homogenized coefficient
  (when plot
    (let ((solution (getf *result* :solution)))
      (plot solution :index 0)
      (sleep 1.0)
      (plot solution :index 1)))
  (format t "The effective tensor is:~%~A~%"
	  (apply #'effective-tensor *result*)))


;;; Testing:
#+(or)(cdr-interior-effective-coeff-demo (porous-cell-problem 2) 3 3)
#+(or)(cdr-interior-effective-coeff-demo (inlay-cell-problem 2 0.1) 5 1)
;;; Gauss-Lobatto points
;;; inexact domain
;     9       450     31597       4.3   7.2318763792d-01
;    36      2250    159768      35.7   7.2310081411d-01
; ;;; exact domain representation:
;     9       450     31597       5.3   7.2315790447d-01
;    36      2250    159768      46.9   7.2310037639d-01
; ;; standard Lagrange points
; ;; inexact domain
;     9       450     31597       8.6   7.2320277654d-01
;    36      2250    159768      50.0   7.2310114601d-01
; ;; exact domain -> Lobatto sind nur fuer Randabb wichtig
;     9       450     31597       5.7   7.2315790447d-01
;    36      2250    159768      48.0   7.2310037639d-01

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

(make-effective-diffusion-porous-domain-demo 2 5 3)
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

(defparameter *result* nil)

(defun cdr-hom-tests ()
;; Warning: The default *lu-solver* in cell-solve does not really work on
;; these singular systems.  Often it appears to work, but introduces large
;; relative errors.

(cdr-interior-effective-coeff-demo (porous-cell-problem 2) 4 2)

(cdr-interior-effective-coeff-demo (inlay-cell-problem 2 0.1) 1 3 :output :all)

;;grid is not inlay-adapted
(plot (simple-ball-inlay-cell-problem 1 0.1) :refinements 3
	:coefficient 'CDR::DIFFUSION :key (rcurry #'mat-ref 0 0))

;; solve cell problem and compute homogenized coefficient
(defparameter *result*
  (cell-solve (simple-ball-inlay-cell-problem 2 0.1)
	      :level 3 :order 1))
(apply #'effective-tensor *result*)
(plot (getf *result* :solution) :index 0)

;; generate 3D cell problem
(simple-ball-inlay-cell-problem 3 0.1)

;; inlay adapted grid
(plot (inlay-cell-problem 2 0.1d0) :refinements 0 :depth 2 :parametric (lagrange-mapping 3)
	:coefficient 'CDR::DIFFUSION :key (rcurry #'mat-ref 0 0))

;; first order with lu-solver
(setq *result*
      (cell-solve (inlay-cell-problem 2 0.1) :level 1
		  :order 1 :parametric (lagrange-mapping 2)))
(apply #'effective-tensor *result*)
(plot (getf *result* :solution) :index 0)

;; first order with full multigrid
(defparameter *result*
  (time
   (let ((level 2) (order 1))
     (cell-solve
      (inlay-cell-problem 2 0.1) :level level
      :order order :parametric (lagrange-mapping order)
      :solver				; *lu-solver* #+(or)
      (make-instance
       '<linear-solver> :iteration
       (geometric-cs
	:fmg t :coarse-grid-iteration
	(make-instance '<multi-iteration> :base *gauss-seidel* :nr-steps 3))
       :output t
       :reduction 1.0d-5
       )))))
(plot (getf *result* :solution) :index 0)
(plot (getf *result* :solution) :index 1)

(defparameter *result*
  (cell-solve (inlay-cell-problem 2 0.1) :level 2
		:order 1 :parametric (lagrange-mapping 2)))
(plot (getf *result* :mesh)) ; :plot :file :format "tiff")

;; higher order
(time
 (defparameter *result*
   (let ((dim 2) (level 1) (order 4))
     (cell-solve
      (inlay-cell-problem dim 0.1) :level level
      :order order :parametric (lagrange-mapping (max order 2))
      :solver ;*lu-solver* #+(or)
      (make-instance
       '<linear-solver> :iteration
       (let ((smoother *gauss-seidel*))
	 (geometric-cs
	  :gamma 2 :fmg t :coarse-grid-iteration
	  (make-instance '<multi-iteration> :base *gauss-seidel* :nr-steps 10)
	  :pre-steps 2 :pre-smooth smoother :post-steps 2 :post-smooth smoother))
       :reduction (expt 0.5 (* level (+ order 2))))))))

;; (l=2, o=4, cgs=10 ILU, 3 smooth, fmg)*2=146.7 URT
(apply #'effective-tensor *result*)
(plot (getf *result* :solution) :index 1)

;;; an adaptive calculation (working?)
(defparameter *result*
  (time
   (let ((dim 2) (order 2))
     (solve-with
      (make-instance
       '<fe-strategy> :fe-class (lagrange-fe order)
       :estimator (make-instance '<projection-error-estimator>)
       :indicator (make-instance '<largest-eta-indicator> :fraction 0.5)
       :appraise #-(or) (global-estimate-smaller-than 1.0d-3) #+(or) (level= 4)
       :solver #-(or)(s1-reduction-amg-solver order)
       #+(or)*lu-solver*
       :output t)
      (inlay-cell-problem dim 0.1)))))

(apply #'effective-tensor *result*)
(getf *result* :global-eta)

;; chequerboard cell
(let ((dim 2) (order 1))
  (apply #'effective-tensor
	 (cell-solve (chequerboard-problem dim 0.1d0):level 3
		     :order order)))

;; examine properties
(let* ((problem (simple-ball-inlay-cell-problem 1 0.1))
       (h-mesh (uniformly-refined-hierarchical-mesh (domain problem) 1)))
  (problem-info problem h-mesh))

)

;;; (cdr-hom-tests)
(tests::adjoin-femlisp-test 'cdr-hom-tests)
