;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;  cdr.lisp
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

(in-package :cl-user)
(defpackage "FL.CDR"
  (:use "COMMON-LISP" "FL.MACROS" "FL.UTILITIES" "FL.MATLISP" "FL.FUNCTION"
	"FL.MESH" "FL.PROBLEM")
  (:export "<CDR-PROBLEM>" "MAP-DOMAIN-TO-CDR-PROBLEM"
	   "SCALAR-DIFFUSION" "IDENTITY-DIFFUSION-TENSOR"
	   "CDR-MODEL-PROBLEM" "CDR-NONLINEAR-RHS-PROBLEM"
	   "BRATU-PROBLEM")
  (:documentation "Defines convection-diffusion-reaction problems"))
(in-package "FL.CDR")

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; <cdr-problem>
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass <cdr-problem> (<pde-problem>)
  ()
  (:documentation "Convection-diffusion-reaction problem."))

(defmethod initialize-instance :after
    ((problem <cdr-problem>) &key &allow-other-keys)
  "The problem is assumed to be coercive if there is a reaction
coefficient."
  (unless (property-set-p problem 'coercive)
    (dohash ((coeffs) (coefficients problem))
      (when (member 'REACTION coeffs)
	(setf (get-property problem 'coercive) T)))))

(defmethod nr-of-components ((problem <cdr-problem>)) 1)

(defmethod interior-coefficients ((problem <cdr-problem>))
  "Interior coefficients for a convection-diffusion-reaction problem."
  (list* 'DIFFUSION 'CONVECTION 'SOURCE 'REACTION 'GAMMA 'FE-RHS
	 (call-next-method)))

(defmethod self-adjoint-p ((problem <cdr-problem>))
  (doskel (patch (domain problem))
    (when (get-coefficient (coefficients-of-patch patch problem)
		'FL.CDR::CONVECTION)
      (return-from self-adjoint-p (values NIL T))))
  (values T T))

(defmethod dual-problem ((problem <cdr-problem>) cell->rhs)
  "Dual problem of a cdr problem.  At the moment it works only for
selfadjoint problems with the functional of interest being the load
functional."
  (make-instance
   '<cdr-problem> :domain (domain problem)
   :patch->coefficients
   #'(lambda (cell)
       ;; check that problem is self adjoint
       (let ((coeffs (copy-seq (coefficients-of-patch cell problem))))
	 (when (get-coefficient coeffs 'CONSTRAINT)
	   (setf (getf coeffs 'CONSTRAINT)
		 (constant-coefficient 0.0)))
	 (assert (not (get-coefficient coeffs 'FL.CDR::CONVECTION)))  ; better: change to negative
	 (unless (eq cell->rhs :load-functional)
	   (let ((dual-rhs (funcall cell->rhs cell)))
	     (setf (getf coeffs 'FL.CDR::SOURCE)
		   (getf dual-rhs 'FL.CDR::SOURCE))
	     (setf (getf coeffs 'FL.CDR::GAMMA)
		   (getf dual-rhs 'FL.CDR::GAMMA))
	     (setf (getf coeffs 'FL.CDR::FE-RHS)
		   (getf dual-rhs 'FL.CDR::FE-RHS))))
	 coeffs))
   :multiplicity (multiplicity problem)))

(defmethod zero-constraints ((problem <cdr-problem>))
  (constant-coefficient 0.0))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Generation of standard cdr problems
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun scalar-diffusion (dim value)
  (constant-coefficient (scal value (eye dim))))

(defun identity-diffusion-tensor (dim)
  (scalar-diffusion dim 1.0))

(defun cdr-model-problem (domain &key (diffusion nil diffusion-p)
			  (source nil source-p) (dirichlet nil dirichlet-p)
			  gamma convection reaction initial evp properties)
  "Generates a convection-diffusion-reaction model problem.  Defaults are
identity diffusion, right-hand-side equal 1, and Dirichlet zero boundary
conditions.  Ordinary function are converted into coefficient functions
depending on a global coordinate.  The first argument can be either a
domain or an integer n which is interpreted as the n-dimensional unit
cube."
  (setq domain (ensure-domain domain))
  (let ((dim (dimension domain)))
    ;; set default values
    (unless diffusion-p (setq diffusion (eye dim)))
    (unless source-p (setq source (if evp 0.0 1.0)))
    (unless dirichlet-p (setq dirichlet 0.0))
    (apply #'fl.amop:make-programmatic-instance
	   (cond (initial '(<time-dependent-problem> <cdr-problem>))
		 (evp '(<cdr-problem> <evp-mixin>))
		 (t '<cdr-problem>))
	   :properties properties
	   :domain domain :patch->coefficients
	   `((:external-boundary
	      ,(when dirichlet
		     `(CONSTRAINT ,(ensure-coefficient dirichlet))))
	     (:d-dimensional
	      ,(append
		(when diffusion
		  `(FL.CDR::DIFFUSION ,(ensure-coefficient diffusion)))
		(when source
		  `(FL.CDR::SOURCE ,(ensure-coefficient source)))
		(when convection
		  `(FL.CDR::CONVECTION ,(ensure-coefficient convection)))
		(when reaction
		  `(FL.CDR::REACTION ,(ensure-coefficient reaction)))
		(when gamma
		  `(FL.CDR::GAMMA ,(ensure-coefficient gamma)))
		(when initial
		  `(FL.CDR::INITIAL ,(ensure-coefficient initial))))))
	   (append
	    (when evp (destructuring-bind (&key lambda mu) evp
			(list :lambda lambda :mu mu))))
	   )))

;;; nonlinear cdr problem

(defun cdr-nonlinear-rhs-problem (domain f &rest args &key source reaction &allow-other-keys)
  "Returns the Newton linearization @math{-\Delta u + F'(u) u = F'(u) u -
F(u)} for the nonlinear problem @math{-\Delta u +F(u) =0}."
  (assert (not (or source reaction)))
  (setq domain (ensure-domain domain))
  (apply #'cdr-model-problem
	 domain
	 :reaction
	 (make-instance '<coefficient> :demands '((:fe-parameters :solution)) :residual nil
			:eval #'(lambda (&key solution &allow-other-keys)
				  (evaluate-gradient f solution)))
	 :source
	 (make-instance '<coefficient> :demands '((:fe-parameters :solution))
			:eval #'(lambda (&key solution &allow-other-keys)
				  (gemm 1.0 (evaluate-gradient f solution) solution
					-1.0 (evaluate f solution))))
	 args))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Bratu problem
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun bratu-problem (dim)
  "Returns a linearization for the Bratu problem @math{-Delta u +e^u =0}."
  (cdr-nonlinear-rhs-problem
   dim (special-1d-function  #'exp #'exp)))

;;; Testing: (test-cdr)

(defun test-cdr ()
  (assert (get-property (cdr-model-problem 1) 'linear-p))
  (coefficients (cdr-model-problem 1))
  (assert (not (get-property (bratu-problem 2) 'linear-p)))
  (check (domain (cdr-model-problem 2)))
  (let* ((domain (n-cell-domain 1))
	 (problem (cdr-model-problem domain)))
    (describe problem))
  (describe (cdr-model-problem 2 :initial (constantly 1.0)))
  )

(fl.tests:adjoin-test 'test-cdr)

