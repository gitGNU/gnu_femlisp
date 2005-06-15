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
  (:use "COMMON-LISP" "FL.MACROS" "FL.UTILITIES" "FL.MATLISP"
	"FL.MESH" "FL.PROBLEM")
  (:export "<CDR-PROBLEM>" "MAP-DOMAIN-TO-CDR-PROBLEM"
	   "SCALAR-DIFFUSION" "IDENTITY-DIFFUSION-TENSOR"
	   "CDR-MODEL-PROBLEM" "BRATU-PROBLEM")
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
    (when (getf (coefficients-of-patch patch problem)
		'FL.CDR::CONVECTION)
      (return-from self-adjoint-p (values NIL T))))
  (values T T))

(defmethod dual-problem ((problem <cdr-problem>) cell->rhs)
  "Dual problem of a cdr problem.  To be improved, at the moment only for
selfadjoint problems with functional of interest being the load
functional."
  (make-instance
   '<cdr-problem> :domain (domain problem)
   :patch->coefficients
   #'(lambda (cell)
       ;; check that problem is self adjoint
       (let ((coeffs (copy-seq (coefficients-of-patch cell problem))))
	 (when (getf coeffs 'FL.CDR::CONSTRAINT)
	   (setf (getf coeffs 'FL.CDR::CONSTRAINT) (constant-coefficient 0.0)))
	 (assert (not (getf coeffs 'FL.CDR::CONVECTION)))  ; better: change to negative
	 (unless (eq cell->rhs :load-functional)
	   (let ((dual-rhs (funcall cell->rhs cell)))
	     (setf (getf coeffs 'FL.CDR::SOURCE) (getf dual-rhs 'FL.CDR::SOURCE))
	     (setf (getf coeffs 'FL.CDR::GAMMA) (getf dual-rhs 'FL.CDR::GAMMA))
	     (setf (getf coeffs 'FL.CDR::FE-RHS) (getf dual-rhs 'FL.CDR::FE-RHS))))
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

(defun cdr-model-problem (dim/domain &key (diffusion nil diffusion-p)
			  (source nil source-p) (dirichlet nil dirichlet-p)
			  gamma convection reaction initial evp properties)
  "Generates a convection-diffusion-reaction model problem.  Defaults are
identity diffusion, right-hand-side equal 1, and Dirichlet zero boundary
conditions.  Ordinary function are converted into coefficient functions
depending on a global coordinate.  The first argument can be either a
domain or an integer n which is interpreted as the n-dimensional unit
cube."
  (let* ((domain (if (numberp dim/domain)
		     (n-cube-domain dim/domain)
		     dim/domain))
	 (dim (dimension domain))
	 (bdry (skeleton-boundary domain)))
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
	   #'(lambda (cell)
	     (let ((coeffs ()))
	       (when (member-of-skeleton? cell bdry)
		 (when dirichlet
		   (setf (getf coeffs 'FL.CDR::CONSTRAINT) (ensure-coefficient dirichlet))))
	       (when (= (dimension cell) dim)
		 (when diffusion
		   (setf (getf coeffs 'FL.CDR::DIFFUSION) (ensure-coefficient diffusion)))
		 (when source
		   (setf (getf coeffs 'FL.CDR::SOURCE) (ensure-coefficient source)))
		 (when convection
		   (setf (getf coeffs 'FL.CDR::CONVECTION) (ensure-coefficient convection)))
		 (when reaction
		   (setf (getf coeffs 'FL.CDR::REACTION) (ensure-coefficient reaction)))
		 (when gamma
		   (setf (getf coeffs 'FL.CDR::GAMMA) (ensure-coefficient gamma)))
		 (when initial
		   (setf (getf coeffs 'FL.CDR::INITIAL) (ensure-coefficient initial))))
	       coeffs))
	   (append
	    (when evp (destructuring-bind (&key lambda mu) evp
			(list :lambda lambda :mu mu))))
	   )))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Bratu problem
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun bratu-problem (dim)
  "Returns the Newton linearization (-Delta + e^u) u = e^u (u-1) for the
Bratu problem -Delta u +e^u =0."
  (cdr-model-problem
   (n-cube-domain dim)
   :diffusion (identity-diffusion-tensor dim)
   :reaction
   (make-instance '<coefficient> :demands '(:solution) :residual nil
		  :eval #'(lambda (&key solution &allow-other-keys)
			    (exp (vref solution 0))))
   :source
   (make-instance '<coefficient> :demands '(:solution)
		  :eval #'(lambda (&key solution &allow-other-keys)
			    (let ((u (vref solution 0)))
			      (* (exp u) (- u 1.0)))))))

;;; Testing: (test-cdr)

(defun test-cdr ()
  (assert (get-property (cdr-model-problem 1) 'linear-p))
  (assert (not (get-property (bratu-problem 2) 'linear-p)))
  (check (domain (cdr-model-problem 2)))
  (let* ((domain (n-cell-domain 1))
	 (problem (cdr-model-problem domain)))
    (describe problem))
  (describe (cdr-model-problem 2 :initial (constantly 1.0)))
  )

(fl.tests:adjoin-test 'test-cdr)

