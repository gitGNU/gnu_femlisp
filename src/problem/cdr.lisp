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
  (:use "COMMON-LISP" "FL.MACROS" "FL.UTILITIES"
	"FL.MATLISP" "FL.ALGEBRA" "FL.FUNCTION"
	"FL.MESH" "FL.PROBLEM" "FL.ELLSYS")
  (:export "<CDR-PROBLEM>"
	   "SCALAR-DIFFUSION" "IDENTITY-DIFFUSION-TENSOR"
	   "SCALAR-SOURCE" "SCALAR-REACTION"
	   "ENSURE-TENSOR-COEFFICIENT" "ENSURE-VECTOR-COEFFICIENT"
	   "ENSURE-DIRICHLET-COEFFICIENT"
	   "CDR-MODEL-PROBLEM" "CDR-NONLINEAR-RHS-PROBLEM"
	   "BRATU-PROBLEM")
  (:documentation "Defines convection-diffusion-reaction problems"))

(in-package "FL.CDR")

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; <cdr-problem>
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass <cdr-problem> (<ellsys-problem>)
  ((components :initform '(u) :initarg :components))
  (:documentation "Convection-diffusion-reaction problem."))

(defmethod self-adjoint-p ((problem <cdr-problem>))
  "The problem is assumed to be self-adjoint if there is no convection
coefficient."
  (doskel (patch (domain problem))
    (when (get-coefficient (coefficients-of-patch patch problem)
			   'FL.ELLSYS::B)
      (return-from self-adjoint-p (values NIL T))))
  (values T T))

(defmethod dual-problem ((problem <cdr-problem>) cell->rhs)
  "Dual problem of a cdr problem.  At the moment it works only for
selfadjoint problems with the functional of interest being the load
functional."
  (let ((mult (multiplicity problem)))
    (make-instance
     '<cdr-problem> :components '(u) :multiplicity mult
     :domain (domain problem)
     :coefficients
     #'(lambda (cell)
	 (lret ((coeffs (coefficients-of-patch cell problem)))
	   ;; check that problem is self adjoint
	   (when (find 'FL.ELLSYS::B coeffs :key #'coefficient-name)
	     (error "NYI")) ; better: change to negative
	   (unless (eq cell->rhs :load-functional)
	     (setf coeffs
		   (append (and cell->rhs (funcall cell->rhs cell))
			   (remove-if (lambda (coeff)
					(member (coefficient-name coeff)
						'(FL.ELLSYS::F 'FL.ELLSYS::G 'FL.ELLSYS::FE-RHS)))
				      coeffs)))))))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Generation of standard cdr problems
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun scalar-diffusion (dim value)
  (constant-coefficient
   'FL.ELLSYS::A
   (diagonal-sparse-tensor (vector (scal value (eye dim))))))

(defun scalar-source (value)
  (constant-coefficient
   'FL.ELLSYS::F
   (vector (if (standard-matrix-p value)
	       value
	       (scal value (ones 1))))))

(defun scalar-reaction (value)
  (constant-coefficient
   'FL.ELLSYS::R
   (diagonal-sparse-tensor (vector (scal value (ones 1))))))

(defun ensure-1-component-tensor (value)
  (when (numberp value)
    (setq value (ensure-matlisp value)))
  (when (standard-matrix-p value)
    (setq value (diagonal-sparse-tensor value 1)))
  value)

(defun ensure-tensor-coefficient (name obj)
  (cond ((typep obj '<coefficient>) obj)
	((or (functionp obj) (typep obj '<function>))
	 (make-instance '<coefficient> :name name :demands '(:global)
			:eval #'(lambda (&key global &allow-other-keys)
				  (ensure-1-component-tensor (evaluate obj global)))))
	(t (constant-coefficient name (ensure-1-component-tensor obj)))))

(defun ensure-1-component-vector (value)
  (when (numberp value)
    (setq value (ensure-matlisp value)))
  (when (standard-matrix-p value)
    (setq value (vector value)))
  value)

(defun ensure-vector-coefficient (name obj)
  (cond ((typep obj '<coefficient>) obj)
	((or (functionp obj) (typep obj '<function>))
	 (make-instance '<coefficient> :name name :demands '(:global)
			:eval #'(lambda (&key global &allow-other-keys)
				  (ensure-1-component-vector (evaluate obj global)))))
	(t (constant-coefficient name (ensure-1-component-vector obj)))))

(defun ensure-dirichlet-coefficient (obj)
  (cond ((typep obj '<coefficient>) obj)
	((or (functionp obj) (typep obj '<function>))
	 (make-instance '<coefficient> :name 'CONSTRAINT :demands '(:global)
			:eval #'(lambda (&key global &allow-other-keys)
				  (let ((result (evaluate obj global)))
				    (values (make-array 1 :initial-element t)
					    (ensure-1-component-vector result))))))
	(t (constant-coefficient
	    'CONSTRAINT (make-array 1 :initial-element t)
	    (ensure-1-component-vector obj)))))

(defun cdr-model-problem (domain &key (diffusion nil diffusion-p)
			  (source nil source-p) (dirichlet nil dirichlet-p)
			  gamma convection reaction sigma initial
			  evp (multiplicity 1)
			  properties)
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
    ;; (setq diffusion (ensure-1-component-tensor diffusion))
    ;;(setq convection (ensure-1-component-tensor convection))
    ;; (setq source (ensure-1-component-vector source))
    ;; (setq gamma (ensure-1-component-vector gamma))
    ;;(setq dirichlet (ensure-dirichlet-coefficient dirichlet))
    (when initial (ensure sigma 1.0))
    (fl.amop:make-programmatic-instance
     (cond ((and initial sigma) '(<time-dependent-problem> <cdr-problem>))
	   (evp '(<cdr-problem> <evp-mixin>))
	   (t '<cdr-problem>))
     :properties properties
     :domain domain :components '(u) :patch->coefficients
     `((:external-boundary
	,(when dirichlet `(,(ensure-dirichlet-coefficient dirichlet))))
       (:d-dimensional
	,(append
	  (when diffusion
	    `(,(ensure-tensor-coefficient 'FL.ELLSYS::A diffusion)))
	  (when source
	    `(,(ensure-vector-coefficient 'FL.ELLSYS::F source)))
	  (when convection
	    `(,(ensure-tensor-coefficient 'FL.ELLSYS::B convection)))
	  (when reaction
	    `(,(ensure-tensor-coefficient 'FL.ELLSYS::R reaction)))
	  (when gamma
	    `(,(ensure-vector-coefficient 'FL.ELLSYS::G gamma)))
	  (when initial
	    `(,(ensure-vector-coefficient 'INITIAL initial)))
	  (when sigma
	    `(,(ensure-tensor-coefficient 'FL.ELLSYS::SIGMA sigma)))
	  )))
     :multiplicity multiplicity)
    ))

;;; nonlinear cdr problem

(defun cdr-nonlinear-rhs-problem (domain f &rest args &key source reaction &allow-other-keys)
  "Returns the Newton linearization @math{-\Delta u + F'(u) u = F'(u) u -
F(u)} for the nonlinear problem @math{-\Delta u +F(u) =0}."
  (assert (not (or source reaction)))
  (setq domain (ensure-domain domain))
  (destructuring-bind (source reaction)
      (linearization f)
    (apply #'cdr-model-problem
	   domain
	   :reaction
	   (make-instance '<coefficient> :name 'FL.ELLSYS::R
			  :demands '((:fe-parameters :solution)) :residual nil
			  :eval #'(lambda (&key solution &allow-other-keys)
				    (let ((solution (first solution)))
				      (funcall reaction solution))))
	   :source
	   (make-instance '<coefficient> :name 'FL.ELLSYS::F
			  :demands '((:fe-parameters :solution))
			  :eval #'(lambda (&key solution &allow-other-keys)
				    (let ((solution (first solution)))
				      (funcall source solution))))
	   args)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; new problem definition support
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmethod make-coefficients-for ((problem fl.cdr::<cdr-problem>)
				  (coeff (eql 'FL.CDR::DIFFUSION))
				  patch args eval)
  (make-coefficients-for problem 'FL.ELLSYS::A patch args
			 (lambda (&rest args)
			   (diagonal-sparse-tensor (apply eval args) 1))))
		    
(defmethod make-coefficients-for ((problem fl.cdr::<cdr-problem>)
				 (coeff (eql 'FL.CDR::ISOTROPIC-DIFFUSION))
				  patch args eval)
  (make-coefficients-for problem 'FL.ELLSYS::A patch args
			 (let ((dim (dimension (domain problem))))
			   (lambda (&rest args)
			     (diagonal-sparse-tensor (scal (apply eval args) (eye dim)) 1)))))

(defmethod make-coefficients-for ((problem fl.cdr::<cdr-problem>)
				 (coeff (eql 'FL.CDR::SCALAR-SOURCE)) patch args eval)
  (make-coefficients-for problem 'FL.ELLSYS::F patch args
			 (lambda (&rest args)
			   (vector (ensure-matlisp (apply eval args) :row-vector)))))
		    
(defmethod make-coefficients-for ((problem fl.cdr::<cdr-problem>)
				 (coeff (eql 'FL.CDR::SCALAR-NONLINEAR-SOURCE))
				  patch args eval)
  (make-coefficients-for problem 'FL.ELLSYS::NONLINEAR-F patch args
			 (lambda (&rest args)
			   (sparse-tensor 
			    `((,(apply eval args) 0))))))

(defmethod make-coefficients-for ((problem fl.cdr::<cdr-problem>)
				 (coeff (eql 'FL.CDR::SCALAR-CONSTRAINT))
				  patch args eval)
  (make-coefficients-for problem 'FL.PROBLEM::CONSTRAINT patch args
			 (lambda (&rest args)
			   (values #(t) (vector (ensure-matlisp (apply eval args) :row-vector))))))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Bratu problem
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun bratu-reaction-jet ()
  "Returns a function which computes the jet (value and derivative) of the
Bratu reaction term."
  (lambda (u)
    (let ((u (vref (vref u 0) 0)))
      (list
       (vector (ensure-matlisp (exp u)))
       (lret ((result (make-instance '<sparse-tensor> :rank 2)))
	 (setf (tensor-ref result 0 0) (ensure-matlisp (exp u))))))))

(defun bratu-problem (dim)
  "Returns a linearization for the Bratu problem @math{-Delta u +e^u =0}."
  (cdr-nonlinear-rhs-problem
   dim (bratu-reaction-jet)))

;;; Testing

(defun test-cdr ()

  #+(or)  ; fl.iteration is not necessarily available
  (let ((problem
	 (create-problem 'FL.CDR::<CDR-PROBLEM>
	     (:domain (n-cube-domain 2) :components '(u) :multiplicity 1)
	   (setup-coefficients (patch)
	     (select-on-patch (patch)
	       (:d-dimensional
		(list
		 (coeff FL.CDR::ISOTROPIC-DIFFUSION () 1.0)
		 (coeff FL.CDR::SCALAR-SOURCE () 1.0)))
	       (:external-boundary
		(coeff FL.CDR::SCALAR-CONSTRAINT () 0.0)))))))
    (solve (blackboard :problem problem :solver (fl.iteration::lu-solver)
                       :success-if '(> :step 3) :output :all)))

  (assert (get-property (cdr-model-problem 1) 'linear-p))
  (coefficients (cdr-model-problem 1))
  (assert (not (get-property (bratu-problem 2) 'linear-p)))
  (check (domain (cdr-model-problem 2)))
  (let* ((domain (n-cell-domain 1))
	 (problem (cdr-model-problem domain)))
    (describe problem))
  (describe (cdr-model-problem 2 :initial (constantly 1.0)))
  )

;;; (test-cdr)
(fl.tests:adjoin-test 'test-cdr)

