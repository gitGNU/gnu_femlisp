;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; evp.lisp - Eigenvalue problems
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;
;;; Copyright (C) 2004 Nicolas Neuss, University of Heidelberg.
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

(in-package :fl.problem)

(file-documentation
 "This file provides definitions for eigenvalue problems.")

(defclass <evp-mixin> ()
  ((multiplicity :reader multiplicity :initform 1 :initarg :multiplicity
		 :documentation "The multiplicity of the eigenspace.")
   (eigenvalues
    :initarg :eigenvalues :initform nil :documentation
    "The current approximation of the eigenvalues.")
   (mu :initform (box 1.0) :initarg :mu :documentation
       "The multiplier for the system matrix."))
  (:documentation "A mixin for eigenvalue problems."))

(defmethod initialize-instance :after ((problem <evp-mixin>)
				       &key &allow-other-keys)
  (with-slots (multiplicity eigenvalues) problem
    (ensure eigenvalues (make-double-vec multiplicity)))
  (setf (getf (properties problem) 'linear-p) nil)
  )

(defclass <evp> (<evp-mixin> <nonlinear-problem>)
  ()
  (:documentation "Standard class for discrete eigenvalue problems."))

(defclass <ls-evp> (<evp>)
  ((stiffness-matrix :reader stiffness-matrix :initarg :stiffness-matrix)
   (mass-matrix :reader mass-matrix :initarg :mass-matrix))
  (:documentation "Generalized eigenvalue problem for matrices."))

(defmethod initialize-instance :after ((lsevp <ls-evp>) &key &allow-other-keys)
  (with-slots (multiplicity linearization solution eigenvalues mu
			    stiffness-matrix mass-matrix) lsevp
    (setf linearization
	  (lambda (solution)
	    (let ((Ax (m* stiffness-matrix solution))
		  (Bx (m* mass-matrix solution)))
	      (loop for i from 0
		 and stiffness across (diagonal (m*-tn solution Ax))
		 and mass across (diagonal (m*-tn solution Bx)) do
		 (setf (aref eigenvalues i) (/ stiffness mass)))
	      ;; equation system for defect correction
	      (lse :matrix stiffness-matrix
		   :rhs (m* Bx (diag eigenvalues))))))
    (ensure solution
	    (fill-random! (make-domain-vector-for stiffness-matrix multiplicity)
			  1.0))))

(defgeneric mass (evp x)
  (:documentation
   "Evaluates the mass bilinear form for a generalized eigenvalue problem.")
  (:method ((evp <evp>) x)
	   (with-slots (mu eigenvalues)
	     evp
	     (fluid-let (((unbox mu) 0.0)
			 ((unbox eigenvalues) -1.0))
	       (let ((lse (linearize evp x)))
		 (dot x (m* (matrix lse) x)))))))

(defgeneric energy (evp x)
  (:documentation
   "Evaluates the energy bilinear form for a generalized eigenvalue problem.")
  (:method ((evp <evp>) x)
	   (with-slots (mu eigenvalues)
	     evp
	     (fluid-let (((unbox mu) 1.0)
			 ((unbox eigenvalues) 0.0))
	       (let ((lse (linearize evp x)))
		 (dot x (m* (matrix lse) x)))))))

(defun test-evp ()
  ;; finding an eigenvalue/eigenfunction by a Wielandt iteration
  (let* ((dim 1) (size 3) (n (expt size dim))
	 (evp (make-instance
	       '<ls-evp> :stiffness-matrix (laplace-full-matrix n)
			 :mass-matrix (eye n)))
	 (bb (blackboard :solution (mrandom n 1))))
    (with-items (&key linearization solution residual residual-p) bb
      (loop initially (ensure-residual evp bb)
	    until (mzerop  residual 1.0e-15) do
	   (aif (ignore-errors (gesv (matrix linearization) solution))
		(setf solution it)
		(return))
	   (let ((mass (mass evp solution))
		 (energy (energy evp solution)))
	     (setf (aref (slot-value evp 'eigenvalues) 0)
		   (/ energy mass))
	     (scal! (/ (sqrt mass)) solution))
	   (setf residual-p nil)
	   (ensure-residual evp bb)
	   (format t "eigenvalues   = ~A~%solution = ~A~%resnorm  = ~A~%~%"
		   (slot-value evp 'eigenvalues) solution (norm residual)))
      (values (slot-value evp 'eigenvalues) solution)))
  )

;;; (test-evp)
(fl.tests:adjoin-test 'test-evp)
