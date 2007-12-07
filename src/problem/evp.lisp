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
  ((lambda :initarg :lambda :initform (box 0.0) :documentation
       "The multiplier for the mass matrix, usually equal to the
eigenvalue.")
   (mu :initform (box 1.0) :initarg :mu :documentation
       "The multiplier for the system matrix."))
  (:documentation "A mixin for eigenvalue problems."))

(defmethod initialize-instance :after ((problem <evp-mixin>)
				       &key &allow-other-keys)
  (setf (getf (properties problem) 'linear-p) nil))

(defclass <evp> (<evp-mixin> <nonlinear-problem>)
  ()
  (:documentation "Standard class for discrete eigenvalue problems."))

(defclass <ls-evp> (<evp>)
  ((A :initarg :A :documentation "(Energy) matrix A.")
   (B :initarg :B :documentation "(Mass) matrix B."))
  (:documentation "Generalized eigenvalue problem for matrices."))

(defmethod initialize-instance :after ((lsevp <ls-evp>) &key &allow-other-keys)
  (with-slots (linearization initial-guess lambda mu A B) lsevp
    (setf linearization
	  #'(lambda (solution)
	      (declare (ignore solution))
	      (lse :matrix (m- (scal (unbox mu) A) (scal (unbox lambda) B))
		   :rhs (make-image-vector-for A))))
    (ensure initial-guess (mrandom (nrows A) 1))))

(defgeneric mass (evp x)
  (:documentation
   "Evaluates the mass bilinear form for a generalized eigenvalue problem.")
  (:method ((evp <evp>) x)
	   (with-slots (mu lambda)
	     evp
	     (fluid-let (((unbox mu) 0.0)
			 ((unbox lambda) -1.0))
	       (let ((lse (linearize evp x)))
		 (dot x (m* (matrix lse) x)))))))

(defgeneric energy (evp x)
  (:documentation
   "Evaluates the energy bilinear form for a generalized eigenvalue problem.")
  (:method ((evp <evp>) x)
	   (with-slots (mu lambda)
	     evp
	     (fluid-let (((unbox mu) 1.0)
			 ((unbox lambda) 0.0))
	       (let ((lse (linearize evp x)))
		 (dot x (m* (matrix lse) x)))))))

(defun test-evp ()
  ;; finding an eigenvalue/eigenfunction by a Wielandt iteration
  (let* ((n 10)
	 (evp (make-instance '<ls-evp> :lambda (box 0.0)
			     :A (scal (expt (+ n 1.0) 2) (laplace-full-matrix n))
			     :B (eye n)))
	 (bb (blackboard :solution (mrandom n 1))))
    (with-items (&key linearization solution residual residual-p) bb
      (loop repeat 20
	 initially (ensure-residual evp bb)
	 until (< (norm residual) 1.0e-10) do
	    (let ((new-sol
		   (ignore-errors
		     (gesv (matrix linearization) solution))))
	      (if new-sol
		  (setf solution new-sol)
		  (return)))
	    (let ((mass (mass evp solution))
		  (energy (energy evp solution)))
	      (setf (unbox (slot-value evp 'lambda)) (/ energy mass))
	      (scal! (/ (sqrt mass)) solution))
	    (setf residual-p nil)
	    (ensure-residual evp bb)
	    (format t "lambda   = ~A~%solution = ~A~%resnorm  = ~A~%~%"
		    (unbox (slot-value evp 'lambda)) solution (norm residual)))
      (let ((lam (unbox (slot-value evp 'lambda))))
	(assert (< (- lam (expt pi 2)) 0.01))
	(values lam solution))))
  )

;;; (test-evp)
(fl.tests:adjoin-test 'test-evp)
