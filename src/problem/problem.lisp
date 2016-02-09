;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; problem.lisp
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

(in-package :fl.problem)

(defclass <problem> (property-mixin)
  ()
  (:documentation "Base class for all problems."))

(defgeneric linear-p (problem)
  (:documentation "Predicate determining if a problem is linear or nonlinear.")
  (:method ((problem <problem>))
    (get-property problem 'linear-p)))

(defgeneric ensure-residual (problem blackboard)
  (:documentation "Ensures that the field :RESIDUAL is computed and that
the flag :RESIDUAL-P is set on the blackboard."))

(defgeneric ensure-solution (problem blackboard)
  (:documentation "Ensures that the field :SOLUTION is set on the
blackboard.")
  (:method (problem blackboard)
    "Default method throws an error."
    (declare (ignore problem blackboard))
    (error "No initial guess for SOLUTION.")))

(defmethod ensure-solution :around (problem blackboard)
  "This :around method takes the field :SOLUTION from the blackboard."
  (declare (ignore problem))
  (or (getbb blackboard :solution)
      (call-next-method)))

(defmethod ensure-residual (problem blackboard)
  "This default method handles nonlinear problems by linearizing them and
computing the residual for it.  The resulting problem is additionally
stored in a field :LINEAR-PROBLEM."
  (when (linear-p problem)
    (error "No method ENSURE-RESIDUAL provided for this linear problem."))
  ;; linearize and compute residual
  (with-items (&key linearization solution) blackboard
    (setq linearization (linearize problem solution))
    (ensure-residual linearization blackboard)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Linear and nonlinear problems
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; Linear problems

(defclass <lse> (<problem>)
  ((matrix :reader matrix :initarg :matrix :initform nil)
   (rhs :reader rhs :initarg :rhs :initform nil))
  (:documentation "Standard form of a linear system of equations."))

(defun lse (&rest args)
  "Constructs a standard LSE."
  (apply #'make-instance '<lse> args))

(defmethod initialize-instance :after ((lse <lse>) &key &allow-other-keys)
  (setf (get-property lse 'linear-p) t))

(defgeneric linearize (problem solution)
  (:documentation "Linearize the nonlinear problem PROBLEM at the point
SOLUTION.  The result should be a linear problem.")
  (:method (problem solution)
    "This default method throws an error for nonlinear problems and is the
identity on linear problems."
    (declare (ignore solution))
    (if (get-property problem 'linear-p)
	problem
	(error "No method LINEARIZE provided for this problem."))))

(defmethod ensure-solution ((lse <lse>) blackboard)
  (setf (getbb blackboard :solution)
	(make-domain-vector-for
	 (matrix lse) (multiplicity (rhs lse)))))

(defmethod ensure-residual ((lse <lse>) blackboard)
  (with-items (&key solution residual residual-p) blackboard
    (ensure residual (make-image-vector-for
		      (matrix lse) (multiplicity (rhs lse))))
    (unless residual-p
      (copy! (rhs lse) residual)
      (gemm! -1.0 (matrix lse) solution 1.0 residual)
      (setq residual-p t))))

;;; Nonlinear problems

(defclass <nonlinear-problem> (<problem>)
  ((linearization :initarg :linearization
    :documentation "A function linearizing the problem.")
   (solution :initform nil :initarg :solution :documentation
		  "An approximation to the solution."))
  (:documentation "Class for nonlinear problems.  The linearization
contains a function returning a linear problem."))

(defmethod initialize-instance :after ((problem <nonlinear-problem>)
				       &key &allow-other-keys)
  (setf (slot-value problem 'linear-p) nil))

(defmethod linearize ((problem <nonlinear-problem>) solution)
  (funcall (slot-value problem 'linearization) solution))

(defmethod ensure-solution ((nlpb <nonlinear-problem>) blackboard)
  (ensure (getbb blackboard :solution)
	  (slot-value nlpb 'solution)))

(defclass <nlse> (<nonlinear-problem>)
  ()
  (:documentation "Class for nonlinear system of equations."))

(defun nlse (&rest args)
  "Constructs a standard NLSE."
  (apply #'make-instance '<nlse> args))

;;;; Testing

(defun test-problem ()
  (describe (lse :matrix #m(1.0) :rhs #m(1.0)))
  ;; problem for solving x^2=2 by Newton's method
  (let* ((nlse
	  (nlse :linearization
		#'(lambda (solution)
		    (assert solution)
		    (let ((x (vref solution 0)))
		      (lse :matrix (make-real-matrix (vector (* 2.0 x)))
			   :rhs (make-real-matrix (vector (+ 2.0 (* x x)))))))
		:solution #m(1.4))))
    (solve (blackboard :problem nlse :output 1)))
  )

;;; (test-problem)
(fl.tests:adjoin-test 'test-problem)
