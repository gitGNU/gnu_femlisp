;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; solve.lisp
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

(in-package :iteration)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; General solver class
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


(defclass <solver> ()
  ((maxsteps :reader maxsteps :initform nil :initarg :maxsteps)
   (threshabs :initform nil :initarg :threshabs)
   (threshold :reader threshold :initform nil :initarg :threshold)
   (reduction :reader reduction :initform nil :initarg :reduction)
   (residual-norm :reader residual-norm :initform #'norm
		  :initarg :residual-norm)
   (success-if :reader success-if :initform nil :initarg :success-if)
   (failure-if :reader failure-if :initform nil :initarg :failure-if)
   (output :reader output :initform nil :initarg :output))
  (:documentation "The base class of linear, nonlinear and whatever
iterative solvers."))

(defgeneric solve (solver &rest parameters)
  (:documentation "Solve a problem specified through the parameter list."))

;;; Output

(defclass <iteration-output> ()
  ((initial :reader output-initial :initarg :initial :initform nil :type (or null function))
   (before-step :reader output-before-step :initarg :before-step
		:initform nil :type (or null function))
   (after-step :reader output-after-step :initarg :after-step
	       :initform nil :type (or null function))
   (final :reader output-final :initarg :final :initform nil :type (or null function)))
  (:documentation "Class for steering the output of an iterative solver."))

#+(or)
(defparameter *solver-output-indentation* 0
  "Indentation.  Should be respected by all instances of solver-output.")

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Testing success and failure
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun test-condition (condition &key defnorm initial-defnorm previous-defnorm step)
  "Condition can be an expression involving the keywords :step, :defnorm,
:reduction, :step-reduction which is evaluated."
  (let ((reduction (and initial-defnorm defnorm (not (zerop initial-defnorm))
			(/ defnorm initial-defnorm)))
	(step-reduction (and defnorm previous-defnorm (not (zerop previous-defnorm))
			     (/ defnorm previous-defnorm))))
    (let ((condition 
	   (sublis `((:step . ,step) (:defnorm . ,defnorm) (:reduction . ,reduction)
		     (:step-reduction . ,step-reduction))
		   condition)))
      (unless (rmember NIL condition)
	;; the condition is valid and its value is returned
	(eval condition)))))

;;;; Testing: (test-solve)
(defun test-solve ()
  (test-condition '(and (> :step 5) (< :defnorm 1.0e-8))
		  :defnorm 1.0e-9 :initial-defnorm 1 :previous-defnorm 0.5 :step 10)
  )
(tests::adjoin-femlisp-test 'test-solve)
