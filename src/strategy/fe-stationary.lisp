;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; stationary.lisp - Solution strategy for stationary PDE problems
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

(in-package :fl.strategy)

(defvar *mentries-observe*
  (list " MENTRIES" "~9D"
	#'(lambda (blackboard)
	    (acond
	     ((getbb blackboard :discretized-problem)
	      (total-entries (matrix it)))
	     (t nil))))
  "Observe entries for the size of the matrix.")

(defvar *stationary-fe-strategy-observe*
  *fe-approximation-observe*
   "Standard observe quantities for stationary fe-strategies.")

(defclass <stationary-fe-strategy> (<fe-approximation>)
  ((fl.iteration::observe :initform *stationary-fe-strategy-observe*)
   (solver :initarg :solver :documentation
	   "The solver to be used for solving the discretized systems."))
  (:documentation "This class describes some iterative finite element
solution strategies for continuous, stationary PDE problems."))

(defmethod initially :after ((strategy <stationary-fe-strategy>) blackboard)
  (unless (slot-boundp strategy 'solver)
    (setf (slot-value strategy 'solver)
	  (select-solver (getbb blackboard :problem) blackboard))))

(defmethod approximate ((fe-strategy <stationary-fe-strategy>) blackboard)
  "Ensures accuracy of the solution and the error estimate."
  (with-items (&key interior-matrix interior-rhs discretized-problem
		    solver-blackboard solution residual output)
      blackboard
    ;; assemble (better would be a local assembly)
    (setf interior-matrix nil interior-rhs nil)
    (fe-discretize blackboard)
    ;; improve approximation by solving
    (setf solver-blackboard (blackboard :problem discretized-problem
					:solution solution :residual residual
					:output output))
    (solve (slot-value fe-strategy 'solver) solver-blackboard)
    (setf solution (getbb solver-blackboard :solution))))

