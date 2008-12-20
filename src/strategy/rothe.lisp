;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; rothe.lisp - Discretization for the Rothe method
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

(in-package :fl.strategy)

(defparameter *rothe-observe*
  (list (list "Step" "~4D"
	      (lambda (blackboard) (getbb blackboard :step)))
	*time-observe*
	(list "    Time" "~8,2F"
	      (lambda (blackboard)
		(getbb blackboard :model-time))))
  "Standard observe quantities for Rothe.")

(defclass <rothe> (<iteration>)
  ((model-time :accessor model-time :initarg :model-time :documentation
	       "Current time in the time-stepping scheme.")
   (time-step :accessor time-step :initarg :time-step)
   (scheme :reader time-stepping-scheme :initform :implicit-euler :initarg :scheme
	   :documentation "Time-stepping scheme,
e.g. @code{:implicit-euler} or @code{:crank-nicolson}.")
   ;; maybe we'll find something more automatic than the following later
   (stationary-success-if :initform nil :initarg :stationary-success-if)
   (stationary-failure-if :initform nil :initarg :stationary-failure-if)
   (plot :initform nil :initarg :plot)
   (fl.iteration::observe :initform *rothe-observe* :initarg :observe
    :documentation "Providing initform for <iteration> slot."))
  (:documentation "Rothe strategy for time-dependent problems.  The idea of
the Rothe method for solving @math{U_t +A U =f} is to do an ODE
time-stepping scheme in an infinite-dimensional function space.  Therefore,
in every time-step, the solution has to be approximated sufficiently well
in the space variable."))

(concept-documentation
 "Projection to the initial value is done by defining a problem performing
an L^2-projection on the finite element ansatz space.  This opens up the
possibility of handling distributional initial values as well.  The
extraction of the interpolation problem as well as the construction of the
time-stepping problem is problem-specific and therefore contained in
separate files.")

(defmethod initially ((rothe <rothe>) blackboard)
  (with-items (&key iv-ip-blackboard) blackboard
    (ensure iv-ip-blackboard (blackboard))
    (transfer-bb blackboard iv-ip-blackboard '(:fe-class :mesh))
    (with-items (&key problem success-if failure-if) iv-ip-blackboard
      (ensure problem (initial-value-interpolation-problem
		       rothe (getbb blackboard :problem)))
      (ensure success-if (slot-value rothe 'stationary-success-if))
      (ensure failure-if (slot-value rothe 'stationary-failure-if))
      (solve iv-ip-blackboard)
      (transfer-bb iv-ip-blackboard blackboard '(:solution :mesh))))
  (call-next-method))

(defmethod intermediate ((rothe <rothe>) blackboard)
  "Plots the solution initially and after each time step."
  (setf (getbb blackboard :model-time) (model-time rothe))
  (when (slot-value rothe 'plot)
    (plot (getbb blackboard :solution)))
  (call-next-method))

(defmethod next-step ((rothe <rothe>) blackboard)
  "Do one step by solving the stationary problem."
  (with-items (&key time-step-blackboard) blackboard
    (ensure time-step-blackboard (blackboard))
    (with-items (&key success-if failure-if fe-class problem mesh
		      ansatz-space solution)
	time-step-blackboard
      (ensure success-if (slot-value rothe 'stationary-success-if))
      (ensure failure-if (slot-value rothe 'stationary-failure-if))
      (transfer-bb blackboard time-step-blackboard '(:mesh :fe-class :plot-mesh))
      (setf problem (time-step-problem rothe (getbb blackboard :problem)))
      (setf ansatz-space (make-fe-ansatz-space fe-class problem mesh))
      (setf (get-property problem :old-solution) (getbb blackboard :solution))
      (setf solution (make-ansatz-space-vector ansatz-space))
      (copy! (getbb blackboard :solution) solution)  ; initial value
      (solve time-step-blackboard)
      (copy! solution (getbb blackboard :solution))
      (incf (model-time rothe) (time-step rothe))
      )))


(defun test-rothe ()
  ;; more test in the problem-adapted files
  (make-instance
   '<rothe> :model-time 0.0 :time-step 0.1
   :stationary-success-if '(> :nr-levels 3)
   :success-if '(>= :step 5)
   :output t :plot t)
  )

;;; (test-rothe)
(fl.tests:adjoin-test 'test-rothe)
