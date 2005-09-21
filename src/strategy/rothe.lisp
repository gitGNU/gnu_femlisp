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
	      #'(lambda (blackboard) (getbb blackboard :step)))
	*time-observe*
	(list "    Time" "~8,2F"
	      #'(lambda (blackboard)
		  (getbb blackboard :model-time))))
  "Standard observe quantities for Rothe.")

(defclass <rothe> (<iteration>)
  ((model-time :accessor model-time :initarg :model-time :documentation
	       "Current time in the time-stepping scheme.")
   (time-step :accessor time-step :initarg :time-step)
   ;; maybe we'll find something more automatic than the following later
   (stationary-success-if :initform nil :initarg :stationary-success-if)
   (stationary-failure-if :initform nil :initarg :stationary-failure-if)
   (plot :initform nil :initarg :plot)
   (fl.iteration::observe :initform *rothe-observe* :initarg :observe
    :documentation "Providing initform for <iteration> slot."))
  (:documentation "Rothe strategy for time-dependent problems.  The idea of
the Rothe method for solving $$U_t +A U =f$$ is to do an ODE time-stepping
scheme in an infinite-dimensional function space.  Therefore, in every
time-step, the solution has to be approximated sufficiently well in the
space variable."))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Interpolation of the initial value
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; We do this by defining a problem performing an L^2-projection on the
;;; finite element ansatz space.  This opens up the possibility of handling
;;; distributional initial values as well.

(defun map-coefficients (func coeffs)
  "Maps a given coefficient list COEFFS into a new coefficient list.  FUNC
takes coefficient name and coefficient and returns two values for new
coefficient name and coefficient.  If the first value returned is NIL, this
coefficient is not collected."
  (loop
   for (coeff-name coeff) on coeffs by #'cddr
   for (new-name new-coeff) = (multiple-value-list (funcall func coeff-name coeff))
   when new-name collect new-name and collect new-coeff))

(defun cdr-initial-value-interpolation-coefficients (rothe coeffs)
  "Return a coefficient list for the initial value interpolation problem to
be solved at the beginning."
  (declare (ignore rothe))
  (map-coefficients
   #'(lambda (coeff-name coeff)
       (case coeff-name
	 (FL.CDR::CONSTRAINT (values coeff-name coeff))
	 (FL.CDR::INITIAL (values 'FL.CDR::SOURCE coeff))
	 (FL.CDR::REACTION (values coeff-name (constant-coefficient 1.0)))
	 (t (values nil nil))))
   coeffs))
   
(defmethod initial-value-interpolation-problem ((rothe <rothe>) (problem <cdr-problem>))
  "Returns a stationary problem for interpolating the initial value."
  (make-instance
   '<cdr-problem>
   :domain (domain problem)
   :multiplicity (multiplicity problem)
   :coefficients
   (map-hash-table
    #'(lambda (patch coeffs)
	(values patch (cdr-initial-value-interpolation-coefficients rothe coeffs)))
    (coefficients problem))
   :linear-p t))

(defun cdr-time-step-coefficients (rothe coeffs)
  "Return a coefficient list for the stationary problem to be solved at
each time step."
  ;; source and reaction modification
  (whereas ((coeff (getf coeffs 'FL.CDR::REACTION)))
    (setf (getf coeffs 'FL.CDR::REACTION)
	  (make-instance
	   '<coefficient> :dimension (dimension coeff) :demands (demands coeff) :eval
	   #'(lambda (&rest args)
	       (+ (evaluate coeff args) (/ (time-step rothe))))))
    (let ((coeff (getf coeffs 'FL.CDR::SOURCE)))
      (setf (getf coeffs 'FL.CDR::SOURCE)
	    (make-instance
	     '<coefficient> :dimension (dimension coeff)
	     :demands (adjoin :solution  (demands coeff)) :eval
	     #'(lambda (&rest args &key solution &allow-other-keys)
		 (axpy (/ (time-step rothe)) solution (evaluate coeff args)))))))
  ;; inclusion of time variable
  (loop for (coeff-name coeff) on coeffs by #'cddr
	for demands = (demands coeff)
	collect coeff-name collect
	(if (find :time demands)
	    (make-instance
	     '<coefficient> :dimension (dimension coeff) :demands demands
	     :eval #'(lambda (&rest args)
		       (evaluate coeff (list* :time (model-time rothe) args))))
	    coeff)))

(defmethod time-step-problem ((rothe <rothe>) (problem <cdr-problem>))
  "Returns a stationary problem corresponding to a step of the Rothe
method."
  (assert (linear-p problem))
  (make-instance
   '<cdr-problem>
   :domain (domain problem)
   :multiplicity (multiplicity problem)
   :coefficients
   (map-hash-table
    #'(lambda (patch coeffs)
	(values patch (cdr-time-step-coefficients rothe coeffs)))
    (coefficients problem))
   :linear-p t))

(defmethod initially ((rothe <rothe>) blackboard)
  (with-items (&key iv-ip-blackboard fe-class) blackboard
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
  (with-items (&key step time-step-blackboard) blackboard
    (ensure time-step-blackboard (blackboard))
    (with-items (&key success-if failure-if fe-class problem mesh
		      ansatz-space solution)
	time-step-blackboard
      (ensure success-if (slot-value rothe 'stationary-success-if))
      (ensure failure-if (slot-value rothe 'stationary-failure-if))
      (transfer-bb blackboard time-step-blackboard '(:mesh :fe-class :plot-mesh))
      (setf problem (time-step-problem rothe (getbb blackboard :problem)))
      (setf ansatz-space (make-fe-ansatz-space fe-class problem mesh))
      (setf solution (make-ansatz-space-vector ansatz-space))
      (copy! (getbb blackboard :solution) solution)
      (solve time-step-blackboard)
      (copy! solution (getbb blackboard :solution))
      (incf (model-time rothe) (time-step rothe))
      )))

(defun test-rothe ()

  (let* ((dim 1) (levels 4) (order 2)
	 (problem (cdr-model-problem
		   dim :initial #'(lambda (x) #I(sin(2*pi*x[0]^^2)))
		   :reaction (constant-coefficient 0.0)
		   :source (constant-coefficient #m(0.0))))
	 (rothe (make-instance
		 '<rothe> :model-time 0.0 :time-step 0.01
		 :stationary-success-if `(> :nr-levels ,levels)
		 :success-if '(>= :step 20)
		 :output t :plot t)))
    (let ((result
	   (iterate rothe (blackboard
			   :problem problem :fe-class (lagrange-fe order)
			   :plot-mesh nil :output t))))
      (assert (< (abs (- 0.3833184697778781  ; value computed by CMUCL
			 (vref (aref (fe-value (getbb result :solution) #d(0.5))
				     0) 0)))
		 1.0e-12)))
    )
  )

;;; (test-rothe)
(fl.tests:adjoin-test 'test-rothe)
