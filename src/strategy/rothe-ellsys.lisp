;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; rothe-ellsys.lisp - Rothe adaption for elliptic systems
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;
;;; Copyright (C) 2005 Nicolas Neuss, University of Heidelberg.
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

(defun ellsys-initial-value-interpolation-coefficients (rothe problem coeffs)
  "Return a coefficient list for the initial value interpolation problem to
be solved at the beginning."
  (declare (ignore rothe))
  (let ((initial (get-coefficient coeffs 'INITIAL))
	(sigma (get-coefficient coeffs 'SIGMA))
	(constraints (remove-if-not #'constraint-p coeffs)))
    (if initial
	(list* (copy-coefficient initial :name 'FL.ELLSYS::F)
	       (if sigma
		   (copy-coefficient sigma :name 'FL.ELLSYS::R)
		   (constant-coefficient 'FL.ELLSYS::R
					 (diagonal-sparse-tensor
					  #m(1.0) (nr-of-components problem))))
	       constraints)
	constraints)))
   
(defmethod initial-value-interpolation-problem ((rothe <rothe>) (problem <ellsys-problem>))
  "Returns a stationary problem for interpolating the initial value."
  (make-instance
   '<ellsys-problem>
   :domain (domain problem)
   :multiplicity (multiplicity problem)
   :components (components problem)
   :coefficients
   (map-hash-table
    #'(lambda (patch coeffs)
	(values patch (ellsys-initial-value-interpolation-coefficients
		       rothe problem coeffs)))
    (coefficients problem))
   :linear-p t))

(defun ellsys-time-step-coefficients (rothe coeffs)
  "Return a coefficient list for the stationary problem to be solved at
each time step."
  (let ((sigma (get-coefficient coeffs 'SIGMA)))
    (if sigma
	(list* (make-instance '<coefficient> :name 'FL.ELLSYS::R :demands ()
			      :eval (lambda (&rest args)
				      (scal (/ (ecase (time-stepping-scheme rothe)
						 (:implicit-euler 1.0)
						 (:crank-nicolson 2.0))
					       (time-step rothe))
					    (evaluate sigma args))))
	       (make-instance '<coefficient> :name 'FL.ELLSYS::F
			      :demands '((:fe-parameters :old-solution))
			      :eval (lambda (&rest args &key old-solution &allow-other-keys)
				      (let ((old-solution (first old-solution)))
					(gemm (/ (time-step rothe)) (evaluate sigma args) old-solution
					      0.0 old-solution))))
	       coeffs)
	coeffs)))

(defmethod time-step-problem ((rothe <rothe>) (problem <ellsys-problem>))
  "Returns a stationary problem corresponding to a step of the Rothe
method."
  (make-instance
   '<ellsys-problem>
   :domain (domain problem)
   :components (components problem)
   :multiplicity (multiplicity problem)
   :coefficients
   (map-hash-table
    #'(lambda (patch coeffs)
	(values patch (ellsys-time-step-coefficients rothe coeffs)))
    (coefficients problem))))

(defvar *u_1/4-observe*
  (list (format nil "~19@A" "u(midpoint)") "~19,10,2E"
	#'(lambda (blackboard)
	    (let* ((sol (getbb blackboard :solution))
		   (dim (dimension (mesh (ansatz-space sol))))
		   (val (fe-value sol (make-double-vec dim 0.25))))
	      (vref (aref val 0) 0)))))

(defun test-rothe-ellsys ()
  ;; the same test as in rothe-cdr
  (let* ((dim 1) (levels 1) (order 4) (end-time 0.1) (steps 64)
	 (problem (ellsys-model-problem
		   dim '((u 1))
		   :initial (lambda (x) (vector (ensure-matlisp #I(sin(2*pi*x[0])))))
		   :sigma (diagonal-sparse-tensor (vector (eye 1)))
		   :a (isotropic-diffusion 1 1.0)
		   :r (diagonal-sparse-tensor (vector #m(0.0)))
		   :f (vector #m(0.0))
		   :dirichlet (constraint-coefficient 1 1)))
	 (rothe (make-instance
		 '<rothe> :model-time 0.0 :time-step (/ end-time steps)
		 :stationary-success-if `(> :nr-levels ,levels)
		 :success-if `(>= :step ,steps)
		 :output t :plot t
		 :observe (append *rothe-observe* (list *u_1/4-observe*)))))
    ;;(time-step-problem rothe problem)
    (let ((result
	   (iterate
	    rothe
	    (blackboard
	     :problem problem :fe-class (lagrange-fe order :nr-comps 1)
	     :plot-mesh nil :output :all
	     ))))
      ;; the correct value is #I(exp(-0.1*4*pi*pi))=0.019296302911016777
      ;; but implicit Euler is extremely bad
      (assert (< (abs (- 2.1669732216e-02  ; value computed by CMUCL
			 (vref (aref (fe-value (getbb result :solution) #d(0.25))
				     0) 0)))
		 1.0e-12))
      ))
  )

;;; (test-rothe-ellsys)
(fl.tests:adjoin-test 'test-rothe-ellsys)
