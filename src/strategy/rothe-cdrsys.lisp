;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; rothe-cdr.lisp - Rothe adaption for cdr problems
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

(defun cdrsys-initial-value-interpolation-coefficients (rothe problem coeffs)
  "Return a coefficient list for the initial value interpolation problem to
be solved at the beginning."
  (declare (ignore rothe))
  (map-coefficients
   #'(lambda (coeff-name coeff)
       (case coeff-name
	 (CONSTRAINT (values coeff-name coeff))
	 (FL.ELLSYS::INITIAL (values 'FL.ELLSYS::F coeff))
	 (FL.ELLSYS::C
	  (values coeff-name (constant-coefficient
			      (diagonal-sparse-tensor #m(1.0) (nr-of-components problem)))))
	 (t (values nil nil))))
   coeffs))
   
(defmethod initial-value-interpolation-problem ((rothe <rothe>) (problem <cdrsys-problem>))
  "Returns a stationary problem for interpolating the initial value."
  (make-instance
   '<cdrsys-problem>
   :domain (domain problem)
   :multiplicity (multiplicity problem)
   :coefficients
   (map-hash-table
    #'(lambda (patch coeffs)
	(values patch (cdrsys-initial-value-interpolation-coefficients rothe problem coeffs)))
    (coefficients problem))
   :linear-p t))

(defun cdrsys-time-step-coefficients (rothe coeffs ncomps)
  "Return a coefficient list for the stationary problem to be solved at
each time step."
  (let ((coeffs (copy-seq coeffs)))
    (let ((sigma (get-coefficient coeffs 'FL.ELLSYS::SIGMA))
	  (reaction (get-coefficient coeffs 'FL.ELLSYS::C))
	  (source (get-coefficient coeffs 'FL.ELLSYS::F)))
      (ensure sigma (constant-coefficient
		     (diagonal-sparse-tensor (ones 1) ncomps)))
      ;; source and reaction modification
      (when reaction
	(setf (getf coeffs 'FL.ELLSYS::C)
	      (make-instance
	       '<coefficient> :dimension (dimension reaction) :demands (demands reaction) :eval
	       #'(lambda (&rest args)
		    (let ((reaction (evaluate reaction args))
			  (sigma (evaluate sigma args))
			  (factor
			   (/ (ecase (time-stepping-scheme rothe)
				(:implicit-euler 1.0)
				(:crank-nicolson 2.0))
			      (time-step rothe))))
		      (axpy factor sigma reaction)))))
	(setf (getf coeffs 'FL.ELLSYS::F)
	      (make-instance
	       '<coefficient> :dimension (dimension source)
	       :demands (add-fe-parameters-demand (demands source) '(:old-solution))
	       :eval
	       #'(lambda (&rest args &key old-solution &allow-other-keys)
		   (let ((sigma (evaluate sigma args))
			 (source (evaluate source args)))
		     (gemm (/ (time-step rothe)) sigma old-solution 1.0 source))))))
      ;; inclusion of time variable
      (loop for (coeff-name coeff) on coeffs by #'cddr
	    for demands = (demands coeff)
	  collect coeff-name collect
	  (if (find :time demands)
	      (make-instance
	       '<coefficient> :dimension (dimension coeff) :demands demands
	       :eval #'(lambda (&rest args)
			 (evaluate coeff (list* :time (model-time rothe) args))))
	      coeff)))))

(defmethod time-step-problem ((rothe <rothe>) (problem <cdrsys-problem>))
  "Returns a stationary problem corresponding to a step of the Rothe
method."
  (make-instance
   '<cdrsys-problem>
   :domain (domain problem)
   :multiplicity (multiplicity problem)
   :coefficients
   (map-hash-table
    #'(lambda (patch coeffs)
	(values patch (cdrsys-time-step-coefficients rothe coeffs (nr-of-components problem))))
    (coefficients problem))))

(defun test-rothe-cdrsys ()

  (let* ((dim 1) (levels 1) (order 4) (end-time 0.1) (steps 64)
	 (problem (cdrsys-model-problem
		   dim 1 :initial #'(lambda (x) (vector (ensure-matlisp #I(sin(2*pi*x[0])))))
		   :reaction (diagonal-sparse-tensor #m(0.0) 1)
		   :source (constant-coefficient (vector #m(0.0)))))
	 (rothe (make-instance
		 '<rothe> :model-time 0.0 :time-step (/ end-time steps)
		 :stationary-success-if `(> :nr-levels ,levels)
		 :success-if `(>= :step ,steps)
		 :output t :plot nil
		 :observe (append *rothe-observe* (list *u_1/4-observe*)))))
    ;;(time-step-problem rothe problem))
    (let ((result
	   (iterate
	    rothe
	    (blackboard
	     :problem problem :fe-class (lagrange-fe order :nr-comps 1)
	     :plot-mesh nil :output t
	     ))))
      ;; the correct value is #I(exp(-0.1*4*pi*pi))=0.019296302911016777
      ;; but implicit Euler is extremely bad
      (assert (< (abs (- 2.1669732216e-02  ; value computed by CMUCL
			 (vref (aref (fe-value (getbb result :solution) #d(0.25))
				     0) 0)))
		 1.0e-12))))
  )

;;; (test-rothe-cdrsys)
(fl.tests:adjoin-test 'test-rothe-cdrsys)
