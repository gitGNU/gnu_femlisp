;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; fe-interpolation.lisp - FE interpolation strategy
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

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Class definition
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass <fe-interpolation> (<fe-approximation>)
  ((coefficient :initarg :coefficient :documentation
	     "A coefficient determining the function to be interpolated."))
  (:documentation "This class implements adaptive finite element
interpolation of the given coefficient function as a variant of finite
element approximation."))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Adaption of the fe-approximation strategy
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun interpolate-function (as func &optional solution)
  "Interpolates @arg{func} with the ansatz-space @arg{as}."
  (ensure solution (make-ansatz-space-vector as))
  (with-slots (mesh) as
    (doskel (cell mesh)
      (let ((fe (get-fe as cell)))
	(unless (zerop (nr-of-inner-dofs fe))
	  (setf (vref solution (cell-key cell mesh))
		(interpolate-on-refcell
		 fe (compose func (curry #'local->global cell)))))))
    solution))

(defmethod approximate ((fe-strategy <fe-interpolation>) blackboard)
  "Interpolates a given function.  Does only Lagrange interpolation at the
moment."
  (with-items (&key mesh problem ansatz-space solution) blackboard
    (doskel (cell mesh)
      (let ((fe (get-fe ansatz-space cell)))
	(unless (zerop (nr-of-inner-dofs fe))
	  (let* ((patch (patch-of-cell cell mesh))
		 (coeff (get-coefficient (coefficients-of-patch patch problem)
					 (slot-value fe-strategy 'coefficient))))
	    (setf (vref solution (cell-key cell mesh))
		  (interpolate-on-refcell
		   fe #'(lambda (x)
			  (evaluate coeff (list :global (local->global cell x))))))))))))

;;;; Testing


(defun test-fe-interpolation ()
  (flet ((test (dim fe-class levels)
	   (let* ((domain (n-cube-domain dim))
		  (problem (make-instance
			    '<interpolation-problem> :domain domain
			    :components `((u ,(nr-of-components fe-class)))
			    :patch->coefficients
			    #'(lambda (patch)
				(princ patch) (terpri)
				(list (function->coefficient
				       'INITIAL
				       #'(lambda (x)
					(let ((phi (* 2 pi (aref x 0))))
					  (vector (cos phi) (sin phi)))))))))
		  (strategy (make-instance
			     '<fe-interpolation> :coefficient 'INITIAL
			     :indicator (make-instance '<uniform-refinement-indicator>)
			     :success-if `(>= :nr-levels ,levels) :plot-mesh nil :output t)))
	     (solve strategy (blackboard :problem problem :fe-class fe-class)))))
    (let ((bb (test 1 (lagrange-fe 4) 3)))
      (plot (getbb bb :solution)))
    (let ((bb (test 1 (lagrange-fe 4 :nr-comps 2) 3)))
      (plot (getbb bb :solution) :component 1)))
  )

;;; (test-fe-interpolation)
(fl.tests:adjoin-test 'test-fe-interpolation)