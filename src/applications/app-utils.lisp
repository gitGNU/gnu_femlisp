;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; app-utils.lisp
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

(in-package :fl.application)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Utilities
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass mutex-object ()
  ((lock :initform (bordeaux-threads:make-recursive-lock "my-lock"))
   (data :initarg :data)))

(defun mutex-wrap (object)
  (make-instance 'mutex-object :data object))

(defun mutex-unwrap (object)
  (slot-value object 'data))

(defmacro accessing-exclusively (bindings &body body)
  (if bindings
      (destructuring-bind ((name object) . rest)
          bindings
        `(let ((,name ,object))
           (bordeaux-threads:with-recursive-lock-held ((slot-value ,name 'lock))
             (with-slots ((,name data)) ,name
               (accessing-exclusively ,rest ,@body)))))
      `(locally ,@body)))

(defmethod average-coefficient (ansatz-space &key coefficient)
  (let* ((problem (problem ansatz-space))
	 (mesh (mesh ansatz-space))
	 (order (discretization-order (fe-class ansatz-space)))
	 (result (mutex-wrap nil)))
    (with-workers
        ((lambda (cell)
           (let* ((coeffs (coefficients-of-cell cell mesh problem))
                  (coeff-function (get-coefficient coeffs coefficient))
                  (factor-dims (mapcar #'dimension (factor-simplices cell)))
                  (qrule (gauss-rule factor-dims (1+ order))))
             (loop for lcoords across (integration-points qrule)
                   and weight across (integration-weights qrule) do
                     (let* ((global (local->global cell lcoords))
                            (Dphi (local->Dglobal cell lcoords))
                            (factor (* weight (abs (det Dphi))))
                            (coeff-ip (evaluate coeff-function (list :global global))))
                       (accessing-exclusively ((tensor result))
                         (if tensor
                             (axpy! factor coeff-ip tensor)
                             (setq tensor (scal factor coeff-ip))))
                       )))))
      (doskel (cell mesh :dimension :highest :where :surface)
        (work-on cell)))
    (mutex-unwrap result)))

(defmethod correction-tensor ((solution <ansatz-space-vector>) (rhs <ansatz-space-vector>))
  (let ((result (make-real-matrix (multiplicity solution) (multiplicity rhs))))
    (dovec ((entry key) solution)
      (gemm! 1.0 entry (vref rhs key) 1.0 result :tn))
    result))

(defun convert-elasticity-correction (mat)
  "Converts the (dim^2)x(dim^2) matrix returned as result of correction
tensor into an (dim x dim)-array with (dim x dim)-matrix entries."
  (let* ((dim (floor (sqrt (nrows mat))))
	 (result (make-array (list dim dim))))
    (dotimes (i dim)
      (dotimes (j dim)
	(setf (mref result i j) (make-real-matrix dim))))
    (dotuple (index (make-list 4 :initial-element dim))
      (setf (mref (mref result (elt index 0) (elt index 1))
		  (elt index 2) (elt index 3))
	    (mref mat
		  (+ (* dim (elt index 0)) (elt index 2))
		  (+ (* dim (elt index 1)) (elt index 3)))))
    result))

(defun effective-tensor (blackboard)
  (with-items (&key ansatz-space problem solution rhs) blackboard
    (and solution rhs
	 (etypecase problem
	   (<cdr-problem>
	    (m- (mref (average-coefficient ansatz-space :coefficient 'FL.ELLSYS::A)
		      0 0)
		(correction-tensor solution rhs)))
	   (<ellsys-problem> 
	    (m- (average-coefficient ansatz-space :coefficient 'FL.ELLSYS::A)
		(convert-elasticity-correction (correction-tensor solution rhs)))))
	 )))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Special routines
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun plot-diffusion (problem &key (refinements 0) (depth 2) parametric)
  "Plots the first component of the diffusion tensor."
  (plot problem :refinements refinements :depth depth :parametric parametric
	:coefficient 'FL.ELLSYS::A 
	:key (lambda (val) (mref (mref val 0 0) 0 0))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Variables
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defparameter *result* nil
  "Special variable used for storing the blackbboard of the last
computation.")

(defmacro storing (&body body)
  "Stores the result of @arg{body} in @var{*result*}."
  `(setq *result* (progn ,@body)))