;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; bl-cell.lisp - Definitions of a boundary layer cell domain
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

(defun bottom-mapping (dim f &optional grad-f)
  (let ((dim-1 (- dim 1)))
    (make-instance
     '<special-function>
     :domain-dimension dim :image-dimension dim
     :evaluator
     #'(lambda (x)
	 (declare (type (simple-array double-float (*)) x))
	 (let* ((x1 (vector-slice x 0 dim-1))
		(value (funcall f x1))
		(result (copy-seq x)))
	   (setf (aref result dim-1) (* (1- (aref x dim-1)) value))
	   result))
     :gradient
     (and grad-f
	  #'(lambda (x)
	      (declare (type (simple-array double-float (*)) x))
	      (let* ((x1 (vector-slice x 0 dim-1))
		     (value (funcall f x1))
		     (gradient (scal! (1- (aref x dim-1))
				      (ensure-matlisp (funcall grad-f x1) :row)))
		     (result (eye dim)))
		(minject gradient result dim-1 0)
		(setf (mref result dim-1 dim-1) value)
		result))))))

(defun oscillating-boundary-domain (dim f &key grad-f (refinements 0) (upper t))
  "Returns a domain with an oscillating lower boundary at $x_n=-1$ where
the oscillation is defined by a scaling function $f$ with values in $\R^+$.
Usually, also $grad-f$ should be provided, because it makes possible an
enhanced domain approximation."
  (assert (zerop refinements) () "Please improve bl-patch-on-lower-boundary
before using a positive value for refinements.")
  (let* ((upper-cell (copy-skeleton
		      (refcell-refinement-skeleton (n-cube dim) refinements)))
	 (lower-cell (transformed-skeleton upper-cell :transformation
					   (bottom-mapping dim f grad-f))))
    (change-class (if upper
		      (skel-add! upper-cell lower-cell)
		      lower-cell)
		  '<domain>)))

(defun boundary-layer-cell-domain (dim f &key grad-f (refinements 0) (extensible t)
				   (upper t) &allow-other-keys)
  "Returns the domain generated by oscillating-boundary-domain with
identified lateral faces."
  (let ((dim-1 (1- dim))
	(domain (identify-unit-cell-faces
		 (oscillating-boundary-domain
		  dim f :grad-f grad-f :refinements refinements :upper upper)
		 :indices (range 0 (- dim 2)))))
    (when (and upper extensible)
      (let* ((upper-cell
	      (find-cell #'(lambda (cell)
			     (every #'(lambda (x) (= x 0.5)) (midpoint cell)))
			 domain))
	     (upper-side (aref (boundary upper-cell) (* 2 dim-1))))
	(multiple-value-bind (extension replacement)
	    (cube-extender upper-cell dim-1)
	  (identify-unit-cell-faces replacement :indices (range 0 (- dim 2)))
	  (setf (get-cell-property upper-side domain 'EXTENSION)
		extension)
	  (fl.mesh::ensure-secondary-information domain))))
    domain))

(defun sinusoidal-bl-cell (dim &rest rest &key (amplitude 0.15) &allow-other-keys)
  "Returns a boundary layer cell with a sinusoidally oscillating lower
boundary."
  (flet ((f (x)
	   (- 1.0 (* amplitude (reduce #'* x :key #'(lambda (xc) #I"sin(2.0*pi*xc)")))))
	 (grad-f (x)
	   (let ((result (make-double-vec (1- dim))))
	     (dotimes (i (1- dim) result)
	       (setf (aref result i)
		     (* (- amplitude)
			(let ((prod 1.0))
			  (dotimes (j (1- dim) prod)
			    (_f * prod (if (= i j)
					   #I"2.0*pi*cos(2.0*pi*x[j])"
					   #I"sin(2.0*pi*x[j])"))))))))))
    (apply #'boundary-layer-cell-domain dim #'f :grad-f #'grad-f rest)))

(defun xsinx-bl-cell (dim &rest rest &key (amplitude 0.4) &allow-other-keys)
  "Returns a boundary layer cell with a sinusoidally oscillating lower
boundary."
  (flet ((f (x)
	   (- 1.0 (* amplitude (reduce #'* x :key #'(lambda (xc) #I"xc*sin(2.0*pi*xc)")))))
	 (grad-f (x)
	   (let ((result (make-double-vec (1- dim))))
	     (dotimes (i (1- dim) result)
	       (setf (aref result i)
		     (* (- amplitude)
			(let ((prod 1.0))
			  (dotimes (j (1- dim) prod)
			    (_f * prod (if (= i j)
					   #I"sin(2.0*pi*x[j])+2.0*pi*x[j]*cos(2.0*pi*x[j])"
					   #I"x[j]*sin(2.0*pi*x[j])"))))))))))
    (apply #'boundary-layer-cell-domain dim #'f :grad-f #'grad-f rest)))

(defun spline-interpolated-bl-cell (heights)
  "Boundary which is interpolated from heights."
  (multiple-value-bind (f Df)
      (cubic-spline heights)
    (boundary-layer-cell-domain
     2 f :grad-f (compose (rcurry #'coerce 'double-vec) Df) :extensible nil)))

(defun bl-patch-on-lower-boundary (bl-domain patch)
  "Returns T if the patch is on the lower oscillating boundary."
  (let ((dim (dimension bl-domain)))
    (and (< (dimension patch) dim)
	 (minusp (aref (midpoint patch) (1- dim)))
	 (or (vertex? patch) (not (identified-p patch bl-domain))))))

(defun bl-patch-on-pellet-boundary (bl-domain patch)
  "Returns T if the patch is on the lower oscillating boundary."
  (let ((dim (dimension bl-domain)))
    (and (< (dimension patch) dim)
	 (notany (rcurry #'identified-p bl-domain) (vertices patch)))))

(defun bl-patch-on-upper-boundary (bl-domain patch)
  "Returns T if the patch is on the upper boundary."
  (let ((dim (dimension bl-domain)))
    (and (< (dimension patch) dim)
	 (plusp (aref (midpoint patch) (1- dim)))
	 (or (vertex? patch) (not (identified-p patch bl-domain))))))

(defun bl-patch-on-artificial-boundary (bl-domain patch)
  "Returns the artificial boundary on which the distributional source
acts."
  (zerop (aref (midpoint patch) (1- (dimension bl-domain)))))


;;;; Testing: (test-bl-cell)
(defun test-bl-cell ()
  (cubic-spline #(1.0 0.8))
  (spline-interpolated-bl-cell #(1.0 0.8))
  (check (spline-interpolated-bl-cell #(1.0 0.8)))
  (check (sinusoidal-bl-cell 2))
  (describe (sinusoidal-bl-cell 2))
  (let* ((domain (sinusoidal-bl-cell 2))
	 (mesh (uniformly-refined-hierarchical-mesh
		domain 0 :parametric :from-domain)))
    (extend mesh)
    (describe mesh)
    (domain-characteristics domain)
    (let ((count 0))
      (doskel (patch domain)
	(when (bl-patch-on-lower-boundary domain patch)
	  (format t "~A~%" patch)
	  (incf count)))
      (assert (= count 3))))
  )

(fl.tests:adjoin-test 'test-bl-cell)
