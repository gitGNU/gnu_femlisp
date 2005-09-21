;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; compat.lisp - Compatibility with old Femlisp
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

(in-package :fl.matlisp)

(defmethod make-analog ((arr array))
  (make-array (array-dimensions arr) :element-type (array-element-type arr)))

;;; Define the Femlisp vector blas commands for matlisp matrices.

;;; new matlisp commands

(defun det-from-lr (lr pivot)
  "This routine computes the determinant using a given LR decomposition."
  (loop with prod of-type double-float = 1.0
	for i below (nrows lr) do
	(setq prod (* prod (mref lr i i)))
	(unless (= (aref pivot i) (1+ i)) (setq prod (- prod)))
	finally (return prod)))

(defun det (mat)
  "Returns the determinant of the square matrix @arg{mat}."
  (declare (values double-float))
  (if (or (zerop (nrows mat)) (zerop (ncols mat)))
      1.0  ; yields correct value for volume in the vertex case
      (multiple-value-bind (lr pivot)
	  (getrf! (copy mat))
	(det-from-lr lr pivot))))

(defun area-of-span (mat)
  "Computes the volume spanned by the columns of @arg{mat}."
  (sqrt (abs (det (m*-tn mat mat)))))

(defmethod mzerop ((mat standard-matrix) &optional (threshold 0.0))
  (every #'(lambda (x) (<= (abs x) threshold))
	 (slot-value mat 'store)))

(defmethod midentity-p ((x number) &optional (threshold 0.0))
  (<= (abs (- x 1.0)) threshold))

(defmethod midentity-p (mat &optional (threshold 0.0))
  (for-each-row-key
   #'(lambda (rk)
       (for-each-key-in-row
	#'(lambda (ck)
	    (let ((entry (mref mat rk ck)))
	      (unless (if (eql rk ck)
			  (midentity-p entry threshold)
			  (mzerop entry threshold))
		(return-from midentity-p
		  (values nil (cons rk ck) entry)))))
	    mat rk))
       mat)
   t)

(defmethod submatrix ((mat standard-matrix) &key row-indices col-indices)
  (unless row-indices (setq row-indices (range< 0 (nrows mat))))
  (unless col-indices (setq col-indices (range< 0 (ncols mat))))
  (let ((result (make-instance (standard-matrix (element-type mat))
			       :nrows (length row-indices)
			       :ncols (length col-indices))))
    (loop for i from 0 and row-ind in row-indices do
      (loop for j from 0 and col-ind in col-indices do
	    (setf (mref result i j) (mref mat row-ind col-ind))))
    result))

(defmethod vector-slice ((vec vector) offset size)
  "Provides a convenient shorthand for constructing a displaced
double-float array."
  (let ((result (make-array size :element-type (array-element-type vec))))
    (dotimes (i size result)
      (setf (aref result i) (aref vec (+ i offset))))))

(defmethod make-image-vector-for ((mat standard-matrix) &optional (multiplicity 1))
  (make-instance (standard-matrix (element-type mat))
		 :nrows (nrows mat) :ncols multiplicity))

(defmethod make-domain-vector-for ((mat standard-matrix) &optional (multiplicity 1))
  (make-instance (standard-matrix (element-type mat))
		 :nrows (ncols mat) :ncols multiplicity))

;;; old interface matrix-vector multiplication (deprecated)

(defmethod x<-Ay (x (A standard-matrix) y)
  (gemm! 1.0 A (ensure-matlisp y) 0.0 (ensure-matlisp x))
  x)
(defmethod x+=Ay (x (A standard-matrix) y)
  (gemm! 1.0 A (ensure-matlisp y) 1.0 (ensure-matlisp x))
  x)
(defmethod x-=Ay (x (A standard-matrix) y)
  (gemm! -1.0 A (ensure-matlisp y) 1.0 (ensure-matlisp x))
  x)

;;; Multiplication of Matlisp matrices and Lisp vectors

(defmethod m*-product-instance ((A standard-matrix) (x vector))
  "Returns a CL vector as result of matrix * CL vector."
  (zero-vector (nrows A) (array-element-type x)))

(defmethod gemm-nn! (alpha (A standard-matrix) (x vector) beta (y vector))
  (gemm-nn! alpha A (ensure-matlisp x) beta (ensure-matlisp y))
  y)

(defmethod getrs! ((lr standard-matrix) (b vector) &optional ipiv)
  (getrs! lr (ensure-matlisp b) ipiv)
  b)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Testing
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun test-matlisp-vector-combination ()
  (let ((A #m((1.0 2.0) (3.0 4.0)))
	(B (make-real-matrix 2 2)))
    (x<-0 A)
    (fill! B 1.0)
    (axpy 0.5 B A)
    (mzerop #(0.0))
    (mzerop #m((0.0)))
    (area-of-span (mrandom 5 3)))
  (let* ((dim 2)
	 (x (make-double-vec dim))
	 (y (make-double-vec dim 1.0)))
    (m* (ensure-matlisp x)
	(transpose (ensure-matlisp y)))
    (mapcar #'class-of
	    (list (eye 2) (ensure-matlisp y) (ensure-matlisp x)))
    (gemm 1.0 (eye 2) (ensure-matlisp x) 1.0 (ensure-matlisp y)))
  )

;;; (fl.matlisp::test-matlisp-vector-combination)
(fl.tests:adjoin-test 'test-matlisp-vector-combination)


