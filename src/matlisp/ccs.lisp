;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; ccs.lisp - Compressed column storage scheme
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

(defclass store-vector (<vector>)
  ((store :reader store :initarg :store :documentation
	  "The vector entries."))
  (:documentation "This mixin yields vector behaviour for a class
containing a store.  The store is a unifom array with elements of a certain
type which can be determined by the funtion @function{element-type}.  It
often is but does not have to be equal to the type of scalars for this
vector which can be obtained by calling the function
@function{scalar-type}."))

(eval-when (:compile-toplevel :load-toplevel :execute)
  (defun store-vector (type)
    (assert (subtypep type 'number))
    (let ((class-name (intern (format nil "~A" (list 'store-vector type))
			      :fl.matlisp)))
      (or (find-class class-name nil)
	  (prog1
	      (eval `(defclass ,class-name (store-vector) ()))
	    (eval `(defmethod element-type ((vector ,class-name))
		    ',type))
	    (eval `(defmethod scalar-type ((vector ,class-name))
		    ',type)))))))

(defmethod initialize-instance :after ((vec store-vector) &key &allow-other-keys)
  "Coerce the store to the correct type."
  (coerce (store vec) `(simple-array ,(element-type vec) (*))))

(defclass ccs-pattern ()
  ((nrows :initarg :nrows :accessor nrows :type fixnum
	  :documentation "Number of rows in the pattern.")
   (ncols :initarg :ncols :accessor ncols :type fixnum
	  :documentation "Number of columns in the pattern.")
   (column-starts :initarg :column-starts
		  :documentation "Vector with start indices of columns.")
   (row-indices :initarg :row-indices
		:documentation "Vector with row indices."))
  (:documentation "A CCS (compressed column storage) pattern.  We work with
integer vectors that do not have to be copied for calls to C programs."))

(defmethod number-nonzero-entries ((pattern ccs-pattern))
  (with-slots (ncols column-starts) pattern
    (aref column-starts ncols)))
  
(defmethod initialize-instance :after ((pattern ccs-pattern) &key &allow-other-keys)
  "Check pattern and coerce vectors, if necessary."
  (with-slots (ncols column-starts row-indices) pattern
    (assert (= (length column-starts) (1+ ncols)))
    (assert (= (length row-indices) (number-nonzero-entries pattern)))
    (setq column-starts (coerce column-starts 'int-vec))
    (setq row-indices (coerce row-indices 'int-vec))))

(defclass ccs-matrix (<matrix>)
  ((pattern :reader pattern :initarg :pattern :type ccs-pattern
	  :documentation "CCS pattern."))
  (:documentation "A CCS (compressed column storage) matrix.  This is an
abstract class which is made concrete by mixing it with a store-vector."))

(defmethod nrows ((ccs ccs-matrix))
  (slot-value (pattern ccs) 'nrows))

(defmethod ncols ((ccs ccs-matrix))
  (slot-value (pattern ccs) 'ncols))

(eval-when (:compile-toplevel :load-toplevel :execute)
  (defun ccs-matrix (type)
    "Construct a CCS matrix with entries of @arg{type}."
    (fl.amop:find-programmatic-class
     (list (find-class 'ccs-matrix) (store-vector type)))))

(defmethod initialize-instance :after ((ccs ccs-matrix) &key &allow-other-keys)
  "Check the length of the store."
  (assert (typep ccs 'store-vector))
  (assert (= (length (store ccs))
	     (number-nonzero-entries (pattern ccs)))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; GEMV!
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmethod gemv-nn! ((alpha number) (x ccs-matrix) (y standard-matrix)
		     (beta number) (z standard-matrix))
  ;;(declare (optimize (speed 3) (space 0) (debug 0) (safety 0)))
  (unless (= beta 1) (scal! beta z))
  (let* ((store (store x))
	 (pattern (pattern x))
	 (column-starts (slot-value pattern 'column-starts))
	 (row-indices (slot-value pattern 'row-indices)))
    (if (= (ncols x) 1)
	(dotimes (i (ncols x))
	  (dotimes (offset (aref column-starts i) (aref column-starts (1+ i)))
	    (let ((k (aref row-indices offset))
		  (factor (* alpha (aref store offset))))
	      (incf (vref z i) (* factor (vref y k))))))
	(dotimes (i (ncols x))
	  (dotimes (offset (aref column-starts i) (aref column-starts (1+ i)))
	    (let ((k (aref row-indices offset))
		  (factor (* alpha (aref store offset))))
	      (dotimes (l (ncols y))
		(incf (mref z i l) (* factor (mref y k l))))))))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; GESV!
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defvar *default-ccs-solver*
  #+superlu 'fl.alien::superlu
  #-superlu #+umfpack 'fl.alien::umfpack
  #-(or superlu umfpack) nil
  "Default solver for the CCS format.  At the moment this can be UMFPACK or
SuperLU.")

(defmethod gesv! ((mat ccs-matrix) (vec standard-matrix))
  "Solve the system by calling an external sparse solver."
  (with-slots (nrows ncols column-starts row-indices)
      (pattern mat)
    (assert (= nrows (nrows vec)))
    (funcall *default-ccs-solver*
	     nrows ncols (number-nonzero-entries (pattern mat))
	     column-starts row-indices (store mat)
	     (ncols vec) (store vec) (store vec)))
  vec)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Test direct solvers on the CCS scheme
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun five-point-stencil-matrix (nx ny)
  "Generate a CCS matrix for the five-point stencil on a grid with the given number of active grid points in each dimension."
  (declare (type (integer 0 10000) nx ny))
  (let* ((nrows (* nx ny))
	 (ncols (* nx ny))
	 (nnz (- (* 5 nrows) (* 2 (+ nx ny))))
	 (column-starts (make-array (1+ nrows) :element-type '(signed-byte 32)))
	 (row-indices (make-array nnz :element-type '(signed-byte 32)))
	 (store (make-array nnz :element-type 'double-float))
	 (pos 0))
    (declare (type (integer 0 100000000) pos))
    (dotimes (j ny)
      (declare (type (integer 0 10000) j))
      (dotimes (i nx)
	(let ((k (+ i (* j nx))))
	  (declare (optimize speed))
	  (setf (aref column-starts k) pos)
	  (flet ((connect (l value)
		   (setf (aref store pos) value)
		   (setf (aref row-indices pos) l)
		   (incf pos)))
	    (when (plusp j) (connect (- k nx) -1.0))
	    (when (plusp i) (connect (- k 1) -1.0))
	    (connect k 4.0)
	    (when (< i (1- nx)) (connect (+ k 1) -1.0))
	    (when (< j (1- ny)) (connect (+ k nx) -1.0))))))
    (assert (= pos nnz))
    (setf (aref column-starts nrows) pos)
    ;; return matrix
    (make-instance
     (ccs-matrix 'double-float)
     :pattern (make-instance 'ccs-pattern :nrows nrows :ncols ncols
			     :column-starts column-starts
			     :row-indices row-indices)
     :store store)))

(defun direct-solver-performance-test (solver n)
  (when (evenp n) (setf n (1- n)))
  (let* ((ccs (five-point-stencil-matrix n n))
	 (rhs (ones (nrows ccs) 1)))
    (time
     (let ((*default-ccs-solver* solver))
       (progn (gesv! ccs rhs)
	      (/ (vref rhs (floor (nrows ccs) 2))
		 (expt (1+ n) 2)))))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Testing
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun test-ccs ()
  #+umfpack (direct-solver-performance-test 'fl.alien::umfpack 400)
  #+superlu (direct-solver-performance-test 'fl.alien::superlu 200)
  (make-instance (store-vector 'double-float) :store #d(0.0))
  (let* ((pattern (make-instance
		   'ccs-pattern :nrows 1 :ncols 1
		   :column-starts #(0 1)
		   :row-indices #(0)))
	 (ccs (make-instance
	       (ccs-matrix 'double-float) :pattern pattern :store #d(2.0)))
	 (rhs #m(1.0)))
    (gesv! ccs rhs)))

(fl.tests:adjoin-test 'test-ccs)


