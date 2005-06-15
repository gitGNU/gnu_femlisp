;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; crs.lisp
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;
;;; Copyright (C) 2003-2005 Nicolas Neuss, University of Heidelberg.
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

(in-package :fl.algebra)

(file-documentation
 "This is rather old code which is almost not used.  It should be
coordinated and combined with @file{ccs.lisp} in the @package{FL.MATLISP}
package.")

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; crs-pattern : sparse matrix pattern
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass crs-pattern ()
  ((nrows :initarg :nrows :accessor nrows :type fixnum
	  :documentation "Number of rows in the pattern.")
   (ncols :initarg :ncols :accessor ncols :type fixnum
	  :documentation "Number of columns in the pattern.")
   (nr-of-entries :initarg :nr-of-entries :type fixnum)
   (store-size :initarg :store-size :type fixnum)
   (row-starts :initarg :row-starts :type fixnum-vec)
   (col-inds :initarg :col-inds :type fixnum-vec)
   (offsets :initarg :offsets :type fixnum-vec))
  (:documentation "A CRS (compact row-ordered storage) pattern allowing for
identification, see @cite{(Neuss 1998)}, for use within a
@class{crs-matrix}."))

(defmethod initialize-instance :after
    ((crs-pat crs-pattern) &key pattern &allow-other-keys)
  "This is the crs-pattern constructor.  A sparse matrix of the form
  | * 0 0 0 a |
  | 0 a 0 0 0 |
can be described by its dimensions nrows=2, ncols=5 together with the pattern
 '( ((* . 0) (a . 4))  ((a . 1)) )
* means a non-identified value.  Other symbols can be used to identify entries."
  (when pattern
    (when (some #'(lambda (lst) (not (apply #'< (mapcar #'cdr lst))))
		pattern)
      (error "unsorted patterns are not allowed"))
    ;; initialize some slots from given pattern-list
    (with-slots (nr-of-entries store-size row-starts col-inds)
	crs-pat
      (let* ((flattened-pattern (flatten-1 pattern))
	     (current-offset 0)
	     (offsets
	      (loop with table = (make-hash-table)
		    for entry in flattened-pattern
		    collect (cond ((eq (car entry) '*)
				   (incf current-offset)
				   (1- current-offset))
				  ((gethash (car entry) table))
				  (t (setf (gethash (car entry) table) current-offset)
				     (incf current-offset)
				     (1- current-offset))))))
	(setf nr-of-entries (length offsets)
	      store-size current-offset
	      row-starts (coerce (cons 0 (loop for row in pattern
					       sum (length row) into rs
					       collect rs))
				 'fixnum-vec)
	      col-inds (map 'fixnum-vec #'cdr flattened-pattern)
	      (slot-value crs-pat 'offsets) (coerce offsets 'fixnum-vec))))))

(defun full-crs-pattern (nrows ncols)
  "Returns trivial rectangular crs-patterns."
  (let ((N (* nrows ncols)))
    (make-instance
     'crs-pattern
     :nrows nrows :ncols ncols
     :nr-of-entries N
     :store-size N
     :row-starts (coerce (loop for i from 0 upto nrows
			       collect (* i ncols))
			 'fixnum-vec)
     :offsets (coerce (range< 0 N) 'fixnum-vec)
     :col-inds (coerce (loop for i from 0 below N
			     collect (mod i ncols))
		       'fixnum-vec))))

(defun pattern->full-pattern (pattern)
  (full-crs-pattern (nrows pattern) (ncols pattern)))

(defgeneric shift-pattern (pattern shift)
  (:documentation "shift-pattern: This function shifts a pattern to
its actual offsets in the sparse graph."))

(defmethod shift-pattern ((pattern crs-pattern) shift)
  (with-slots (nrows ncols nr-of-entries store-size row-starts col-inds offsets)
      pattern
    (make-instance
     'crs-pattern :nrows nrows :ncols ncols
     :nr-of-entries nr-of-entries :store-size store-size 
     :row-starts row-starts :col-inds col-inds
     :offsets (map 'fixnum-vec #'(lambda (offset) (+ offset shift))
		   offsets))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; crs-matrix
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass crs-matrix (store-vector <matrix>)
  ((pattern :reader pattern :initarg :pattern :type crs-pattern
	    :documentation "The pattern of the CRS matrix.  This is kept
separate such that the same pattern can be used for many matrices."))
  (:documentation "The class @class{crs-matrix} combines a crs-pattern and
a value vector."))

(defmethod nrows ((A crs-matrix))
  (slot-value (pattern A) 'nrows))

(defmethod ncols ((A crs-matrix))
  (slot-value (pattern A) 'ncols))

(defmethod nr-of-entries ((A crs-matrix))
  (slot-value (pattern A) 'nr-of-entries))

(eval-when (:compile-toplevel :load-toplevel :execute)
  (defun crs-matrix (type)
    "Construct a CRS matrix with entries of @arg{type}."
    (fl.amop:find-programmatic-class
     (list (find-class 'crs-matrix) (store-vector type)))))

(defmethod initialize-instance :after ((crs crs-matrix) &key &allow-other-keys)
  (assert (typep crs 'store-vector))
  (with-slots (store pattern) crs
    (unless (slot-boundp crs 'store)
      (setf store (make-array (slot-value pattern 'store-size)
			      :element-type (element-type crs)
			      :initial-element (coerce 0 (element-type crs)))))))

(defun make-crs-matrix (pattern store)
  "make-crs-matrix: crs-matrix constructor."
  (unless (= (slot-value pattern 'store-size) (length store))
    (error "pattern does not fit with value vector"))
  (make-instance (crs-matrix 'double-float) :pattern pattern :store store))

(defun make-full-crs-matrix (nrows ncols)
  (make-crs-matrix (full-crs-pattern nrows ncols)
		   (make-double-vec (* nrows ncols))))

(defmethod mref ((A crs-matrix) i j)
  (with-slots (pattern store) A
    (with-slots (row-starts col-inds offsets) pattern
    (loop for k from (aref row-starts i) below (aref row-starts (1+ i))
	  do (if (= (aref col-inds k) j)
		 (return (aref store (aref offsets k))))
	  finally (return 0.0)))))

(defmethod (setf mref) (val (A crs-matrix) i j)
  (with-slots (pattern store) A
    (with-slots (row-starts col-inds offsets) pattern
      (loop for k from (aref row-starts i) below (aref row-starts (1+ i))
	    do (if (= (aref col-inds k) j)
		   (return
		     (setf (aref store (aref offsets k)) val)))
	    finally (error "i/j not in pattern")))))

;;; transformation to matlisp format
(defun crs->matlisp (A)
  (let ((mlmat (make-real-matrix (nrows A) (ncols A))))
    (dotimes (i (nrows A))
      (dotimes (j (ncols A))
	(setf (mref mlmat i j) (mref A i j))))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; vector blas operations for crs-matrix
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; These are given because crs-matrix is a store-vector.

(defmethod make-analog ((mat crs-matrix))
  (make-instance (crs-matrix (element-type mat))
		 :pattern (pattern mat)))

;;;; Testing

(defun test-crs ()
  (let* ((crs-pat (make-instance 'crs-pattern
				 :nrows 2 :ncols 5
				 :pattern '( ((* . 0) (a . 4))  ((a . 1)) )))
	 (A (make-instance (crs-matrix 'double-float)
			   :pattern crs-pat :store #d(2.0 2.0)))
	 (B (make-instance (crs-matrix 'double-float)
			   :pattern crs-pat :store #d(2.0 1.0))))
    (describe (m+ A B))
    (describe (m- A B))
    (describe A)
    (x<-0 A)
    (describe A)
    )
  (describe (full-crs-pattern 2 2))
  nil)

;;; (fl.algebra::test-crs)
(fl.tests:adjoin-test 'test-crs)