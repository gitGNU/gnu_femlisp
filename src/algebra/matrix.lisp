;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; matrix.lisp - Matrix class and interface
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

(in-package :algebra)

(defclass <matrix> ()
  ()
  (:documentation "Every matrix should belong to this class.
Unfortunately, most matrices don't because they are defined in separate
packages."))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Matrix and matrix/vector interface
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; cell access
(defgeneric* mat-ref (mat i j))
(defgeneric* (setf mat-ref) (val mat i j))

(defgeneric mzerop (mat &optional threshold)
  (:documentation "Returns t if the matrix is the zero matrix."))
(defgeneric midentity-p (mat &optional threshold)
  (:documentation "Returns t if the matrix is the zero matrix."))
(defgeneric x<-id (mat))
(defgeneric clear-row (mat row &optional row2))
(defgeneric clear-column (mat col &optional col2))
(defgeneric matrix-slice (mat &key from-row from-col nrows ncols))

(defgeneric make-row-vector-for (mat &optional multiplicity))
(defgeneric make-column-vector-for (mat &optional multiplicity))

;;; Operating on separate rows and columns.  For these routines, inlining
;;; may help a lot.

(defgeneric* for-each-row-key (func mat)
  (:documentation "Loop through row keys."))

(defgeneric* for-each-col-key (func mat)
  (:documentation "Loop through column keys."))

(defgeneric* for-each-key-in-row (func mat row-key)
  (:documentation "Loop through col-keys in row."))

(defgeneric* for-each-key-in-col (func mat col-key)
  (:documentation "Loop through row-keys in column col."))

(defgeneric* for-each-entry-in-row (func mat row-key)
  (:documentation "Loop through col-keys in row."))

(defgeneric* for-each-entry-in-col (func mat col-key)
  (:documentation "Loop through entries in column col."))

(defgeneric* for-each-key-and-entry-in-row (func mat row-key)
  (:documentation "Loop through col-keys and entries in row."))

(defgeneric* for-each-key-and-entry-in-col (func mat col-key)
  (:documentation "Loop through row-keys and entries in col."))

(defmacro dorows ((loop-vars mat) &body body)
  (let* ((loop-vars (if (consp loop-vars) loop-vars (list loop-vars)))
	 (mat-for-each-row
	  (if (null (car loop-vars))
	      (if (null (cdr loop-vars))
		  (error "no loop variable")
		  'for-each-row)
	      (if (or (null (cdr loop-vars)) (null (cadr loop-vars)))
		  'for-each-row-key
		  'for-each-key-and-row))))
    `(,mat-for-each-row #'(lambda ,(remove nil loop-vars) ,@body) ,mat)))

(defmacro dorow ((loop-vars row) &body body)
  (let* ((loop-vars (if (consp loop-vars) loop-vars (list nil loop-vars)))
	 (row-for-each (if (null (car loop-vars))
			   (if (null (cdr loop-vars))
			       (error "no loop variable")
			       'for-each-entry-in-row)
			   (if (or (null (cdr loop-vars)) (null (cadr loop-vars)))
			       'for-each-key-in-row
			       'for-each-key-and-entry-in-row))))
    `(,row-for-each #'(lambda ,(remove nil loop-vars) ,@body) ,row)))

(defmethod mzerop (mat &optional (threshold 0.0))
  "Default method."
  (for-each-entry
   #'(lambda (entry)
       (unless (mzerop entry threshold)
	 (return-from mzerop nil)))
   mat))


;;; vector-matrix operations
(defgeneric* x<-Ay (x A y))
(defgeneric* x+=Ay (x A y))
(defgeneric* x-=Ay (x A y))
(defgeneric* x+=s*Ay (x s A y))

;;; Additionally, most reasonable matrix classes will provide some subset
;;; of the vector operations defined in vector.lisp.

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Submatrices
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass <submatrix> (<matrix>)
  ((matrix :initarg :matrix)
   (row-keys :reader row-keys :initarg :row-keys :type vector)
   (col-keys :reader col-keys :initarg :col-keys :type vector))
  (:documentation "Describes an ordered submatrix of a matrix.  Only a
restricted set of operations is allowed for these matrices and element
access is slow.  They are indexed with ordinary integers."))

(defmethod matrix-ref-2d ((submat <submatrix>) i j)
  (matrix-ref (slot-value submat 'matrix)
	      (aref (row-keys submat) i)
	      (aref (col-keys submat) j)))

(defmethod (setf matrix-ref-2d) (value (submat <submatrix>) i j)
  (setf (matrix-ref (slot-value submat 'matrix)
		    (aref (row-keys submat) i)
		    (aref (col-keys submat) j))
	value))

(defmethod nrows ((mat <submatrix>)) (length (row-keys mat)))

(defmethod ncols ((mat <submatrix>)) (length (col-keys mat)))

(defmethod nr-of-entries ((submat <submatrix>))
  (* (nrows submat) (ncols submat)))

(defmethod for-each-key ((fn function) (submat <submatrix>))
  (dotimes (i (nrows submat))
    (dotimes (j (ncols submat))
      (funcall fn i j))))

(defmethod for-each-key-and-entry ((fn function) (submat <submatrix>))
  (dotimes (i (nrows submat))
    (dotimes (j (ncols submat))
      (funcall fn i j (matrix-ref submat i j)))))

(defmethod for-each-entry ((fn function) (submat <submatrix>))
  (dotimes (i (nrows submat))
    (dotimes (j (ncols submat))
      (funcall fn (matrix-ref submat i j)))))


