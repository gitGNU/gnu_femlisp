;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; sparse.lisp
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

;;;; This module defines a graph sparse matrix format.  It is somewhat
;;;; similar to the format used in UG and outlined in [Neuss1998], but it
;;;; uses hash-table sparse matrices which allow for O(1) random access and
;;;; indexing over arbitrary key sets.

;;;(declaim (optimize (safety 3) (debug 3)))

(in-package :fl.algebra)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;
;;; <sparse-vector>
;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass <sparse-vector> ()
  ((blocks :initform (make-hash-table :test #'eq) :type hash-table
	   :documentation "Table of blocks.")
   (key->size :reader key->size :initarg :key->size :type function
	      :documentation "Function determining the dimension of a block.")
   (print-key :reader print-key :initarg :print-key :initform #'princ
	      :documentation "Print function for a key.")
   (multiplicity :reader multiplicity :initform 1 :initarg :multiplicity :type fixnum
		 :documentation "Multiplicity of the sparse vector.  A
multiplicity different from 1 is used when handling multiple right-hand
sides and solutions simultaneously."))
  (:documentation "Sparse block vector class indexed with general keys."))

(defmethod element-type ((svec <sparse-vector>))
  (standard-matrix 'double-float))

(defmethod scalar-type ((svec <sparse-vector>))
  'double-float)

(defmethod nr-of-entries ((svec <sparse-vector>))
  (hash-table-count (slot-value svec 'blocks)))

(defmethod make-analog ((svec <sparse-vector>))
  (make-instance
   '<sparse-vector> :key->size (key->size svec) :print-key (print-key svec)
   :multiplicity (multiplicity svec)))

(defmethod vref ((svec <sparse-vector>) key)
  "Random access to vector components.  Fast version."
  ;;(declare (values real-matrix))
  (let ((blocks (slot-value svec 'blocks)))
    (or (gethash key blocks)
	(setf (gethash key blocks)
	      (make-real-matrix (funcall (key->size svec) key)
				(multiplicity svec))))))

(defmethod (setf vref) (value (svec <sparse-vector>) key)
  "Inserts a vector block into the sparse vector."
  ;;(declare (type real-matrix value))
  (setf (gethash key (slot-value svec 'blocks)) value))

(defmethod for-each-key ((fn function) (svec <sparse-vector>))
  (loop for key being the hash-keys of (slot-value svec 'blocks) do
	(funcall fn key)))
(defmethod for-each-key-and-entry ((fn function) (svec <sparse-vector>))
  (maphash fn (slot-value svec 'blocks)))
(defmethod for-each-entry ((fn function) (svec <sparse-vector>))
  (loop for val being the hash-values of (slot-value svec 'blocks) do
	(funcall fn val)))
(defmethod for-each-entry-of-vec1 ((fn function) (svec1 <sparse-vector>) (svec2 <sparse-vector>))
  (maphash #'(lambda (key val) (funcall fn val (vref svec2 key)))
	   (slot-value svec1 'blocks)))
(defmethod for-each-entry-of-vec2 ((fn function) (svec1 <sparse-vector>) (svec2 <sparse-vector>))
  (maphash #'(lambda (key val) (funcall fn (vref svec1 key) val))
	   (slot-value svec2 'blocks)))

(defmethod keys ((svec <sparse-vector>))
  (let ((keys ()))
    (dovec ((key) svec) (push key keys))
    (nreverse keys)))

(defmethod show ((svec <sparse-vector>) &key keys (zeros t) &allow-other-keys)
  (format t "~&Sparse vector with ~D allocated components:~%"
	  (nr-of-entries svec))
  (flet ((keyout (key vblock)
	   (when (or zeros (not (mzerop vblock)))
	     (funcall (print-key svec) key)
	     (format t " --> ~A~%" vblock))))
    (dolist (key (or keys (keys svec)))
      (keyout key (vref svec key))))
  svec)

(defmethod remove-key ((svec <sparse-vector>) key)
  (remhash key (slot-value svec 'blocks)))

(defmethod remove-keys (sobj keys)
  (loop for key in keys do (remove-key sobj key)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; transformation to matlisp matrices
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmethod sparse-vector->matlisp ((svec <sparse-vector>) &optional keys ranges)
  "Transforms all or a part of a <sparse-vector> corresponding to the keys
in 'keys' and maybe the ranges in 'ranges' to a matlisp matrix."
  (setq keys (coerce (or keys (keys svec)) 'vector))
  (setq ranges (and ranges (coerce ranges 'vector)))
  (let* ((n (if ranges
		(reduce #'+ ranges :key #'(lambda (x) (- (cdr x) (car x))))
		(reduce #'+ keys :key (key->size svec))))
	 (multiplicity (multiplicity svec))
	 (mm (make-real-matrix n multiplicity)))
    (loop for key across keys and k from 0
	  and offset of-type fixnum = 0 then (+ offset (- end-comp start-comp))
	  for start-comp of-type fixnum = (if ranges (car (aref ranges k)) 0)
	  for end-comp of-type fixnum = (if ranges
					    (cdr (aref ranges k))
					    (funcall (key->size svec) key))
	  for entry = (vref svec key) do
	  (loop for i of-type fixnum from start-comp below end-comp do
		(dotimes (j multiplicity)
		  (setf (mref mm (+ offset i) j)
			(if entry (mref entry i j) 0.0)))))
    mm))

(defun combine-svec-block (svec local-vec keys ranges operation)
  "Puts a local block in matlisp format into a <sparse-vector>."
  (setq keys (coerce (or keys (keys svec)) 'vector))
  (setq ranges (and ranges (coerce ranges 'vector)))
  (let ((multiplicity (multiplicity svec)))
    (loop for key across keys and k from 0
	  and offset fixnum = 0 then (+ offset (- end-comp start-comp))
	  for start-comp of-type fixnum = (if ranges (car (aref ranges k)) 0)
	  for end-comp of-type fixnum = (if ranges
					    (cdr (aref ranges k))
					    (funcall (key->size svec) key))
	  for entry = (vref svec key) do
	  (loop for i of-type fixnum from start-comp below end-comp do
		(dotimes (j multiplicity)
		  (let ((local-entry (mref local-vec (+ offset i) j)))
		    (ecase operation
		      (:add (incf (mref entry i j) local-entry))
		      (:set (setf (mref entry i j) local-entry)))))))))
  
(defmethod set-svec-to-local-block ((svec <sparse-vector>) local-vec &optional keys ranges)
  "Copies a local block in matlisp format into a <sparse-vector>."
  (combine-svec-block svec local-vec keys ranges :set))

(defmethod add-svec-to-local-block ((svec <sparse-vector>) local-vec &optional keys ranges)
  "Copies a local block in matlisp format into a <sparse-vector>."
  (combine-svec-block svec local-vec keys ranges :add))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; matlisp operations for the <sparse-vector> class
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmethod copy ((x <sparse-vector>))
  (let ((result (make-analog x)))
    (copy! x result)
    result))

(defmethod extract-if ((test function) (svec <sparse-vector>) &key &allow-other-keys)
  "Extracts a sub-matrix from a sparse matrix."
  (let ((sub-vec (make-analog svec)))
    (dovec ((key entry) svec)
      (when (funcall test key entry)
	(setf (vref sub-vec key) entry)))
    sub-vec))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; vector blas operations for the <svec> class
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; We use the general operations for block vectors.  Later on these will
;;; be inlined.

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;
;;; <sparse-matrix>
;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass <sparse-matrix> (<matrix>)
  ((row-table :accessor row-table :initarg :row-table
	      :initform (make-hash-table :test #'eq) :type hash-table
	      :documentation "Table of rows.")
   (column-table :accessor column-table :initarg :column-table
		 :initform (make-hash-table :test #'eq) :type hash-table
		 :documentation "Table of columns.")
   (print-row-key :reader print-row-key :initarg :print-row-key
		  :initform #'princ :type function
		  :documentation "Print function for row keys.")
   (print-col-key :reader print-col-key :initarg :print-col-key
		  :initform #'princ :type function
		  :documentation "Print function for column keys.")
   (row-key->size :reader row-key->size :initarg :row-key->size
		  :type (function (t) positive-fixnum) :documentation
		  "Function determining the number of rows of a block.")
   (col-key->size :reader col-key->size :initarg :col-key->size
		  :type (function (t) positive-fixnum) :documentation
		  "Function determining the number of columns of a block.")
   (keys->pattern :reader keys->pattern :initarg :keys->pattern :type function :documentation
		  "Function determining the pattern of a block."))
  (:documentation "The <sparse-matrix> represents an unordered matrix graph
indexed by general keys."))

(defmethod element-type ((smat <sparse-matrix>))
  (standard-matrix 'double-float))

(defmethod scalar-type ((smat <sparse-matrix>))
  'double-float)

(defun make-sparse-matrix (&key (print-row-key #'princ) (print-col-key #'princ)
			   row-key->size col-key->size keys->pattern)
  (make-instance
   '<sparse-matrix>
   :print-row-key print-row-key :print-col-key print-col-key :keys->pattern keys->pattern
   :row-key->size row-key->size :col-key->size col-key->size))

(defmethod nr-of-entries ((smat <sparse-matrix>))
  (loop for entry being the hash-values of (row-table smat)
	summing (hash-table-count entry)))

(defmethod nrows ((mat <sparse-matrix>))
  (hash-table-count (row-table mat)))

(defmethod ncols ((mat <sparse-matrix>))
  (hash-table-count (column-table mat)))

(defmethod total-nrows ((mat <sparse-matrix>))
  (let ((sum 0))
    (for-each-row-key
     #'(lambda (row-key)
	 (incf sum (funcall (row-key->size mat) row-key)))
     mat)
    sum))

(defmethod total-entries ((mat <sparse-matrix>))
  (let ((entries 0))
    (for-each-entry
     #'(lambda (entry)
	 (unless entry (break))
	 (incf entries (* (nrows entry) (ncols entry))))
     mat)
    entries))

(defmethod nr-nonempty-rows (mat) (nrows mat))
(defmethod nr-nonempty-rows ((smat <sparse-matrix>))
  (hash-table-count (row-table smat)))

(defmethod nr-nonempty-columns (mat) (ncols mat))
(defmethod nr-nonempty-columns ((smat <sparse-matrix>))
  (hash-table-count (column-table smat)))

(defmethod mzerop ((smat <sparse-matrix>) &optional (threshold 0.0))
  (declare (ignore threshold))
  (or (zerop (nr-nonempty-rows smat))
      (call-next-method)))

(defmethod make-analog ((smat <sparse-matrix>))
  (make-sparse-matrix
   :print-row-key (print-row-key smat)
   :print-col-key (print-col-key smat)
   :keys->pattern (keys->pattern smat)
   :row-key->size (row-key->size smat)
   :col-key->size (col-key->size smat)))

(defun make-sparse-automorphism (&key (print-key #'princ) key->size keys->pattern)
  (make-instance
   '<sparse-matrix> :print-row-key print-key :print-col-key print-key
   :keys->pattern keys->pattern :row-key->size key->size :col-key->size key->size))

(defun automorphism? (A)
  (and (eql (print-row-key A)
	    (print-col-key A))
       (eql (row-key->size A)
	    (col-key->size A))))

(defmethod matrix-row ((mat <sparse-matrix>) key)
  (values (gethash key (row-table mat))))

(defmethod (setf matrix-row) ((row-table hash-table) (mat <sparse-matrix>) key)
  (setf (gethash key (row-table mat)) row-table))

(defmethod matrix-column ((mat <sparse-matrix>) key)
  (values (gethash key (column-table mat))))

(defmethod (setf matrix-column) ((col-table hash-table) (mat <sparse-matrix>) key)
  (setf (gethash key (column-table mat)) col-table))

(defmethod keys-of-row ((mat <sparse-matrix>) key)
  (aand (matrix-row mat key) (hash-table-keys it)))
(defmethod keys-of-column ((mat <sparse-matrix>) key)
  (aand (matrix-column mat key) (hash-table-keys it)))

(defun make-full-block-analog (sm)
  (declare (type <sparse-matrix> sm))
  (make-sparse-matrix
   :print-row-key (print-row-key sm)
   :print-col-key (print-col-key sm)
   :row-key->size (row-key->size sm)
   :col-key->size (col-key->size sm)
   :keys->pattern
   #'(lambda (row-key col-key)
       (full-crs-pattern (funcall (row-key->size sm) row-key)
			 (funcall (row-key->size sm) col-key)))))

(defmethod matrix-block ((smat <sparse-matrix>) row-key col-key)
  (aand (matrix-row smat row-key)
	(gethash col-key it)))

(defmethod mref ((smat <sparse-matrix>) row-key col-key)
  "If the matrix-block indexed by row-key/col-key exists, it is returned.
Otherwise, a new matrix block is created according to the pattern specified
for this index."
  (or (matrix-block smat row-key col-key)
      (setf (mref smat row-key col-key)
	    (let ((pattern (funcall (keys->pattern smat) row-key col-key)))
	      (with-slots (nrows ncols store-size) pattern
		(cond
		  ((or (zerop nrows) (zerop ncols))
		   (error "No sparse matrix entry allowed here!"))
		  ((= store-size (* nrows ncols))
		   (make-real-matrix nrows ncols))
		  (t (make-instance 'crs-matrix :pattern pattern
				    :store (make-double-vec store-size)))))))))

(defmethod (setf mref) (value (smat <sparse-matrix>) row-key col-key)
  "Inserts a matrix block into the sparse matrix.  Note: using this routine
means a lot of consing."
  (setf (gethash col-key
		 (or (gethash row-key (row-table smat))
		     (setf (gethash row-key (row-table smat))
			   (make-hash-table :size 30 :test #'eq))))
	value)
  (setf (gethash row-key
		 (or (gethash col-key (column-table smat))
		     (setf (gethash col-key (column-table smat))
			   (make-hash-table :size 30 :test #'eq))))
	value))

(defmethod remove-entry ((smat <sparse-matrix>) row-key col-key)
  (let ((row (matrix-row smat row-key))
	(column (matrix-column smat col-key)))
    (when row
      (remhash col-key row)
      (when (zerop (hash-table-count row))
	(remhash row-key (row-table smat)))
      (remhash row-key column)
      (when (zerop (hash-table-count column))
	(remhash col-key (column-table smat))))))

(defmethod remove-row ((smat <sparse-matrix>) row-key)
  "Warning: destructive operation."
  (awhen (matrix-row smat row-key)
    (loop for col-key being the hash-keys of it do
	  (remhash row-key (matrix-column smat col-key))
	  finally
	  (remhash row-key (row-table smat)))))

(defmethod remove-column ((smat <sparse-matrix>) col-key)
  "Warning: destructive operation."
  (awhen (matrix-column smat col-key)
    (loop for row-key being the hash-keys of it do
	  (remhash col-key (matrix-row smat row-key))
	  finally
	  (remhash col-key (column-table smat)))))

(defmethod remove-key ((smat <sparse-matrix>) key)
  (remove-row smat key)
  (remove-column smat key))

(defmethod row<-id ((A <sparse-matrix>) key)
  (remove-row A key)
  (setf (mref A key key)
	(eye (funcall (row-key->size A) key)
	     (funcall (col-key->size A) key))))

(defmethod column<-id ((A <sparse-matrix>) key)
  (remove-column A key)
  (setf (mref A key key)
	(eye (funcall (row-key->size A) key)
	     (funcall (col-key->size A) key))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Vector operation support
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmethod for-each-key ((fn function) (smat <sparse-matrix>))
  (maphash
   #'(lambda (row-key row-dic)
       (loop for col-key being each hash-key of row-dic do
	     (funcall fn (cons row-key col-key))))
   (row-table smat)))

(defmethod for-each-entry ((fn function) (smat <sparse-matrix>))
  (loop for row-dic being the hash-values of (row-table smat) do
	(loop for entry being the hash-values of row-dic do
	      (funcall fn entry))))

(defmethod for-each-key-and-entry ((fn function) (smat <sparse-matrix>))
  (maphash #'(lambda (row-key row-dic)
	       (maphash #'(lambda (col-key entry)
			    (funcall fn (cons row-key col-key) entry))
			row-dic))
	   (row-table smat)))

(defmethod for-each-row-key ((fn function) (smat <sparse-matrix>))
  "Loop through row keys."
  (loop for key being the hash-keys of (row-table smat) do
	(funcall fn key)))
(defmethod for-each-col-key ((fn function) (smat <sparse-matrix>))
  "Loop through column keys."
  (loop for key being the hash-keys of (column-table smat) do
	(funcall fn key)))

;;; Operating on separate rows and columns.  For these routines, inlining
;;; may help a lot.
(defmethod for-each-entry-in-row ((fn function) (smat <sparse-matrix>) key)
  (awhen (matrix-row smat key)
    (loop for mblock being the hash-values of it
	  do (funcall fn mblock))))
(defmethod for-each-entry-in-col ((fn function) (smat <sparse-matrix>) key)
  (awhen (matrix-column smat key)
    (loop for mblock being the hash-values of it
	  do (funcall fn mblock))))

(defmethod for-each-key-in-row ((fn function) (smat <sparse-matrix>) key)
  (whereas ((row (matrix-row smat key)))
    (loop for col-key being the hash-keys of row
	  do (funcall fn col-key))))
(defmethod for-each-key-in-col ((fn function) (smat <sparse-matrix>) key)
  (whereas ((col (matrix-column smat key)))
    (loop for row-key being the hash-keys of col
	  do (funcall fn row-key))))

(defmethod for-each-key-and-entry-in-row ((fn function) (smat <sparse-matrix>) key)
  (awhen (matrix-row smat key)
    (loop for col-key being the hash-keys of it
	  and mblock being the hash-values of it
	  do (funcall fn col-key mblock))))
(defmethod for-each-key-and-entry-in-col ((fn function) (smat <sparse-matrix>) key)
  (awhen (matrix-column smat key)
    (loop for row-key being the hash-keys of it
	  and mblock being the hash-values of it
	  do (funcall fn row-key mblock))))

(defmethod clear-row ((mat <sparse-matrix>) row-key1 &optional row-key2)
  (for-each-entry-in-row
   #'(lambda (mblock)
       (if row-key2
	   (clear-row mblock row-key2)
	   (x<-0 mblock)))
   mat row-key1))

(defmethod clear-column ((mat <sparse-matrix>) col-key1 &optional col-key2)
  (for-each-entry-in-col
   #'(lambda (mblock)
       (if col-key2
	   (clear-column mblock col-key2)
	   (x<-0 mblock)))
   mat col-key1))

(defmethod row-keys ((mat <sparse-matrix>))
  (let ((row-keys ()))
    (for-each-row-key #'(lambda (key) (push key row-keys)) mat)
    (nreverse row-keys)))

(defmethod col-keys ((mat <sparse-matrix>))
  (let ((col-keys ()))
    (for-each-col-key #'(lambda (key) (push key col-keys)) mat)
    (nreverse col-keys)))

(defmethod symmetric-p ((smat <sparse-matrix>) &key (threshold 0.0) output)
  (let ((flag t))
    (for-each-key-and-entry
     #'(lambda (keys entry)
	 (let ((entry2 (mref smat (cdr keys) (car keys))))
	   (when (> (norm (m- entry2 (transpose entry))) threshold)
	     (when output
	       (format t "~&Mismatch~%(i,j)=~A:~%Aij=~A~%Aji=~A~%~%"
		       keys entry entry2))
	     (setq flag nil))))
     smat)
    flag))

(defmethod matrix-transpose-instance ((smat <sparse-matrix>))
  (make-sparse-matrix
   :print-row-key (print-col-key smat)
   :print-col-key (print-row-key smat)
   :keys->pattern #'(lambda (rk ck) (funcall (keys->pattern smat) ck rk))
   :row-key->size (col-key->size smat)
   :col-key->size (row-key->size smat)))


(defmethod transpose! ((x <sparse-matrix>) (y <sparse-matrix>))
  (for-each-key-and-entry
     #'(lambda (keys entry)
	 (setf (mref y (cdr keys) (car keys)) (transpose entry)))
     x)
  y)

(defmethod show ((smat <sparse-matrix>) &key keys (zeros t) &allow-other-keys)
  (assert (null keys))
  (format t "~&Sparse matrix:~%")
  (for-each-row-key
   #'(lambda (row-key)
       (format t " ~%*** ")
       (funcall (print-row-key smat) row-key)
       (format t " ***~%")
       (maphash
	#'(lambda (col-key vblock)
	    (when (or zeros (not (mzerop vblock)))
	      (funcall (print-col-key smat) col-key)
	      (format t " -> ~A~%" vblock)))
	(matrix-row smat row-key)))
   smat)
  smat)

(defmethod display ((smat <sparse-matrix>) &key row-order col-order order)
  (unless row-order
    (setq row-order (or order (row-keys smat))))
  (let ((row-indices (let ((result (make-hash-table)))
		       (loop for i from 0 and key in row-order do
			     (setf (gethash key result) i)
			     finally (return result)))))
    (unless col-order (setq col-order order))
    (unless col-order
      (setq col-order (col-keys smat))
      (when (some (rcurry #'gethash row-indices) col-order)
	(assert (every (rcurry #'gethash row-indices) col-order))
	(setq col-order row-order)))
    (flet ((inter-line (key1 key2)
	     (format t "~&+")
	     (loop for col-key in col-order do
		   (format t (if (or (eql col-key key1) (eql col-key key2))
				 "~V,,,'=<~>+"
				 "~V,,,'-<~>+")
			   (* 10 (funcall (col-key->size smat) col-key))))
	     (format t "~%"))
	   (inter-mark (row-key col-key1 col-key2)
	     (format t (if (or (eql row-key col-key1) (eql row-key col-key2))
			   "I" "|"))))
      (loop
       initially (inter-line (car row-order) (car row-order))
       for (row-key . row-keys) on row-order and row from 0 do
       (dotimes (i (funcall (row-key->size smat) row-key))
	 (inter-mark row-key (car col-order) (car col-order))
	 (loop for (col-key . col-keys) on col-order and col from 0
	       for entry = (matrix-block smat row-key col-key)
	       for ncols-block = (funcall (col-key->size smat) col-key) do
	       (if entry
		   (dotimes (j ncols-block)
		     (format t "~9,2,,,,,'Eg " (mref entry i j)))
		   (format t "~V,,,' <~>" (* ncols-block 10)))
	       (inter-mark row-key col-key (car col-keys)))
	 (format t "~%"))
       (inter-line row-key (car row-keys))))))

(defmethod mat-diff ((smat1 <sparse-matrix>) (smat2 <sparse-matrix>))
  (format t "Missing in [2]~%")
  (for-each-key-and-entry
   #'(lambda (keys entry)
       (unless (matrix-block smat2 (car keys) (cdr keys))
	 (format t "~A:~%~A~%" keys entry)))
   smat1)
  (format t "Missing in [1]~%")
  (for-each-key-and-entry
   #'(lambda (keys entry)
       (unless (matrix-block smat1 (car keys) (cdr keys))
	 (format t "~A:~%~A~%" keys entry)))
   smat2)
  (format t "Differences~%")
  (for-each-key-and-entry
   #'(lambda (keys entry)
       (whereas ((entry2 (matrix-block smat2 (car keys) (cdr keys))))
	 (unless (mzerop (m- entry entry2))
	   (format t "~A :~%" keys)
	   (format t "~A~%~A~%" entry entry2))))
   smat1))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Transformation of unknowns
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; This is needed by hanging nodes, identified boundaries, ...

(defun index-range-disjoint-p (mat1 mat2)
  "Checks if the range of indices of two sparse matrices is disjoint."
  (for-each-row-key #'(lambda (row)
	 		(when (matrix-row mat2 row)
			  (return-from index-range-disjoint-p nil)))
 		    mat1)
  t)

(defun range-and-domain-disjoint-p (mat)
  "Checks if index range and index domain of some matrix are disjoint."
  (for-each-row-key #'(lambda (row)
	 		(when (matrix-column mat row)
			  (return-from range-and-domain-disjoint-p nil)))
 		    mat)
  t)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Matrix-vector stuff
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


(defmethod make-image-vector-for ((A <sparse-matrix>) &optional (multiplicity 1))
  (let ((result (make-instance '<sparse-vector> :key->size (row-key->size A)
			       :print-key (print-row-key A) :multiplicity multiplicity)))
    (for-each-row-key (curry #'vref result) A)
    result))

(defmethod make-domain-vector-for ((A <sparse-matrix>) &optional (multiplicity 1))
  (let ((result (make-instance '<sparse-vector> :key->size (col-key->size A)
			       :print-key (print-col-key A) :multiplicity multiplicity)))
    (for-each-col-key (curry #'vref result) A)
    result))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; vector blas operations for the <smat> class
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; those are defined via the general methods in vector.lisp

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; matrix-vector blas operations <smat>/<svec>
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(definline sparse-vecmat-loop (operation x A y)
  (maphash
   #'(lambda (row-key row-dic)
       (let ((x-values (vref x row-key)))
	 (maphash
	  #'(lambda (col-key mblock)
	      (let ((y-values (vref y col-key)))
		(funcall operation x-values mblock y-values)))
	  row-dic)))
   (row-table A)))

(defmethod x<-Ay ((x <sparse-vector>) (A <sparse-matrix>) (y <sparse-vector>))
  (x-on-range-of-A<-0 x A)
  (sparse-vecmat-loop #'x+=Ay x A y))
(defmethod x+=Ay ((x <sparse-vector>) (A <sparse-matrix>) (y <sparse-vector>))
  (sparse-vecmat-loop #'x+=Ay x A y))
(defmethod x-=Ay ((x <sparse-vector>) (A <sparse-matrix>) (y <sparse-vector>))
  (sparse-vecmat-loop #'x-=Ay x A y))

(defmethod x-on-range-of-A<-0 ((x <sparse-vector>) (A <sparse-matrix>))
  (for-each-row-key #'(lambda (key) (x<-0 (vref x key))) A))

(defmethod x-on-domain-of-A<-0 ((x <sparse-vector>) (A <sparse-matrix>))
  (for-each-col-key #'(lambda (key) (x<-0 (vref x key))) A))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; matlisp operations for the <smat> class
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmethod copy ((x <sparse-matrix>))
  (let ((result (make-analog x)))
    (copy! x result)
    result))

(defmethod m*-product-instance ((A <sparse-matrix>) (y <sparse-vector>))
  (make-image-vector-for A (multiplicity y)))

(defun generate-sparse-matrix-vector-gemm!-template (job)
  "Generates the GEMM-XX! routine defined by JOB."
  (assert (member job '(:nn :nt :tn :tt)))
  (let ((gemm-job (symconc "GEMM-" (symbol-name job) "!"))
	(table-accessor (if (member job '(:nn :nt))
			    'row-table
			    'column-table)))
    (eval
     `(defmethod ,gemm-job
       (alpha (A <sparse-matrix>) (y <sparse-vector>) beta (x <sparse-vector>))
       (declare (optimize (speed 3) (space 2) (safety 1)))
       (maphash
	#'(lambda (x-key row-or-col-dic)
	    (let ((x-values (vref x x-key)))
	      (scal! beta x-values)
	      (maphash
	       #'(lambda (y-key mblock)
		   (let ((y-values (vref y y-key)))
		     (,gemm-job alpha mblock y-values 1.0 x-values)))
	       row-or-col-dic)))
	(,table-accessor A))
       x))))

(mapc #'generate-sparse-matrix-vector-gemm!-template '(:nn :nt :tn :tt))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; sparse matrix multiplication
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defgeneric sparse-m* (A B &key job sparsity)
  (:documentation "Sparse matrix-matrix or matrix-vector multiplication.
Usually, m* should be used.  But in situations, where A or B are very
sparse, the complexity of this routine is much lower."))

(defmethod sparse-m* ((A <sparse-matrix>) (B <sparse-vector>)
		      &key (job :nn) (sparsity :A))
  (let ((result (ecase job
		  ((:nn :tn) (make-image-vector-for A (multiplicity B)))
		  ((:nt :tt) (make-domain-vector-for A (multiplicity B))))))
    (ecase sparsity
      (:A (maphash
	   #'(lambda (i row-or-column)
	       (maphash
		#'(lambda (j Aij)
		    (gemm! 1.0 Aij (vref B j) 1.0 (vref result i) job))
		row-or-column))
	   (ecase job
	     ((:nn :nt) (row-table A))
	     ((:tn :tt) (column-table A)))))
      (:B (for-each-key-and-entry
	   #'(lambda (i bi)
	       (whereas ((row-or-column (ecase job
					  ((:nn :nt) (matrix-column A i))
					  ((:tn :tt) (matrix-row A i)))))
		 (maphash #'(lambda (j Aij) (gemm! 1.0 Aij bi 1.0 (vref result j) job))
			  row-or-column)))
	   B)))
    result))


(defmethod sparse-m* ((A <sparse-matrix>) (B <sparse-matrix>)
		      &key (job :nn) (sparsity :A))
  (let* ((transpose-A? (or (eq job :tn) (eq job :tt)))
	 (transpose-B? (or (eq job :nt) (eq job :tt)))
	 (row-key->size (if transpose-A? (col-key->size A) (row-key->size A)))
	 (col-key->size (if transpose-B? (row-key->size B) (col-key->size B)))
	 (C (make-sparse-matrix
	     :print-row-key (if transpose-A? (print-col-key A) (print-row-key A))
	     :print-col-key (if transpose-B? (print-row-key B) (print-col-key B))
	     :row-key->size row-key->size :col-key->size col-key->size
	     :keys->pattern
	     #'(lambda (row-key col-key) ; at the moment this works only for full blocks
		 (full-crs-pattern (funcall row-key->size row-key)
				   (funcall col-key->size col-key))))))
  (ecase sparsity
    ((:A)
     (maphash
      #'(lambda (i Ai)
	  (maphash
	   #'(lambda (j Aij)
	       (whereas ((row-or-column (ecase job
					  ((:nn :tn) (matrix-row B j))
					  ((:nt :tt) (matrix-column B j)))))
		 (maphash #'(lambda (k Bjk)
			      (gemm! 1.0 Aij Bjk 1.0 (mref C i k) job))
			  row-or-column)))
	   Ai))
      (ecase job
	((:nn :nt) (row-table A))
	((:tn :tt) (column-table A)))))
    ((:B)
     (maphash
      #'(lambda (k Bk)
	  (maphash
	   #'(lambda (j Bjk)
	       (whereas ((row-or-column (ecase job
					  ((:nn :nt) (matrix-column A j))
					  ((:tn :tt) (matrix-row A j)))))
		 (maphash #'(lambda (i Aij)
			      (gemm! 1.0 Aij Bjk 1.0 (mref C i k) job))
			  row-or-column)))
	   Bk))
      (ecase job
	((:nn :tn) (column-table B))
	((:nt :tt) (row-table B))))))
  C))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; submatrices / extraction
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; Remark: Several operations can only be performed when an ordering is
;;; specified.  This must augment the sparse-matrix structure which is only
;;; an encoding for the linear mapping itself - on the other hand, maybe we
;;; shouldn't call it matrix then...

(defgeneric extended-extract (mat keys &key row? col?)
  (:documentation "Extracts a sub-matrix from a sparse matrix."))

(defmethod extract-if ((test function) (smat <sparse-matrix>) &key &allow-other-keys)
  "Extracts a sub-matrix from a sparse matrix."
  (let ((sub-mat (make-analog smat)))
    (for-each-key-and-entry
     #'(lambda (keys entry)
	 (when (funcall test (car keys) (cdr keys) entry)
	   (setf (vref sub-mat keys) entry)))
     smat)
    sub-mat))

(defmethod extended-extract ((svec <sparse-vector>) (keys hash-table) &key &allow-other-keys)
  (extract-if #'(lambda (key entry)
		  (declare (ignore entry))
		  (gethash key keys))
	      svec))

(defmethod extended-extract ((smat <sparse-matrix>) (keys hash-table)
			     &key (row? t) (col? t))
  "Extracts a sub-matrix from a sparse matrix.  This routine could be
accelerated by taking member-checks out of column or row loop."
  (extract-if #'(lambda (row-key col-key entry)
		  (declare (ignore entry))
		  (and (or (not row?) (gethash row-key keys))
		       (or (not col?) (gethash col-key keys))))
	      smat))


(defmethod extract-matrix-block ((smat <sparse-matrix>) row-keys col-keys)
  "Extracts a sub-matrix from a sparse matrix.  row-keys=nil or
col-keys=nil means to allow every key."
  (extract-if #'(lambda (row-key col-key entry)
		  (declare (ignore entry))
		  (or (and row-keys (not (gethash row-key row-keys)))
		      (and col-keys (not (gethash col-key col-keys)))))
	      smat))

(defmethod submatrix ((smat <sparse-matrix>) &key row-indices col-indices)
  (extract-matrix-block smat row-indices col-indices))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; transformation to matlisp matrices
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmethod sparse-matrix->matlisp ((A <sparse-matrix>)
				   &key keys row-keys col-keys
				   ranges row-ranges col-ranges)
  (setq row-keys (coerce (or row-keys keys (row-keys A)) 'vector))
  (setq col-keys (coerce (or col-keys keys (col-keys A)) 'vector))
  (setq row-ranges (aand (or row-ranges ranges) (coerce it 'vector)))
  (setq col-ranges (aand (or col-ranges ranges) (coerce it 'vector)))
  (let* ((n (if row-ranges
		(reduce #'+ row-ranges :key #'(lambda (x) (- (cdr x) (car x))))
		(reduce #'+ row-keys :key (row-key->size A))))
	 (m (if col-ranges
		(reduce #'+ col-ranges :key #'(lambda (x) (- (cdr x) (car x))))
		(reduce #'+ col-keys :key (col-key->size A))))
	 (mm (make-real-matrix n m)))
    ;; copy the matrix values
    (loop for row-key across row-keys and k from 0
	  and row-offset of-type fixnum = 0 then (+ row-offset (- row-b row-a))
	  for row-a of-type fixnum = (if row-ranges (car (aref row-ranges k)) 0)
	  for row-b of-type fixnum = (if row-ranges
					 (cdr (aref row-ranges k))
					 (funcall (row-key->size A) row-key))
	  do
	  (loop for col-key across col-keys and l from 0
		and col-offset of-type fixnum = 0 then (+ col-offset (- col-b col-a))
		for col-a of-type fixnum = (if col-ranges (car (aref col-ranges l)) 0)
		for col-b of-type fixnum = (if col-ranges
					       (cdr (aref col-ranges l))
					       (funcall (col-key->size A) col-key))
		and mblock = (matrix-block A row-key col-key) do
		(loop for i of-type fixnum from row-a below row-b do
		      (loop for j of-type fixnum from col-a below col-b do
			    (setf (mref mm (+ row-offset i) (+ col-offset j))
				  (if mblock
				      (mref mblock i j)
				      0.0))))))
    mm))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Local matrix manipulation
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun combined-projection (P1 P2)
  "Returns a projection to the range of the given projections."
  (let ((sum (m+ P1 P2)))
    (m- sum (sparse-m* P1 P2))))


(defmethod remove-projection-range ((vec <sparse-vector>) projection &key &allow-other-keys)
  (for-each-row-key
   #'(lambda (i)
       (let* ((P (mref projection i i))
	      (S (m- (eye (nrows P)) P)))
	 (setf (vref vec i)
	       (m* S (vref vec i)))))
   projection)
  vec)

(defmethod remove-projection-range ((mat <sparse-matrix>) projection &key row-p column-p)
  (for-each-row-key
   #'(lambda (i)
       (let* ((P (mref projection i i))
	      (S (m- (eye (nrows P)) P)))
	 (when column-p
	   (for-each-key-in-col
	    #'(lambda (j)
		(setf (mref mat j i)
		      (m* (mref mat j i) S)))
	    mat i))
	 (when row-p
	   (for-each-key-in-row
	    #'(lambda (j)
		(setf (mref mat i j)
		      (m* S (mref mat i j))))
	    mat i))))
   projection)
  mat)

(defun extend-by-identity (mat extend &key ignore (copy t))
  "Extends A such that the keys in extend which are not in ignore are
mapped to identity."
  (when copy (setq mat (copy mat)))
  (loop for key being each hash-key of extend
	unless (and ignore (gethash key ignore))
	do
	(assert (not (matrix-row mat key)))
	(row<-id mat key)
	finally (return mat)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Special matrices
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun laplace-sparse-matrix (n)
  "Generates a sparse matrix for a 1-dimensional Laplace problem
discretized with the 3-point stencil on a structured mesh."
  (let* ((A (make-sparse-automorphism
	     :key->size (constantly 1)
	     :keys->pattern (constantly (full-crs-pattern 1 1)))))
    (dotimes (i n)
      (setf (mref A i i) #m(2.0))
      (when (> i 0) (setf (mref A i (1- i)) #m(-1.0)))
      (when (< i (1- n)) (setf (mref A i (1+ i)) #m(-1.0))))
    A))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Tests
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun test-sparse ()
  "Test the hash-table based sparse-matrix implementation."
  (flet ((constantly-1 (x) (declare (ignore x)) 1)
	 (constantly-2 (x) (declare (ignore x)) 2))
    (let* ((x (make-instance '<sparse-vector> :key->size #'constantly-1))
	   (y (make-instance '<sparse-vector> :key->size #'constantly-1))
	   (A (make-sparse-matrix
	       :row-key->size #'constantly-1 :col-key->size #'constantly-1
	       :keys->pattern (constantly (full-crs-pattern 1 1))))
	   (B (make-analog A))
	   (AA (make-sparse-matrix
		:row-key->size #'constantly-2 :col-key->size #'constantly-2
		:keys->pattern (constantly (full-crs-pattern 2 2)))))
      
      (setf (mref AA 0 1) #m((1.0 2.0e-15) (3.0 -4.0)))
      (display AA :order '(0 1))
      (assert (mzerop A))
      (setf (vref x 1) #m((1.0)))
      (setf (vref y 1) #m((1.0)))
      (copy! y x)
      (assert (= (mref (vref x 1) 0 0) (mref (vref y 1) 0 0)))
      (x<-0 x)
      (assert (mzerop (vref x 1)))
      (axpy! 2.0 y x)
      (assert (= (norm x 1) 2.0))
      ;;
      (setf (mref B 0 0) #m((1.0)))
      (terpri)
      (show B)
      (show x)
      (show (m* B x))
  
      (x<-0 A)
      (show B)
      (show (copy! B A))
      (assert (= (vref (mref A 0 0) 0) 1.0))
      ;;
      (sparse-matrix->matlisp A :keys '(0 1 2))
      (setf (mref AA 0 1) (make-full-crs-matrix 2 2))
      (fill! (mref AA 0 1) 1.0)
      (sparse-matrix->matlisp AA :row-keys '(0 1) :col-keys '(1))
      ;;
      (setf (mref A 0 1) #m((2.0)))
      (show (sparse-m* A A))
      (show (sparse-m* A A :job :nt :sparsity :B))
      (show (sparse-m* A A :job :tn :sparsity :A))
      (show (sparse-m* A A :job :tt))
      ;; end of testing environment
      ))
  ;; end of test procedure
  )

;;; (test-sparse)
(fl.tests:adjoin-test 'test-sparse)

