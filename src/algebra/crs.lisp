;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; crs.lisp
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

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; <crs-pattern> : sparse matrix pattern
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; We use a compact row-ordered scheme (CRS) with identification, see
;;; [Neuss1998] for more details.
(defstruct (<crs-pattern> (:conc-name crs-))
  "<crs-pattern>: This class defines a sparse pattern for use within <crs-matrix>."
  (nrows 0 :type positive-fixnum)
  (ncols 0 :type positive-fixnum)
  (nr-of-entries 0 :type positive-fixnum)
  (store-size 0 :type positive-fixnum)
  (row-starts (required-argument) :type fixnum-vec)
  (col-inds (required-argument) :type fixnum-vec)
  (offsets (required-argument) :type fixnum-vec))

(defun make-crs-pattern (nrows ncols pattern)
  "make-crs-pattern: This is the crs-pattern constructor.  A sparse matrix of the form
  | * 0 0 0 a |
  | 0 a 0 0 0 |
can be described by its dimensions nrows=2, ncols=5 together with the pattern
 '( ((* . 0) (a . 4))  ((a . 1)) )
* means a non-identified value.  Other symbols can be used to identify entries."
  (when (some #'(lambda (lst) (not (apply #'< (mapcar #'cdr lst))))
	      pattern)
    (error "unsorted patterns are not allowed"))
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
    (make-<crs-pattern>
     :nrows nrows :ncols ncols
     :nr-of-entries (length offsets)
     :store-size current-offset
     :row-starts (list->fixnum-vec
		  (cons 0 (loop for row in pattern
				sum (length row) into rs
				collect rs)))
     :offsets (list->fixnum-vec offsets)
     :col-inds (map 'fixnum-vec #'cdr flattened-pattern))))

(defun full-crs-pattern (nrows ncols)
  "Returns trivial rectangular crs-patterns."
  (let ((N (* nrows ncols)))
    (make-<crs-pattern>
     :nrows nrows :ncols ncols
     :nr-of-entries N
     :store-size N
     :row-starts (list->fixnum-vec
		  (loop for i from 0 upto nrows
			collect (* i ncols)))
     :offsets (list->fixnum-vec (range 0 (1- N)))
     :col-inds (list->fixnum-vec
		(loop for i from 0 below N
		      collect (mod i ncols))))))

(defun pattern->full-pattern (pattern)
  (full-crs-pattern (nrows pattern) (ncols pattern)))

(defgeneric shift-pattern (pattern shift)
  (:documentation "shift-pattern: This function shifts a pattern to
its actual offsets in the sparse graph."))

(defmethod shift-pattern ((pattern <crs-pattern>) shift)
  (make-<crs-pattern>
   :nrows (crs-nrows pattern)
   :ncols (crs-ncols pattern)
   :nr-of-entries (crs-nr-of-entries pattern)
   :store-size (crs-store-size pattern)
   :row-starts (crs-row-starts pattern)
   :col-inds (crs-col-inds pattern)
   :offsets (map 'fixnum-vec #'(lambda (offset) (+ offset shift))
		 (crs-offsets pattern))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; <crs-matrix>
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defstruct (<crs-matrix> (:conc-name crs-))
  "<crs-matrix>: This class is a combination of <crs-pattern> and a value vector."
  (pattern (required-argument) :type <crs-pattern>)
  (store (required-argument) :type double-vec))

(defmethod nrows ((A <crs-matrix>))
  (crs-nrows (crs-pattern A)))

(defmethod ncols ((A <crs-matrix>))
  (crs-ncols (crs-pattern A)))

(defmethod nr-of-entries ((A <crs-matrix>))
  (crs-nr-of-entries (crs-pattern A)))

(defun make-crs-matrix (pattern store)
  "make-crs-matrix: <crs-matrix> constructor."
  (if (not (= (crs-store-size pattern) (length store)))
      (error "pattern does not fit with value vector"))
  (make-<crs-matrix>
   :pattern pattern
   :store store))

(defun make-full-crs-matrix (nrows ncols)
  (make-crs-matrix (full-crs-pattern nrows ncols)
		   (make-double-vec (* nrows ncols))))

(defmethod* mat-ref ((A <crs-matrix>) (i fixnum) (j fixnum))
  (declare (values double-float))
  (let* ((pattern (crs-pattern A))
	 (row-starts (crs-row-starts pattern))
	 (col-inds (crs-col-inds pattern)))
    (loop for k from (aref row-starts i) below (aref row-starts (1+ i))
	  do (if (= (aref col-inds k) j)
		 (return (aref (crs-store A) (aref (crs-offsets pattern) k))))
	  finally (return 0.0d0))))

(defmethod* (setf mat-ref) ((val double-float) (A <crs-matrix>) (i fixnum) (j fixnum))
  (declare (values double-float))
  (let* ((pattern (crs-pattern A))
	 (row-starts (crs-row-starts pattern))
	 (col-inds (crs-col-inds pattern)))
    (loop for k from (aref row-starts i) below (aref row-starts (1+ i))
	  do (if (= (aref col-inds k) j)
		 (return
		   (setf (aref (crs-store A) (aref (crs-offsets pattern) k)) val)))
	  finally (error "i/j not in pattern"))))

(defmethod matrix-ref ((A <crs-matrix>) i &optional j)
  (mat-ref A i j))

(defmethod (setf matrix-ref) (val (A <crs-matrix>) i &optional j)
  (setf (mat-ref A i j) val))

;;; read-only access
(defmethod row-fold ((A <crs-matrix>) (i fixnum) (fn function))
  (let* ((pattern (crs-pattern A))
	 (row-starts (crs-row-starts pattern)))
    (loop for k from (aref row-starts i) below (aref row-starts (1+ i))
	  do (funcall fn (aref (crs-store A) (aref (crs-offsets pattern) k))))))

;;; transformation to matlisp format
(defun crs->matlisp (A)
  (let ((mlmat (make-float-matrix (nrows A) (ncols A))))
    (dotimes (i (nrows A))
      (dotimes (j (ncols A))
	(setf (mat-ref mlmat i j) (mat-ref A i j))))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; <ring> methods (matrix arithmetic)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; <vector> methods
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(definline crs-apply (operator mat1 mat2)
  (if (not (eq (crs-pattern mat1) (crs-pattern mat1)))
      (error "The addition of different crs is not yet implemented."))
  (make-crs-matrix
    (crs-pattern mat1)
    (funcall operator (crs-store mat1) (crs-store mat2))))

(defmethod vec+ ((mat1 <crs-matrix>) (mat2 <crs-matrix>))
  (crs-apply #'vec+ mat1 mat2))
  
(defmethod vec- ((mat1 <crs-matrix>) (mat2 <crs-matrix>))
  (crs-apply #'vec- mat1 mat2))

(defmethod vec-s* ((val double-float) (mat <crs-matrix>))
  (make-instance '<crs-matrix>
    :pattern (crs-pattern mat)
    :store (vec-s* val (crs-store mat))))
(defmethod vec-s* ((mat <crs-matrix>) (val double-float)) (vec-s* val mat))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; vector blas operations for <crs-matrix>
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun crs-blas-caller (proc) (list proc))
(defun crs-vector-transformer (obj) (list 'crs-store obj))
(initialize-vector-blas-methods
 '<crs-matrix> 'array #'crs-blas-caller
 :vector-transformer #'crs-vector-transformer)

;;; we want to have the blas operations also for setting the whole
;;; store to a constant
(initialize-vector-blas-methods
 '<crs-matrix> 'double-float #'crs-blas-caller
 :vector-transformer #'crs-vector-transformer)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; matrix-vector blas operations
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(definline x-op-Ay (x op A y)
  (let* ((store (crs-store A))
	 (pattern (crs-pattern A))
	 (offsets (crs-offsets pattern))
	 (col-inds (crs-col-inds pattern))
	 (row-starts (crs-row-starts pattern)))
    (dotimes (row (crs-nrows pattern) x)
      (setf (aref x row)
	    (let ((acc (if (eq op '=) 0.0d0 (aref x row))))
	      (loop for k from (aref row-starts row) below (aref row-starts (1+ row))
		    for val = (* (aref y (aref col-inds k))
				 (aref store (aref offsets k)))
		    do (if (eq op '-)
			   (decf acc val)
			 (incf acc val))
		    finally (return acc)))))))

(defmethod x<-Ay ((x array) (A <crs-matrix>) (y array))
  (x-op-Ay x '= A y))
(defmethod x+=Ay ((x array) (A <crs-matrix>) (y array))
  (x-op-Ay x '+ A y))
(defmethod x-=Ay ((x array) (A <crs-matrix>) (y array))
  (x-op-Ay x '- A y))

#|
;;; Testing
(let* ((crs-pat (make-crs-pattern 2 5 '( ((* . 0) (a . 4))  ((a . 1)) )))
       (store1 (list->double-vec '(2.0d0 2.0d0)))
       (A (make-crs-matrix crs-pat store1))
       (store2 (list->double-vec '(1.0d0 1.0d0)))
       (B (make-crs-matrix crs-pat store2))
       (x (list->double-vec '(2.0d0 2.0d0)))
       (y (list->double-vec '(1.0d0 1.0d0 1.0d0 1.0d0 1.0d0))))
  ;(describe (vec+ A B))
  ;(describe (vec- A B))
  ;(describe (x+=Ay x A y))
  (describe A)
  (x<-0 A)
  (describe A)
  )
|#
