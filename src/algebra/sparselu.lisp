;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; sparselu.lisp
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

(in-package :fl.algebra)

(defclass <ldu-sparse> ()
  ((lower-left :reader lower-left :initarg :lower-left)
   (diagonal :reader diagonal :initarg :diagonal)
   (upper-right :reader upper-right :initarg :upper-right)
   (ordering :reader ordering :initarg :ordering :type vector)))

(defun modify-diagonal (D R omega)
  "Implements the diagonal modification for ILU_mod."
  (dorows (rk R)
    (let ((sum 0.0))
      (for-each-key-in-row
       #'(lambda (ck)
	   (incf sum (mref R rk ck)))
       R rk)
      (decf (mref D rk rk) (* omega sum))))
  D)

(defun shift-diagonal-inverter (eta)
  "Can be used for obtainint a diagonal modification to get ILU_mod."
  #'(lambda (mat)
      (dorows (i mat)
	(incf (mref mat i i) eta))
      (m/ mat)))

(defmethod sparse-ldu ((A <sparse-matrix>) &key ordering incomplete (omega 0.0) diagonal-inverter)
  (declare (type double-float omega))
  
  ;; default is all keys of A
  (setq ordering (coerce (or ordering (row-keys A)) '(simple-array * (*))))
  
  (let ((L (make-full-block-analog A))
	(D (make-full-block-analog A))
	(U (make-full-block-analog A)))
    (declare (type <sparse-matrix> L D U))
    
    ;; Extract part out of A into U.  Entries are transformed into matlisp
    ;; matrices.  Zero entries are dropped.
    (loop for row-key across ordering
	  for entry = (mref A row-key row-key) do
	  (assert entry)
	  (setf (mref U row-key row-key)
		(typecase entry
		  (crs-matrix (crs->matlisp entry))
		  (t (copy entry)))))
    (loop for row-key across ordering do
	  (for-each-key-and-entry-in-row
	   #'(lambda (col-key entry)
	       (when (and (matrix-row U col-key) (not (mzerop entry)))
		 (setf (mref U row-key col-key)
		       (typecase entry
			 (crs-matrix (crs->matlisp entry))
			 (t (copy entry))))))
	   A row-key))

    ;; decomposition
    (loop
     for k across ordering do
     (let* ((U_kk (mref U k k))
	    (D_kk (if diagonal-inverter
		      (funcall diagonal-inverter U_kk)
		      (m/ U_kk))))
     (setf (mref D k k) D_kk)	; store D
     (remove-entry U k k)		; clean U
     (let ((col-k (matrix-column U k))
	   (row-k (matrix-row U k)))
       (when col-k
	 (loop for i being each hash-key of col-k
	       and Uik being each hash-value of col-k
	       unless (matrix-row D i) do
	       (let ((factor (m* Uik D_kk))
		     (row-i (matrix-row U i)))
		 (setf (mref L i k) factor) ; store L
		 (remove-entry U i k)	; clean U
		 (when row-k
		   (loop for j being each hash-key of row-k
			 and Ukj being each hash-value of row-k do
			 (let ((mblock (gethash j row-i)))
			   (cond ((or mblock (not incomplete))
				  (setq mblock (or mblock (mref U i j)))
				  ;; (decf (vref mblock 0) (* (vref factor 0) (vref Ukj 0)))
				  (gemm! -1.0 factor Ukj 1.0 mblock))
				 (t (unless (zerop omega)
				      (modify-diagonal (mref U i i) (m* factor (mref U i j)) omega)))
				 ))))))))))
    ;; finally return the result
    (make-instance '<ldu-sparse> :lower-left L :diagonal D
		   :upper-right U :ordering ordering)))


;;; Matlisp/LAPACK interface

(defmethod getrf! ((A <sparse-matrix>) &optional ipiv)
  (assert (null ipiv))
  (values (sparse-ldu A :ordering nil) nil t))

(defmethod getrs! ((ldu <ldu-sparse>) (result <sparse-vector>) &optional ipiv)
  (assert (null ipiv))
  ;; solve (I + L) result~ = rhs
  (loop with L of-type <sparse-matrix> = (lower-left ldu)
	for i across (ordering ldu)
	for result-i = (vref result i)
	when (matrix-row L i) do
	(for-each-key-and-entry-in-row
	 #'(lambda (j Lij) (gemm! -1.0 Lij (vref result j) 1.0 result-i))
	 L i))
  ;; solve (D + U) result = result~ (where D^{-1} is stored)
  (loop with U of-type <sparse-matrix> = (upper-right ldu)
	with D of-type <sparse-matrix> = (diagonal ldu)
	for i across (reverse (ordering ldu))
	for result-i = (vref result i) do
	(when (matrix-row U i)
	  (for-each-key-and-entry-in-row
	   #'(lambda (j Uij) (gemm! -1.0 Uij (vref result j) 1.0 result-i))
	   U i))
	(copy! (m* (mref D i i) result-i) result-i))
  
  ;; return the result
  (values result ipiv t))

(defmethod m* ((ldu <ldu-sparse>) (rhs <sparse-vector>))
  "An ldu-decomposition is considered as an inverse."
  (getrs ldu rhs))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Conversion to CCS format
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun numbering (keys &optional (count 0))
  "Returns a hash-table mapping each element of @arg{keys} to an index."
  (let ((numbering (make-hash-table)))
    (map nil #'(lambda (x)
		 (setf (gethash x numbering) count)
		 (incf count))
	 keys)
    numbering))

(defun sparse-matrix->ccs (A &key keys row-keys col-keys ranges row-ranges col-ranges)
  "Converts the sparse matrix @arg{A} to CCS format.  @arg{row-keys} and
@arg{col-keys} may denote a submatrix, @arg{col-ranges} and
@arg{row-ranges} may be used for extracting even subblocks of the entries.
This is a rather difficult routine, which might suggest switching to CCS
completely."
  (ensure row-keys keys) (ensure col-keys keys)
  (unless row-keys
    (setq row-keys (coerce (row-keys A) 'vector))
    (ensure col-keys (if (automorphism? A)
			 row-keys
			 (coerce (col-keys A) 'vector))))
  (ensure row-ranges ranges) (ensure col-ranges ranges)
  (let ((row-numbering (numbering row-keys))
	row-sizes col-sizes row-offsets col-offsets)
    (flet ((setup (keys ranges key->size)
	     (let ((sizes (if ranges
			      (vector-map #'(lambda (x) (- (cdr x) (car x))) ranges)
			      (vector-map key->size keys))))
	       (values sizes (concatenate 'vector #(0) (partial-sums sizes))))))
      (multiple-value-setq (row-sizes row-offsets)
	(setup row-keys ranges (row-key->size A)))
      (multiple-value-setq (col-sizes col-offsets)
	(if (and (automorphism? A) (equalp row-keys col-keys) (equalp row-ranges col-ranges))
	    (values row-sizes row-offsets)
	    (setup col-keys col-ranges (col-key->size A)))))
    (let* ((nrows (vector-last row-offsets))
	   (ncols (vector-last col-offsets))
	   (n (loop for ck across col-keys
		 and cs across col-sizes summing
		 (* cs (let ((count 0))
			 (for-each-key-in-col
			  #'(lambda (rk)
			      (whereas ((row-index (gethash rk row-numbering)))
				(incf count (aref row-sizes row-index))))
			  A ck)
			 count))))
	   (column-starts (make-int-vec (1+ ncols)))
	   (row-indices (make-int-vec n))
	   (store (make-double-vec n))
	   (column-index 0) (pos 0))
      (setf (aref column-starts ncols) n)
      ;; loop through columns
      (loop
	 for ci from 0 and ck across col-keys do
	 ;; sort the block columns ...
	 (let* ((sorted-block-column
		 (sort (remove nil (mapcar #'(lambda (x) (gethash x row-numbering))
					   (hash-table-keys (matrix-column A ck))))
		       #'<))
		(column-width (aref col-sizes ci))
		(col-range-start (if col-ranges
				     (car (aref col-ranges ci))
				     0))
		(column-length
		 (loop for ri in sorted-block-column
		    summing (aref row-sizes ri))))
	   (declare (type fixnum column-width col-range-start column-length))
	   ;; and put in the ccs store
	   (loop with pos1 of-type fixnum = pos
		 for ri in sorted-block-column do
		 (let* ((row-range-start
			 (if row-ranges (car (aref row-ranges ri)) 0))
			(row-offset (aref row-offsets ri))
			(entry (mref A (aref row-keys ri) ck))
			(entry-store (store entry)))
		   (declare (type fixnum row-range-start row-offset)
			    (type double-vec entry-store))
		   (dotimes (i (the fixnum (aref row-sizes ri)))
		     (declare (type fixnum i))
		     (declare (optimize speed (safety 0)))
		     (dotimes (j column-width)
		       (declare (type fixnum j))
		       (let ((k (the fixnum (+ pos1 (the fixnum (* j column-length))))))
			 (declare (type fixnum k))
			 (setf (aref row-indices k) (+ i row-offset))
			 (setf (aref store k)
			       (aref entry-store (fl.matlisp::standard-matrix-indexing
						  (the fixnum (+ i row-range-start))
						  (the fixnum (+ j col-range-start))
						  (the fixnum (nrows entry)))))))
		     (setf pos1 (the fixnum (+ pos1 1))))))
	   ;; set column-starts
	   (loop for i below column-width do
		(setf (aref column-starts column-index)
		      (+ pos (* i column-length)))
		(incf column-index))
	   ;; next block column
	   (incf pos (* column-length column-width))))
      ;; return the result
      (values
       (let ((pattern (make-instance
		       'ccs-pattern :nrows nrows :ncols ncols
		       :column-starts column-starts :row-indices row-indices)))
	 (make-instance (ccs-matrix 'double-float) :pattern pattern :store store))
       row-keys col-keys))))

(defun key->index (mat keys type)
  "Returns a hash-table mapping keys to CCS offsets."
  (let ((numbering (make-hash-table))
	(k 0))
    (map nil #'(lambda (ck)
		 (setf (gethash ck numbering) k)
		 (incf k (funcall (ecase type
				    (:row (row-key->size mat))
				    (:column (col-key->size mat)))
				  ck)))
	 keys)
    (values numbering k)))

#+(or)  ; old: only for testing purposes
(defun sparse-matrix-to-ccs (smat keys)
  "Converts the sparse-matrix @arg{smat} to ccs format.  @arg{keys} is a
list of index keys, @arg{ranges} can be nil or a list of ranges of the same
length as @arg{keys} which specify subblocks of the block unknowns
associated with each key."
  (let* ((n (total-entries smat))
	 (nrows (total-nrows smat))
	 (column-starts (make-uint-vec (1+ nrows)))
	 (row-indices (make-uint-vec n))
	 (store (make-double-vec n))
	 (block-numbering (key->index smat keys :row))
	 (column-index 0) (pos 0))
    (setf (aref column-starts nrows) n)
    ;; loop through columns (hopefully in increasing order)
    (dolist (ck keys)
      ;; convert the column in a sorted list
      (let* ((sorted-block-column
	      (sort (hash-table-keys (matrix-column smat ck))
		    #'< :key (lambda (x) (gethash x block-numbering))))
	     column-width
	     (column-length
	      (loop for rk in sorted-block-column
		    for entry = (mref smat rk ck)
		    do (ensure column-width (ncols entry))
		    summing (nrows entry))))
	(loop with pos1 = pos
	      for rk in sorted-block-column
	      for entry = (mref smat rk ck) do
	      (dotimes (i (nrows entry))
		(dotimes (j (ncols entry))
		  (let ((k (+ pos1 (* j column-length))))
		    (setf (aref row-indices k) (+ i (gethash rk block-numbering)))
		    (setf (aref store k) (mref entry i j))))
		(incf pos1)))
	;; set column-starts
	(loop for i below column-width do
	      (setf (aref column-starts column-index)
		    (+ pos (* i column-length)))
	      (incf column-index))
	;; next block column
	(incf pos (* column-length column-width))))
    (let ((pattern (make-instance
		    'ccs-pattern :nrows nrows :ncols nrows 
		    :column-starts column-starts :row-indices row-indices)))
      (make-instance (ccs-matrix 'double-float) :pattern pattern :store store))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; GESV!
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

#+(or umfpack superlu)  ; otherwise there is no GESV! for ccs-matrices
(defmethod gesv! ((smat <sparse-matrix>) (svec <sparse-vector>))
  "Solve the system by calling an external sparse solver."
  (let* ((keys (row-keys smat))
	 (ccs (sparse-matrix->ccs smat :keys (coerce keys 'vector)))
	 (cvec (sparse-vector->matlisp svec keys)))
    (gesv! ccs cvec)
    (set-svec-to-local-block svec cvec keys)
    svec))

;;; Testing:

(defun test-sparselu ()
  
  (let ((A (make-sparse-matrix
	    :row-key->size (constantly 1) :col-key->size (constantly 1)
	    :keys->pattern (constantly (full-crs-pattern 1 1)))))
    (describe A)
    (setf (mref A 0 0) #m((2.0)))
    (setf (mref A 1 0) #m((1.0)))
    (setf (mref A 1 1) #m((1.0)))
    (setf (mref A 0 1) #m((1.0)))
    (describe (sparse-matrix->ccs A :row-keys #(0) :col-keys #(0 1)))
    (describe (sparse-matrix->ccs A :keys #(1)))
    (display A)
    (display (lower-left (sparse-ldu A)))
    (display (diagonal (sparse-ldu A)))
    (display (upper-right (sparse-ldu A)))
    (ordering (sparse-ldu A))

    (let ((rhs (make-instance '<sparse-vector> :key->size (constantly 1))))
      (setf (vref rhs 0) #m((3.0)))
      (setf (vref rhs 1) #m((1.0)))
      (show rhs)
      (show (getrs (sparse-ldu A) rhs))))

  (let ((A (make-sparse-matrix
	    :row-key->size (constantly 2) :col-key->size (constantly 2)
	    :keys->pattern (constantly (full-crs-pattern 2 2)))))
    (describe A)
    (setf (mref A 0 0) #m((4.0 -1.0) (-1.0 4.0)))
    (setf (mref A 1 0) #m((-1.0 0.0) (0.0 -1.0)))
    (setf (mref A 1 1) #m((4.0 -1.0) (-1.0 4.0)))
    (setf (mref A 0 1) #m((-1.0 0.0) (0.0 -1.0)))
    
    (display A)
    (display (lower-left (sparse-ldu A)))
    (display (diagonal (sparse-ldu A)))
    (display (upper-right (sparse-ldu A)))
    (ordering (sparse-ldu A))
    (let ((rhs (make-instance '<sparse-vector> :key->size (constantly 2))))
      (setf (vref rhs 0) #m((1.0) (2.0)))
      (setf (vref rhs 1) #m((3.0) (4.0)))
      (show rhs)
      (mzerop (m- (m* A (getrs (sparse-ldu A) rhs)) rhs) 1.0e-15)
      (gesv! A rhs)
      (show rhs)
      (show (m* A rhs))))
  )