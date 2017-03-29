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

(in-package :fl.matlisp)

(defclass <ldu-sparse> ()
  ((lower-left :reader lower-left :initarg :lower-left)
   (diagonal :reader ldu-diagonal :initarg :diagonal)
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
                 (compressed-matrix (compressed->matlisp entry))
                 (t (copy entry)))))
    (loop for row-key across ordering do
         (for-each-key-and-entry-in-row
          #'(lambda (col-key entry)
              (when (and (matrix-row U col-key) (not (mzerop entry)))
                (setf (mref U row-key col-key)
                      (typecase entry
                        (compressed-matrix (compressed->matlisp entry))
                        (t (copy entry))))))
          A row-key))

    ;; decomposition
    (dovec (k ordering)
      (let* ((U_kk (mref U k k))
             (D_kk (if diagonal-inverter
                       (funcall diagonal-inverter U_kk)
                       (m/ U_kk))))
        (setf (mref D k k) D_kk)	; store D
        (remove-key U k k)		; clean U
        (let ((col-k (matrix-column U k))
              (row-k (matrix-row U k)))
          (flet
              ((work-on (i)
                 (let ((factor (mref L i k))
                       (row-i (matrix-row U i)))
                   (when row-k
                     (dodic ((j Ukj) row-k)
                       (let ((mblock (gethash j row-i)))
                         (cond ((or mblock (not incomplete))
                                (setq mblock (or mblock (mref U i j)))
                                ;; (decf (vref mblock 0) (* (vref factor 0) (vref Ukj 0)))
                                (gemm! -1.0 factor Ukj 1.0 mblock))
                               (t (unless (zerop omega)
                                    (modify-diagonal (mref U i i) (m* factor (mref U i j)) omega)))
                               )))))))
            (when col-k
              (dodic (i col-k)
                (unless (matrix-row D i)
                  (let* ((Uik (dic-ref col-k i))
                         (factor (m* Uik D_kk)))
                    (setf (mref L i k) factor) ; store L
                    (remove-key U i k)	; clean U
                    (work-on i)))))))))
    
    ;; finally return the result
    (make-instance '<ldu-sparse> :lower-left L :diagonal D
		   :upper-right U :ordering ordering)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; nested disection numbering
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun distance-numbering (mat key-table key)
  "Result is a table mapping the keys from @arg{key-table}
to their distance from @arg{key}."
  (lret ((result (make-instance 'sorted-hash-table)))
    (when (dic-ref key-table key)
      (setf (dic-ref result key) 0)
      (dic-for-each
       (lambda (key index)
         (for-each-key-in-row
          (lambda (rk)
            (when (dic-ref key-table rk)
              (unless (dic-ref result rk)
                (setf (dic-ref result rk) (1+ index)))))
          mat key))
       result))))

(defun last-key (sorted-ht)
  "Returns the last key in a sorted hash table.
This function should be in the interface of the sorted hash table."
  (car (dll-peek-last (slot-value sorted-ht 'fl.dictionary::store))))

(defun find-ordering (mat keys)
  "Finds a nested disection ordering of @arg{keys}
with respect to the graph given by @arg{mat}."
  (if (<= (length keys) 2)
      keys
      (let* ((key-table (map-list-in-hash-table (_ (values _ t)) keys))
             (key (first keys))
             (middle-distance (distance-numbering mat key-table key))
             (left (last-key middle-distance))
             (left-distance (distance-numbering mat key-table left))
             (right (last-key left-distance))
             (right-distance (distance-numbering mat key-table right)))
        (unless (= (length keys)
                   (length (keys left-distance))
                   (length (keys middle-distance))
                   (length (keys right-distance)))
          ;; algorithm does not work
          (dbg :sparselu "Giving up finding a good ordering - need better algorithm for unsymmetric matrices")
          (return-from find-ordering keys))
        (let (interface left right)
          (dic-for-each-key
           (lambda (key)
             (let ((k (dic-ref left-distance key))
                   (l (dic-ref right-distance key)))
               (cond 
                 ((or (= k l) (= (1+ k) l)) (push key interface))
                 ((< k l) (push key left))
                 (t (push key right)))))
           key-table)
          ;;
          (append (find-ordering mat left)
                  (find-ordering mat right)
                  interface)))))

;;; Matlisp/LAPACK interface

(defmethod getrf! ((A <sparse-matrix>) &optional ipiv)
  (ensure ipiv (find-ordering A (row-keys A)))
  (dbg-show :sparselu A)
  (values (sparse-ldu A :ordering ipiv) ipiv t))

(defmethod getrs! ((ldu <ldu-sparse>) (result <ht-sparse-vector>) &optional ipiv)
  (with-slots (lower-left upper-right diagonal ordering) ldu
    ;; solve (I + L) result~ = rhs
    (loop with L of-type <sparse-matrix> = lower-left
          for i across ordering
          for result-i = (vref result i)
          when (matrix-row L i) do
            (for-each-key-and-entry-in-row
             #'(lambda (j Lij) (gemm! -1.0 Lij (vref result j) 1.0 result-i))
             L i))
    ;; solve (D + U) result = result~ (where D^{-1} is stored)
    (loop with U of-type <sparse-matrix> = upper-right
          with D of-type <sparse-matrix> = diagonal
          for i across (reverse ordering)
          for result-i = (vref result i) do
            (when (matrix-row U i)
              (for-each-key-and-entry-in-row
               #'(lambda (j Uij) (gemm! -1.0 Uij (vref result j) 1.0 result-i))
               U i))
            (copy! (m* (mref D i i) result-i) result-i))
  
    ;; return the result
    (values result ipiv t)))

(defmethod m* ((ldu <ldu-sparse>) (rhs <ht-sparse-vector>))
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

This is a rather complicated routine which has not yet been parallelized.
Theoretically, this might be a bottleneck for some applications, but,
practically, an appropriate case has not yet appeared."
  (ensure row-keys keys) (ensure col-keys keys)
  (unless row-keys
    (setq row-keys (coerce (row-keys A) 'vector))
    (ensure col-keys (if (automorphism-p A)
			 row-keys
			 (coerce (col-keys A) 'vector))))
  (ensure row-ranges ranges) (ensure col-ranges ranges)
  (let ((row-numbering (numbering row-keys))
	row-sizes col-sizes row-offsets col-offsets)
    ;; we ensure that row-sizes and col-sizes are set to their correct values
    (flet ((setup (keys ranges key->size)
	     (let ((sizes (if ranges
			      (vector-map #'(lambda (x) (- (cdr x) (car x))) ranges)
			      (vector-map key->size keys))))
	       (values sizes (concatenate 'vector #(0) (partial-sums sizes))))))
      (multiple-value-setq (row-sizes row-offsets)
	(setup row-keys ranges (row-key->size A)))
      (multiple-value-setq (col-sizes col-offsets)
	(if (and (automorphism-p A) (equalp row-keys col-keys) (equalp row-ranges col-ranges))
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
		       'compressed-pattern :sizes (vector ncols nrows)
		       :starts column-starts :indices row-indices
		       :orientation :column)))
	 (make-instance (compressed-matrix 'double-float) :pattern pattern :store store))
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

#+(or) ; old: only for testing purposes
(defun sparse-matrix-to-ccs-old (smat keys)
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
		    'compressed-pattern :sizes (vector ncols nrows)
		    :column-starts column-starts :row-indices row-indices)))
      (make-instance (compressed-matrix 'double-float) :pattern pattern :store store))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; GESV!
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

#+(or umfpack superlu)  ; otherwise there is no GESV! for ccs-matrices
(defmethod gesv! ((smat <sparse-matrix>) (svec <ht-sparse-vector>))
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
    (display (ldu-diagonal (sparse-ldu A)))
    (display (upper-right (sparse-ldu A)))
    (ordering (sparse-ldu A))

    (let ((rhs (make-instance '<ht-sparse-vector> :key->size (constantly 1))))
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
    (display (ldu-diagonal (sparse-ldu A)))
    (display (upper-right (sparse-ldu A)))
    (ordering (sparse-ldu A))
    (let ((rhs (make-instance '<ht-sparse-vector>)))
      (setf (vref rhs 0) #m((1.0) (2.0)))
      (setf (vref rhs 1) #m((3.0) (4.0)))
      (show rhs)
      (show (copy rhs))
      (mzerop (m- (m* A (getrs (sparse-ldu A) rhs)) rhs) 1.0e-15)
      (gesv! A rhs)
      (show rhs)
      (show (m* A rhs))))

  (let ((A (make-sparse-matrix
	    :row-key->size (constantly 2) :col-key->size (constantly 2)
	    :keys->pattern (constantly (full-crs-pattern 2 2)))))
    (setf (mref A 0 0) #m((1.0 2.0) (3.0 4.0)))
    (let* ((block #(0))
	   (ranges nil)
	   (ccs (sparse-matrix->ccs A :keys block :ranges ranges))
	   (m1 (compressed->matlisp ccs))
	   (m2 (sparse-matrix->matlisp A :keys block :ranges ranges)))
      (unless (mequalp m1 m2)
	(mat-diff m1 m2)
	(error "Not equal"))))
  )
