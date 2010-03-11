;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; sparse-matrix.lisp
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;
;;; Copyright (C) 2003-2006 Nicolas Neuss, University of Heidelberg.
;;; Copyright (C) 2007- Nicolas Neuss, University of Karlsruhe.
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

;;; (declaim (optimize (safety 3) (debug 3)))

(in-package :fl.matlisp)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; <sparse-matrix>
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass <sparse-matrix> (<matrix>) ()
  (:documentation "Abstract class for sparse matrices."))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; sparse matrix format definition
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmethod element-type ((smat <sparse-matrix>))
  (standard-matrix 'double-float))

(defmethod scalar-type ((smat <sparse-matrix>))
  'double-float)

(defgeneric row-key->size (smat)
  (:method ((smat <sparse-matrix>)) nil))
(defgeneric col-key->size (smat)
  (:method ((smat <sparse-matrix>)) nil))
(defgeneric keys->pattern (smat)
  (:method ((smat <sparse-matrix>)) nil))

(defgeneric access-type (mat &key suggest require)
  (:documentation "If @arg{suggest} and @arg{require} are NIL, returns which of
  @symbol{:row} or @symbol{:column} is prefered for @arg{mat}.  Otherwise
  determine if the order in @arg{suggest} or @arg{require} is acceptable without
  serious performance hit.")
  (:method (mat &key suggest require)
    "The default allows the suggested access or selects row-wise access."
    (declare (ignore mat))
    (or require suggest :row)))

(defclass <block-definition-mixin> ()
  ((row-key->size :reader row-key->size :initarg :row-key->size :type function)
   (col-key->size :reader col-key->size :initarg :col-key->size :type function)
   (keys->pattern :reader keys->pattern :initarg :keys->pattern
                  :type function :documentation
		  "Function determining the pattern of a block.")))

;;; also fl.multiprocessing::locked-region-mixin can be included here

(defmethod make-analog ((smat <block-definition-mixin>))
  (copy-slots (call-next-method) smat
              '(row-key->size col-key->size keys->pattern)))

;;; make-sparse-matrix is defined later

(defun make-sparse-automorphism (&key key->size keys->pattern)
  (make-sparse-matrix
   :keys->pattern keys->pattern :row-key->size key->size :col-key->size key->size))

(defun make-full-block-analog (sm)
  (make-sparse-matrix
   :row-key->size (row-key->size sm)
   :col-key->size (col-key->size sm)
   :keys->pattern
   #'(lambda (row-key col-key)
       (full-crs-pattern (funcall (row-key->size sm) row-key)
			 (funcall (row-key->size sm) col-key)))))

(defgeneric automorphism-p (smat)
  (:documentation "Returns T, if domain and range are equal.")
  (:method ((mat <matrix>))
    "The default method checks only if the matrix is quadratic."
    (= (nrows mat) (ncols mat)))
  (:method ((A <block-definition-mixin>))
    (and (call-next-method)
         (eql (row-key->size A) (col-key->size A)))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; sparse matrix entry access
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defgeneric matrix-block (smat row-key col-key)
  (:documentation "Low-level block lookup for a sparse block matrix."))

(defgeneric (setf matrix-block) (value smat row-key col-key)
  (:documentation "Low-level block insertion into a sparse block matrix."))

(defmethod in-pattern-p ((smat <sparse-matrix>) &rest indices)
  (apply #'matrix-block smat indices))

(defmethod entry-allowed-p ((smat <sparse-matrix>) &rest indices)
  (destructuring-bind (row-key col-key) indices
    (plusp (number-nonzero-entries (funcall (keys->pattern smat)
					    row-key col-key)))))

(defun make-matrix-block (smat row-key col-key)
  (declare (type <sparse-matrix> smat))
  (let* ((pattern (funcall (keys->pattern smat) row-key col-key))
         (sizes (slot-value pattern 'sizes))
         (store-size (number-nonzero-entries pattern)))
    (cond
      ((some #'zerop sizes)
       (error "Request for an empty matrix."))
      ((= store-size (reduce #'* sizes))
       (make-real-matrix (aref sizes 0) (aref sizes 1)))
      (t
       ;; warning: the following use of make-instance with
       ;; keyword parameters on a runtime-computed class might be
       ;; a performance problem
       (make-instance (compressed-matrix 'double-float)
                      :pattern pattern)))))

(defmethod mref ((smat <sparse-matrix>) row-key col-key)
  "If the matrix-block indexed by row-key/col-key exists, it is returned.
Otherwise, a new matrix block is created according to the pattern specified
for this index."
  (or (matrix-block smat row-key col-key)
      (setf (matrix-block smat row-key col-key)
            (make-matrix-block smat row-key col-key))))

(defmethod (setf mref) (entry (smat <sparse-matrix>) row-key col-key)
  (let ((rks (row-key->size smat))
        (cks (col-key->size smat)))
    (assert (and (or (null rks) (= (nrows entry) (funcall rks row-key)))
                 (or (null cks) (= (ncols entry) (funcall cks col-key))))))
  (setf (matrix-block smat row-key col-key) entry))

(defgeneric remove-row (smat row-key)
  (:documentation "Removes a row of @arg{smat}.")
  (:method ((smat <sparse-matrix>) row-key)
    (loop for ck in (keys-of-row smat row-key)
       do (remove-key smat row-key ck))))

(defgeneric remove-column (smat col-key)
  (:documentation "Removes a column of @arg{smat}.")
  (:method ((smat <sparse-matrix>) col-key)
    (loop for rk in (keys-of-column smat col-key)
       do (remove-key smat rk col-key))))

(defmethod row<-id ((A <sparse-matrix>) key)
  (remove-row A key)
  (unless (and (row-key->size A) (col-key->size A))
    "Cannot determine size of diagonal block.")
  (setf (mref A key key)
	(eye (funcall (row-key->size A) key)
	     (funcall (col-key->size A) key))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; derived methods
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmethod show ((smat <sparse-matrix>) &key keys (zeros t) &allow-other-keys)
  (format t "~&Sparse matrix:~%")
  (dolist (row-key (or keys (row-keys smat)))
    (format t "~&*** ~A ***~%" row-key)
    (for-each-key-and-entry-in-row 
     (lambda (col-key vblock)
       (when (or zeros (not (mzerop vblock)))
         (format t "~A -> ~A~%" col-key vblock)))
     smat row-key)))

(defgeneric total-nrows (mat)
  (:documentation "Total number of rows for a matrix (works also for block
  matrices).")
  (:method (mat) (nrows mat))
  (:method ((mat <sparse-matrix>))
    (lret ((sum 0))
      (for-each-row-key
       #'(lambda (row-key)
           (incf sum
                 (aif (row-key->size mat) 
                      (funcall (row-key->size mat) row-key)
                      (total-nrows (mapper-select-first
                                    #'for-each-entry-in-row
                                    mat row-key)))))
       mat))))

(defmethod symmetric-p ((smat <sparse-matrix>) &key (threshold 0.0) output)
  (and (automorphism-p smat)
       (let ((flag t))
         (dovec ((entry rk ck) smat)
           (let ((entry2 (mref smat ck rk)))
             (when (> (norm (m- entry2 (transpose entry))) threshold)
               (when output
                 (format t "~&Mismatch~%(i,j)=(~A,~A):~%Aij=~A~%Aji=~A~%~%"
                         rk ck entry entry2))
               (setq flag nil))))
         flag)))

(defmethod transpose! ((x <sparse-matrix>) (y <sparse-matrix>))
  (dovec ((entry rk ck) x y)
    (setf (mref y ck rk) (transpose entry))))

(defmethod extract-value-blocks ((smat <sparse-matrix>) row-keys &optional col-keys)
  (lret ((result (make-array (list (length row-keys)
                                   (length col-keys)))))
    (loop+ (i (rk row-keys)) do
       (loop+ (j (ck col-keys))
          do (setf (aref result i j)
                   (and (entry-allowed-p smat rk ck)
                        (mref smat rk ck)))))))

(defmethod mat-diff ((smat1 <sparse-matrix>) (smat2 <sparse-matrix>))
  (format t "Missing in [2]~%")
  (dovec ((entry i j) smat1)
    (unless (matrix-block smat2 i j)
      (format t "(~A,~A):~%~A~%" i j entry)))
  (format t "Missing in [1]~%")
  (dovec ((entry i j) smat2)
    (unless (matrix-block smat1 i j)
      (format t "(~A,~A):~%~A~%" i j entry)))
  (format t "Differences~%")
  (dovec ((entry i j) smat1)
    (whereas ((entry2 (matrix-block smat2 i j)))
      (unless (mzerop (m- entry entry2))
	(format t "(~A,~A) :~%" i j)
	(format t "~A~%~A~%" entry entry2)))))

(defmethod display ((smat <sparse-matrix>) &key row-order col-order order)
  (ensure row-order (or order (row-keys smat)))
  (let ((row-indices (let ((result (make-hash-table)))
		       (loop for i from 0 and key in row-order do
			     (setf (gethash key result) i)
			     finally (return result)))))
    (unless col-order
      (setq col-order (or order (col-keys smat)))
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

;;; Matrix-vector stuff

(defmethod make-image-vector-for ((A <sparse-matrix>) &optional (multiplicity 1))
  (lret ((result (make-sparse-vector :key->size (row-key->size A)
                                     :multiplicity multiplicity)))
    (for-each-row-key (curry #'vref result) A)))

(defmethod make-domain-vector-for ((A <sparse-matrix>) &optional (multiplicity 1))
  (lret ((result (make-sparse-vector :key->size (col-key->size A)
                                     :multiplicity multiplicity)))
    (for-each-col-key (curry #'vref result) A)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Transformation of unknowns
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; This is needed by hanging nodes, identified boundaries, ...

(defun matrix-row-p (mat row-key)
   (mapper-select-first #'for-each-key-in-row mat row-key))

(defun matrix-col-p (mat col-key)
   (mapper-select-first #'for-each-key-in-col mat col-key))

(defun index-range-disjoint-p (mat1 mat2)
  "Checks if the range of indices of two sparse matrices is disjoint."
  (not (mapper-some (curry #'matrix-row-p mat2)
                    #'for-each-row-key mat1)))

(defun range-and-domain-disjoint-p (mat)
  "Checks if index range and index domain of some matrix are disjoint."
  (not (mapper-some (curry #'matrix-col-p mat)
                    #'for-each-row-key mat)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; vector blas operations for the <smat> class
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; those are defined via the general methods in vector.lisp

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; matlisp operations for the <smat> class
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmethod m*-product-instance ((A <sparse-matrix>) (y <sparse-vector>))
  (make-image-vector-for A (multiplicity y)))

#+(or)
(defun generate-sparse-matrix-vector-gemm!-template (job)
  "Generates the GEMM-XX! routine defined by JOB."
  (assert (member job '(:nn :nt :tn :tt)))
  (let ((gemm-job (symconc "GEMM-" (symbol-name job) "!")))
    (eval
     `(defmethod ,gemm-job
       (alpha (A <sparse-matrix>) (y <sparse-vector>) beta (x <sparse-vector>))
       (declare (optimize (speed 3) (space 2) (safety 1)))
       (with-workers ((lambda (x-key)
			(let ((x-values (vref x x-key)))
			  (scal! beta x-values)
			  (,(if (member job '(:nn :nt))
                                'for-each-key-and-entry-in-row
                                'for-each-key-and-entry-in-col)
                            #'(lambda (y-key mblock)
                                (let ((y-values (vref y y-key)))
                                  (,gemm-job alpha mblock y-values 1.0 x-values)))
                            A x-key))))
	 (,(if (member job '(:nn :nt))
               'for-each-row-key
               'for-each-col-key)
           #'work-on A))
       x))))

(defun generate-sparse-matrix-vector-gemm!-template (job)
  "Generates the GEMM-XX! routine defined by JOB."
  (assert (member job '(:nn :nt :tn :tt)))
  (let ((gemm-job (symconc "GEMM-" (symbol-name job) "!")))
    (eval
     `(defmethod ,gemm-job
       (alpha (A <sparse-matrix>) (y <sparse-vector>) beta (x <sparse-vector>))
       (declare (optimize (speed 3) (space 2) (safety 1)))
       ,(ecase
         job
         ((:nn :nt)
          `(ecase (access-type A :suggest :row)
             (:row (dorows (rk A)
                     (let ((x-values (vref x rk)))
                       (scal! beta x-values)
                       (dorow ((ck mblock) A rk)
                         (let ((y-values (vref y ck)))
                           (,gemm-job alpha mblock y-values 1.0 x-values))))))
             (:column
              (scal! beta x)
              (docols (ck A)
                (let ((y-values (vref y ck)))
                  (dorow ((rk mblock) A ck)
                    (,gemm-job alpha mblock y-values 1.0 (vref x rk))))))))
         ((:tn :tt)
          `(ecase (access-type A :suggest :column)
             (:column
              (docols (ck A)
                (let ((x-values (vref x ck)))
                  (scal! beta x-values)
                  (docol ((rk mblock) A ck)
                    (let ((y-values (vref y rk)))
                      (,gemm-job alpha mblock y-values 1.0 x-values))))))
             (:row
                 (scal! beta x)
               (dorows (rk A)
                 (let ((y-values (vref y rk)))
                   (dorow ((ck mblock) A rk)
                     (,gemm-job alpha mblock y-values 1.0 (vref x ck))))))
             )))
       ;; the routine returns x!
       x
       ))))

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
  (lret ((result (ecase job
                   ((:nn :tn) (make-image-vector-for A (multiplicity B)))
                   ((:nt :tt) (make-domain-vector-for A (multiplicity B))))))
    (ecase sparsity
      (:A (funcall
           (ecase job
	     ((:nn :nt) #'for-each-row-key)
	     ((:tn :tt) #'for-each-col-key))
           (lambda (i)
             (funcall (ecase job
                        ((:nn :nt) #'for-each-key-and-entry-in-row)
                        ((:tn :tt) #'for-each-key-and-entry-in-col))
                      (lambda (j Aij)
                        (gemm! 1.0 Aij (vref B j) 1.0 (vref result i) job))
                      A i))
           A))
      (:B (for-each-entry-and-key
           (lambda (bi i)
             (funcall (ecase job
                        ((:nn :nt) #'for-each-key-and-entry-in-col)
                        ((:tn :tt) #'for-each-key-and-entry-in-row))
                      (lambda (j Aij)
                        (gemm! 1.0 Aij bi 1.0 (vref result j) job))
                      A i))
           B)))))

(defmethod matrix-transpose-instance ((smat <sparse-matrix>))
  (make-sparse-matrix
   :keys->pattern #'(lambda (rk ck)
                      (transposed-pattern (funcall (keys->pattern smat) ck rk)))
   :row-key->size (col-key->size smat)
   :col-key->size (row-key->size smat)))

(defmethod sparse-m* ((A <sparse-matrix>) (B <sparse-matrix>)
		      &key (job :nn) (sparsity :A))
  (lret* ((transpose-A? (member job '(:tn :tt)))
          (transpose-B? (member job '(:nt :tt)))
          (row-key->size (if transpose-A? (col-key->size A) (row-key->size A)))
          (col-key->size (if transpose-B? (row-key->size B) (col-key->size B)))
          (C (make-sparse-matrix
              :row-key->size row-key->size :col-key->size col-key->size
              :keys->pattern
              #'(lambda (row-key col-key) ; at the moment this works only for full blocks
                  (full-crs-pattern (funcall row-key->size row-key)
                                    (funcall col-key->size col-key))))))
    (ecase sparsity
      ((:A)
       (funcall
        (ecase job
          ((:nn :nt) #'for-each-row-key)
          ((:tn :tt) #'for-each-col-key))
        #'(lambda (i)
            (funcall (ecase job
                       ((:nn :nt) #'for-each-key-and-entry-in-row)
                       ((:tn :tt) #'for-each-key-and-entry-in-col))
                     #'(lambda (j Aij)
                         (funcall (ecase job
                                    ((:nn :tn) #'for-each-key-and-entry-in-row)
                                    ((:nt :tt) #'for-each-key-and-entry-in-col))
                                  #'(lambda (k Bjk)
                                      (gemm! 1.0 Aij Bjk 1.0 (mref C i k) job))
                                  B j))
                     A i))
        A))
      ((:B)
       (funcall
        (ecase job
          ((:nn :tn) #'for-each-col-key)
          ((:nt :tt) #'for-each-row-key))
        #'(lambda (k)
            (funcall (ecase job
                       ((:nn :tn) #'for-each-key-and-entry-in-col)
                       ((:nt :tt) #'for-each-key-and-entry-in-row))
                     #'(lambda (j Bjk)
                         (funcall (ecase job
                                    ((:nn :nt) #'for-each-key-and-entry-in-col)
                                    ((:tn :tt) #'for-each-key-and-entry-in-row))
                                  #'(lambda (i Aij)
                                      (gemm! 1.0 Aij Bjk 1.0 (mref C i k) job))
                                  A j))
                     B k))
        B)))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; submatrices / extraction
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmethod extract-if ((test function) (smat <sparse-matrix>) &key &allow-other-keys)
  "Extracts a sub-matrix from a sparse matrix."
  (lret ((sub-mat (make-analog smat)))
    (dovec ((entry i j) smat)
      (when (funcall test entry i j)
	(setf (mref sub-mat i j) entry)))))

#+(or) ;deprecated
(defun extract-matrix-block (smat row-keys col-keys)
  "Extracts a sub-matrix from @arg{smat}.  @arg{row-keys} and @arg{col-keys} are
hash-tables of keys (or NIL, which means to allow every key)."
  (extract-if #'(lambda (entry row-key col-key)
		  (declare (ignore entry))
		  (or (and row-keys (not (gethash row-key row-keys)))
		      (and col-keys (not (gethash col-key col-keys)))))
	      smat))

#+(or) ;deprecated
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


(defmethod remove-projection-range ((vec <vector>) projection &key &allow-other-keys)
  (for-each-row-key
   #'(lambda (i)
       (let* ((P (mref projection i i))
	      (S (m- (eye (nrows P)) P)))
	 (setf (vref vec i)
	       (m* S (vref vec i)))))
   projection)
  vec)

(defmethod remove-projection-range ((mat <matrix>) projection &key row-p column-p)
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
;;;; GPS choice of solver
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defvar *maxrows-for-direct-solving* 5000
  "Maximum number of rows for which direct solving is applied.")

(defmethod select-linear-solver :around ((asa <matrix>) blackboard)
  "Select a suitable solver depending on size of the matrix and the pde
problem."
  (declare (ignore blackboard))
  (if (<= (nrows asa) *maxrows-for-direct-solving*)
      (make-instance '<linear-solver> :success-if `(>= :step 1)
		     :iteration (make-instance '<lu> :store-p nil))
      (call-next-method)))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; <sparse-dictionary-matrix>
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass <sparse-dictionary-matrix> (<sparse-matrix>)
  ((row-table :accessor row-table :initarg :row-table
	      :documentation "Table of row dictionaries.")
   (column-table :accessor column-table :initarg :column-table
		 :documentation "Table of column dictionaries."))
  (:documentation "Sparse matrices defined via dictionaries."))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Accessing specific rows and columns as dictionaries
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmethod access-type ((csm <sparse-dictionary-matrix>) &key suggest require)
  (cond ((and (slot-boundp csm 'row-table) (slot-boundp csm 'column-table))
         (or require suggest :row))
        ((slot-boundp csm 'row-table)
         (unless (eql require :column)
           :row))
        ((slot-boundp csm 'column-table)
         (unless (eql require :row)
           :column))))

(defgeneric matrix-row (mat row-key)
  (:documentation "Returns a dictionary mapping col-key to entry.")
  (:method ((mat <sparse-dictionary-matrix>) key)
    (dic-ref (row-table mat) key)))

(defgeneric matrix-column (mat col-key)
  (:documentation "Returns a dictionary mapping row-key to entry.")
  (:method  ((mat <sparse-dictionary-matrix>) key)
    (dic-ref (column-table mat) key)))

(defmethod (setf matrix-row) (row-table (mat <sparse-dictionary-matrix>) key)
  (setf (dic-ref (row-table mat) key)
        row-table))

(defmethod (setf matrix-column) (col-table (mat <sparse-dictionary-matrix>) key)
  (setf (dic-ref (column-table mat) key) col-table))

(defmethod matrix-block ((smat <sparse-dictionary-matrix>) row-key col-key)
  (dic-ref (matrix-row smat row-key) col-key))

(defmethod row-keys ((mat <sparse-dictionary-matrix>))
  (keys (row-table mat)))

(defmethod col-keys ((mat <sparse-dictionary-matrix>))
  (keys (column-table mat)))

;;; derived methods

(defmethod for-each-entry-in-row (fn (smat <sparse-dictionary-matrix>) key)
  (dic-for-each-value fn (matrix-row smat key)))

(defmethod for-each-entry-in-col (fn (smat <sparse-dictionary-matrix>) key)
  (dic-for-each-value fn (matrix-column smat key)))

(defmethod for-each-key-in-row (fn (smat <sparse-dictionary-matrix>) key)
  (dic-for-each-key fn (matrix-row smat key)))

(defmethod for-each-key-in-col (fn (smat <sparse-dictionary-matrix>) key)
  (dic-for-each-key fn (matrix-column smat key)))

(defmethod for-each-key-and-entry-in-row (fn (smat <sparse-dictionary-matrix>) key)
  (dic-for-each fn (matrix-row smat key)))

(defmethod for-each-key-and-entry-in-col (fn (smat <sparse-dictionary-matrix>) key)
  (dic-for-each fn (matrix-column smat key)))

(defmethod for-each-row-key (fn (smat <sparse-dictionary-matrix>))
  "Loop through row keys."
  (assert (access-type smat :require :row))
  (dic-for-each-key fn (row-table smat)))

(defmethod parallel-for-each-row-key (fn (smat <sparse-dictionary-matrix>))
  (with-workers (fn)
    (dic-for-each-key #'work-on (row-table smat))))

(defmethod for-each-col-key (fn (smat <sparse-dictionary-matrix>))
  "Loop through column keys."
  (assert (access-type smat :require :column))
  (dic-for-each-key fn (column-table smat)))

(defmethod remove-key ((smat <sparse-dictionary-matrix>) &rest indices)
  (ecase (length indices)
    (2 (destructuring-bind (row-key col-key) indices
         (let ((row (matrix-row smat row-key))
               (column (matrix-column smat col-key)))
           (when row
             (dic-remove col-key row)
             (when (dic-empty-p row)
               (dic-remove row-key (row-table smat)))
             (dic-remove row-key column)
             (when (dic-empty-p column)
               (dic-remove col-key (column-table smat)))))))
    (1 (destructuring-bind (key) indices
         (loop for ck in (keys-of-row smat key)
            do (remove-key smat key ck))
         (loop for rk in (keys-of-column smat key)
            do (remove-key smat rk key))))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; hash-table-based sparse-matrix
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass <ht-sparse-matrix>
    (fl.multiprocessing:locked-region-mixin <sparse-dictionary-matrix>)
  ((row-table :accessor row-table :initarg :row-table
	      :initform (make-hash-table) :type hash-table
	      :documentation "Table of rows.")
   (column-table :accessor column-table :initarg :column-table
		 :initform (make-hash-table) :type hash-table
		 :documentation "Table of columns."))
  (:documentation "The <ht-sparse-matrix> represents an unordered matrix graph
indexed by general keys."))

(defun make-sparse-matrix
    (&rest args
     &key (type '(<block-definition-mixin> <ht-sparse-matrix>))
     &allow-other-keys)
  (apply #'fl.amop::make-programmatic-instance type
                        (sans args '(:type))))

(defmethod (setf matrix-block) (value (smat <ht-sparse-matrix>) row-key col-key)
  (setf (gethash col-key
		 (or (gethash row-key (row-table smat))
		     (setf (gethash row-key (row-table smat))
			   (make-hash-table :size 30))))
	value)
  (setf (gethash row-key
		 (or (gethash col-key (column-table smat))
		     (setf (gethash col-key (column-table smat))
			   (make-hash-table :size 30))))
	value))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Tests
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun test-sparse-matrix ()
  (flet ((constantly-1 (x) (declare (ignore x)) 1)
	 (constantly-2 (x) (declare (ignore x)) 2))
    (let* ((x (make-sparse-vector :key->size #'constantly-1))
	   (A (make-sparse-matrix
	       :row-key->size #'constantly-1 :col-key->size #'constantly-1
	       :keys->pattern (constantly (full-crs-pattern 1 1))))
	   (B (make-analog A))
	   (AA (make-sparse-matrix
		:row-key->size #'constantly-2 :col-key->size #'constantly-2
		:keys->pattern (constantly (full-crs-pattern 2 2)))))

      (setf (mref AA 0 1) #m((1.0 2.0e-15) (3.0 -4.0)))
      (assert (mzerop A))
      (setf (vref x 1) #m((1.0)))
      ;;
      (setf (mref B 0 0) #m((1.0)))
      (terpri)
      (show B)
      (show x)
      (show (m* B x))
      (sparse-matrix->matlisp AA :row-keys '(0) :col-keys '(1))
  
      (x<-0 A)
      (show B)
      (show (copy! B A))
      (assert (= (vref (mref A 0 0) 0) 1.0))
      ;;
      (sparse-matrix->matlisp A :keys '(0 1 2))
      (setf (mref AA 0 1) (make-full-crs-matrix 2 2))
      (fill! (mref AA 0 1) 1.0)
      ;;(sparse-matrix->matlisp AA :row-keys '(0 1) :col-keys '(1))
      ;;
      (setf (mref A 0 1) #m((2.0)))
      (show (sparse-m* A A))
      (show (sparse-m* A A :job :nt :sparsity :B))
      (show (sparse-m* A A :job :tn :sparsity :A))
      (show (sparse-m* A A :job :tt))
      (extract-value-blocks A '(0 1) '(0 1))
      (nr-of-entries A)
      ;; end of testing environment
      ))
  ;; sparse matrices with crs entries
  (let* ((offdiag (make-instance 'compressed-pattern
				 :sizes #(2 2) :orientation :row
				 :pattern '( ((a . 0))  ((a . 1)) )))
	 (diag (full-crs-pattern 2 2))
	 (n 3)
	 (key->size (constantly 2))
	 (A (make-sparse-matrix
	     :row-key->size key->size :col-key->size key->size
	     :keys->pattern (lambda (i j)
			      (if (eql i j) diag offdiag)))))
    ;; generate a 1d discretization
    (loop for i from 1 below n do
	  (setf (mref A i i) (diag (double-vec 2.0 2.0)))
	  (when (> i 1)
	    (fill! (mref A i (1- i)) -1.0))
	  (when (< i (- n 1))
	    (fill! (mref A i (1+ i)) -1.0)))
    (describe (mref A 1 2)))
  ;; end of test procedure
  )

;;; (test-sparse-matrix)
(fl.tests:adjoin-test 'test-sparse-matrix)

