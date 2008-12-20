;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; standard-matrix-lr.lisp - LAPACK routines getrf/getrs
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

;;; Helpers

(declaim (inline swap-rows))
(defun swap-rows (j k mat-store m n)
  "Swaps rows j and k in an m x n-matrix."
  (loop for pos1 of-type fixnum
	from (indexing j 0 m n) below (indexing j n m n) by m
	and pos2 of-type fixnum
	from (indexing k 0 m n) by m do
	(rotatef (aref mat-store pos1) (aref mat-store pos2))))

(declaim (inline trsm))
(defun trsm (lr-store b-store m n &key side unit-diagonal)
  (macrolet ((b-index (i j) `(indexing ,i ,j m n))
	     (lr-index (i j) `(indexing ,i ,j m m))
	     (b-ref (i) `(aref b-store ,i))
	     (lr-ref (i) `(aref lr-store ,i)))
    (dotimes (j n)
      (declare (fixnum j))
      (ecase side
	(:lower
	 (loop for b-off-1 of-type fixnum from (b-index 0 j) below (b-index m j)
	       and lr-off-1 of-type fixnum from (lr-index 0 0) by (1+ m) do
	       (unless unit-diagonal
		 (setf (b-ref b-off-1)
		       (/ (b-ref b-off-1) (lr-ref lr-off-1))))
	       (loop for b-off-2 of-type fixnum from (+ b-off-1 1) below (b-index m j)
		     and lr-off-2 of-type fixnum upfrom (+ lr-off-1 1) do
		     (decf (b-ref b-off-2)
			   (* (b-ref b-off-1) (lr-ref lr-off-2))))))
	(:upper
	 (loop for b-off-1 of-type fixnum from (b-index (1- m) j) downto (b-index 0 j)
	       and lr-off-1 of-type fixnum downfrom (lr-index (1- m) (1- m)) by (1+ m) do
	       (unless unit-diagonal
		 (setf (b-ref b-off-1)
		       (/ (b-ref b-off-1) (lr-ref lr-off-1))))
	       (loop for b-off-2 of-type fixnum from (- b-off-1 1) downto (b-index 0 j)
		     and lr-off-2 of-type fixnum downfrom (- lr-off-1 1) do
		     (decf (b-ref b-off-2)
			   (* (b-ref b-off-1) (lr-ref lr-off-2))))))))))


(define-blas-template getrf! ((mat standard-matrix) &optional ipiv)
  "Computes the PLR decomposition with column pivoting of the matrix MAT."
  (let* ((m (slot-value mat 'nrows))
	 (n (slot-value mat 'ncols))
	 (k (min m n))
	 (ipiv (cond ((null ipiv)
		      (make-array k :element-type 'fixnum :initial-element -1))
		     ((eq ipiv :none) nil)
		     (t ipiv))))
    (declare (type fixnum m n k)
	     (type (or null (simple-array fixnum (*))) ipiv))
    (with-blas-data (mat)
      (loop
	 for j of-type fixnum from 0 below k
	 for offset-jj of-type fixnum = (indexing j j m n) do
	 ;; column pivoting
	 (when ipiv
	   (let ((pivot-index j)
		 (pivot-value (abs (aref mat-store offset-jj))))
	     (declare (type fixnum pivot-index))
	     ;; find pivot
	     (loop for i from (1+ j) below m
		for off = (indexing i j m n)
		for value of-type element-type = (abs (aref mat-store off))
		do (when (> value pivot-value)
		     (setq pivot-index i
			   pivot-value value)))
	     (setf (aref ipiv j) pivot-index)
	     (unless (= pivot-index j)
	       (swap-rows j pivot-index mat-store m n))))
	 ;;
	 (when (zerop (aref mat-store offset-jj))
	   (return-from getrf! (values mat ipiv j)))
	   
	 ;; compute elements of column J of L
	 (loop with factor = (/ (aref mat-store offset-jj))
	    for off of-type fixnum
	    from (1+ offset-jj) below (indexing n j m n) do
	    (setf (aref mat-store off) (* (aref mat-store off) factor)))
	   
	 ;; update trailing submatrix of R
	 (loop
	    for off of-type fixnum
	    from (indexing j (1+ j) m n) below (indexing j n m n) by m
	    for factor of-type element-type = (aref mat-store off)
	    unless (zerop factor) do
	    (loop for pos1 of-type fixnum
	       from (1+ off) below (+ off (- n j))
	       and pos2 of-type fixnum from (1+ offset-jj) do
	       (decf (aref mat-store pos1)
		     (* factor (aref mat-store pos2)))))))
    (values mat ipiv t)))

(define-blas-template getrs! ((LR standard-matrix) (b standard-matrix) &optional ipiv)
  "Uses the LR decomposition computed by getrf! to solve a linear system
with rhs B.  LR must be a n x n - matrix, b must be a n x m matrix."
  (with-blas-data (LR b)
    (unless (= b-nrows LR-nrows LR-ncols)
      (error "Matrix LR is not quadratic or does not fit to right-hand side."))
    (when ipiv    ; swap rhs according to ipiv
      (locally
	  (declare (type (simple-array fixnum (*)) ipiv))
	(unless (= b-nrows (length ipiv))
	  (error "Matrix and pivot vector do not fit."))
	(dotimes (i b-nrows)
	  (declare (type fixnum i))
	  (let ((k (aref ipiv i)))
	    (unless (= i k)
	      (swap-rows i k b-store b-nrows b-ncols))))))
    
    ;; solve Lx'=b
    (trsm LR-store b-store b-nrows b-ncols :side :lower :unit-diagonal t)
    ;; solve Rx=x'
    (trsm LR-store b-store b-nrows b-ncols :side :upper :unit-diagonal nil))
  b)

#+(or)
(multiple-value-bind (lr ipiv)
    (getrf! #m((1.0 1.0) (1.0 0.0)))
  (getrs! lr (eye 2) ipiv))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; derived routines (special for standard-matrix because of eye)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(definline m/ (x)
  "Returns the inverse of X."
  (if (numberp x)
      (/ x)
      (gesv! x (eye (nrows x) (ncols x) (element-type x)))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; unoptimized
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defparameter *cholesky-threshold* 1.0e-8)

(defun cholesky (A &key inverse compactify (threshold *cholesky-threshold*))
  (let* ((n (nrows A))
	 (G (zeros n)))
    (loop for j below n do
	 (let ((v (matrix-slice A :from-row j :from-col j :nrows 1)))
	   (dotimes (k j)
	     (axpy! (- (mref G k j))
		    (matrix-slice G :from-row k :from-col j :nrows 1)
		    v))
	   (let ((s (vref v 0)))
	     (when (minusp s) (setf s (- s)))
	     (assert (not (minusp s)))
	     (when (> s threshold)
	       (scal! (/ (sqrt s)) v)
	       (minject! v G j j)))))
    (when inverse
      ;; invert G
      (let ((F (zeros n n)))
	(loop for i below n
	   do
	   (when (plusp (mref G i i))
	     (assert (not (minusp (mref G i i))))
	     (setf (mref F i i) (/ (mref G i i)))
	     (loop for j from (1- i) downto 0 do
		  (setf (mref F j i)
			;; use 0=delta_ji=sum_{k=j}^i F_jk G_ki, 
			(/ (- (loop for k from j below i
				 sum (* (mref F j k) (mref G k i))))
			   (mref G i i))))))
	(setf G F)))
    ;; compactify G by skipping columns with zeros on the diagonals
    (when compactify
      (let* ((m (loop for i below n count (not (zerop (mref G i i)))))
	     (F (zeros n m))
	     (i 0))
	(dotimes (j n)
	  (unless (zerop (mref G j j))
	    (dotimes (k n)
	      (setf (mref F k i) (mref G k j)))
	    (incf i)))
	(setf G F)))
    ;; return the result
    G))



;;;; Testing

(defun test-standard-matrix-lr ()
  (multiple-value-bind (lr pivot)
      (getrf! (scal! 0.5 (eye 1)))
    (getrs! lr #m((1.0 2.0)) pivot))

  (multiple-value-bind (lr pivot)
      (getrf! #m((0.0 -1.0) (1.0 0.0)))
    (getrs! lr #m((1.0) (2.0)) pivot))

  (let ((A #m((4.0 -1.0) (-1.0 4.0))))
    (assert (mzerop (m- (m* (m/ A) A) (eye 2)))))
  
  #+(or)
  (time
   (let* ((mat (eye 3))
	  (mat2 (copy mat)))
     (loop repeat 100000 do
	   (copy! mat mat2)
	   (getrf! mat2))))
  #+(or)
  (time
   (let* ((mat (matlisp::eye 1))
	  (mat2 (matlisp::copy mat)))
     (loop repeat 10000 do
	   (matlisp::copy! mat mat2)
	 (matlisp::getrf! mat2))))
  
  ;; we can handle also general matrices, however we comment this out at
  ;; the moment, due to an optimization warning of CMUCL/SBCL which we do
  ;; not want to see when testing -- allowed again 11.12.2008
  (mequalp (m/ (make-instance (standard-matrix t)
			      :content '((1 2) (3 4))))
	   (make-instance (standard-matrix t)
			  :content '((-2 1) (3/2 -1/2))))
  )

;;; (fl.matlisp::test-standard-matrix-lr)
(fl.tests:adjoin-test 'test-standard-matrix-lr)
