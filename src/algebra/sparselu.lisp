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

(in-package :algebra)

(defclass <ldu-sparse> ()
  ((lower-left :reader lower-left :initarg :lower-left :type <sparse-matrix>)
   (diagonal :reader diagonal :initarg :diagonal :type <sparse-matrix>)
   (upper-right :reader upper-right :initarg :upper-right :type <sparse-matrix>)
   (ordering :reader ordering :initarg :ordering :type vector)))

(defmethod getrf! ((A <sparse-matrix>) &optional ipiv)
  (declare (ignore ipiv))
  (values (sparse-ldu A) nil t))


(defun modify-diagonal (D R omega)
  "Implements the diagonal modification for ILU_mod."
  (dorows (rk R)
    (let ((sum 0.0d0))
      (for-each-key-in-row
       #'(lambda (ck)
	   (incf sum (mat-ref R rk ck)))
       R rk)
      (decf (mat-ref D rk rk) (* omega sum))))
  D)

(defun shift-diagonal-inverter (eta)
  "Can be used for obtainint a diagonal modification to get ILU_mod."
  #'(lambda (mat)
      (dorows (i mat)
	(incf (mat-ref mat i i) eta))
      (m/ mat)))

(defmethod sparse-ldu ((A <sparse-matrix>) &key ordering incomplete (omega 0.0d0) diagonal-inverter)
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
	  for entry = (mat-ref A row-key row-key) do
	  (assert entry)
	  (setf (mat-ref U row-key row-key)
		(typecase entry
		  (<crs-matrix> (crs->matlisp entry))
		  (t (copy entry)))))
    (loop for row-key across ordering do
	  (for-each-key-and-entry-in-row
	   #'(lambda (col-key entry)
	       (when (and (matrix-row U col-key) (not (mzerop entry)))
		 (setf (mat-ref U row-key col-key)
		       (typecase entry
			 (<crs-matrix> (crs->matlisp entry))
			 (t (copy entry))))))
	   A row-key))

    ;; decomposition
    (loop
     for k across ordering do
     (let* ((U_kk (mat-ref U k k))
	    (D_kk (if diagonal-inverter
		      (funcall diagonal-inverter U_kk)
		      (m/ U_kk))))
     (setf (mat-ref D k k) D_kk)	; store D
     (remove-entry U k k)		; clean U
     (let ((col-k (matrix-column U k))
	   (row-k (matrix-row U k)))
       (when col-k
	 (loop for i being each hash-key of col-k
	       and Uik being each hash-value of col-k
	       unless (matrix-row D i) do
	       (let ((factor (m* Uik D_kk))
		     (row-i (matrix-row U i)))
		 (setf (mat-ref L i k) factor) ; store L
		 (remove-entry U i k)	; clean U
		 (when row-k
		   (loop for j being each hash-key of row-k
			 and Ukj being each hash-value of row-k do
			 (let ((mblock (gethash j row-i)))
			   (cond ((or mblock (not incomplete))
				  (setq mblock (or mblock (mat-ref U i j)))
				  ;; (decf (vec-ref mblock 0) (* (vec-ref factor 0) (vec-ref Ukj 0)))
				  (gemm! -1.0d0 factor Ukj 1.0d0 mblock))
				 (t (unless (zerop omega)
				      (modify-diagonal (mat-ref U i i) (m* factor (mat-ref U i j)) omega)))
				 ))))))))))
    ;; finally return the result
    (make-instance '<ldu-sparse> :lower-left L :diagonal D
		   :upper-right U :ordering ordering)))


;;; Matlisp/LAPACK interface

(defmethod getrf! ((A <sparse-matrix>) &optional ipiv)
  (declare (ignore ipiv))
  (values (sparse-ldu A :ordering nil) nil t))

(defmethod getrs! ((ldu <ldu-sparse>) ipiv (result <sparse-vector>) &key trans)
  (declare (ignore trans))
  
  ;; solve (I + L) result~ = rhs
  (loop with L of-type <sparse-matrix> = (lower-left ldu)
	for i across (ordering ldu)
	for result-i = (vec-ref result i)
	when (matrix-row L i) do
	(for-each-key-and-entry-in-row
	 #'(lambda (j Lij) (x-=Ay result-i Lij (vec-ref result j)))
	 L i))
  
  ;; solve (D + U) result = result~ (where D^{-1} is stored)
  (loop with U of-type <sparse-matrix> = (upper-right ldu)
	with D of-type <sparse-matrix> = (diagonal ldu)
	for i across (reverse (ordering ldu))
	for result-i = (vec-ref result i) do
	(when (matrix-row U i)
	  (for-each-key-and-entry-in-row
	   #'(lambda (j Uij) (x-=Ay result-i Uij (vec-ref result j)))
	   U i))
	(m*! (mat-ref D i i) result-i))
  
  ;; return the result
  (values result ipiv t))

(defmethod getrs ((ldu <ldu-sparse>) ipiv (rhs <sparse-vector>) &key trans)
  (declare (ignore ipiv trans))
  (let ((result (copy rhs)))
    (getrs! ldu nil result)))


;;; Alternative interface

(defmethod x<-Ay ((result <sparse-vector>) (ldu <ldu-sparse>) (rhs <sparse-vector>))
  "Performs a matrix multiplication with U^-1 D^-1 L^-1."
  (declare (optimize (speed 3) (safety 0)))
  
  ;; copy rhs to result, further we work only on result
  (copy! rhs result)
  (getrs! ldu nil result))

(defmethod m* ((ldu <ldu-sparse>) (rhs <sparse-vector>))
  (declare (optimize (speed 3) (safety 0)))
  (let ((result (copy rhs)))
    (x<-Ay result ldu rhs)))


;;; Testing:

(defun test-sparselu ()
  
  (let ((A (make-sparse-matrix
	    :row-key->size (constantly 1) :col-key->size (constantly 1)
	    :keys->pattern (constantly (full-crs-pattern 1 1)))))
    (describe A)
    (setf (mat-ref A 0 0) [[2.0]])
    (setf (mat-ref A 1 0) [[1.0]])
    (setf (mat-ref A 1 1) [[1.0]])
    (setf (mat-ref A 0 1) [[1.0]])
    
    (display A)
    (display (lower-left (sparse-ldu A)))
    (display (diagonal (sparse-ldu A)))
    (display (upper-right (sparse-ldu A)))
    (ordering (sparse-ldu A))

    (let ((rhs (make-instance '<sparse-vector> :key->size (constantly 1))))
      (setf (vec-ref rhs 0) [3.0])
      (setf (vec-ref rhs 1) [1.0])
      (show rhs)
      (show (m* (sparse-ldu A) rhs))))

  (let ((A (make-sparse-matrix
	    :row-key->size (constantly 2) :col-key->size (constantly 2)
	    :keys->pattern (constantly (full-crs-pattern 2 2)))))
    (describe A)
    (setf (mat-ref A 0 0) [[4.0 -1.0]' [-1.0 4.0]'])
    (setf (mat-ref A 1 0) [[-1.0 0.0]' [0.0 -1.0]'])
    (setf (mat-ref A 1 1) [[4.0 -1.0]' [-1.0 4.0]'])
    (setf (mat-ref A 0 1) [[-1.0 0.0]' [0.0 -1.0]'])
    
    (display A)
    (display (lower-left (sparse-ldu A)))
    (display (diagonal (sparse-ldu A)))
    (display (upper-right (sparse-ldu A)))
    (ordering (sparse-ldu A))
    (let ((rhs (make-instance '<sparse-vector> :key->size (constantly 2))))
      (setf (vec-ref rhs 0) [1.0 2.0]')
      (setf (vec-ref rhs 1) [3.0 4.0]')
      (show rhs)
      (show (m* (sparse-ldu A) rhs))
      (show (m* A (m* (sparse-ldu A) rhs)))))
  )