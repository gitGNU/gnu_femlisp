;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; hegv.lisp - Hermitian generalized eigenvalue problem
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;
;;; Copyright (C) 2007 Nicolas Neuss, University of Karlsruhe.
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
;;; NO EVENT SHALL THE AUTHOR, THE UNIVERSITY OF KARLSRUHE, OR OTHER
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

(defgeneric hegv (a b &optional job)
  (:documentation
   "Syntax: (HEGV A B [job])

Purpose: Computes the generalized eigenvalues and left/right eigenvectors
of @math{A - s B} for Hermitian matrices A and B.

  1. (HEGV A B :N) => lambda

     Computes the generalized eigenvalues of @math{A - s B}.

  2. (HEGV A B :V) => lambda, V

     Computes generalized eigenvalues and eigenvectors of (A - sB).

     @math{ A*V = B*V*diag(lambda), \\ W'*A = diag(lambda)*W'*B}
            
     with V and W orthogonal (unitary).

Remark: The non-symmetric counterpart of this routine is
@function{ggev}."))

(defun symmetric-packed-store (A &optional (part :lower))
  (assert (= (nrows A) (ncols A)))
  (lret* ((n (nrows A))
	  (store (zero-vector (/ (* n (1+ n)) 2) (element-type A))))
    (ecase part
      (:lower
       (loop for i below n do
	    (loop for j upto i do
		 (setf (vref store
                             ;; #I"i+j*(2*n-j-1)/2"
                             (+ i (* 1/2 j (- (* 2 n) j 1))))
		       (mref A i j)))))
      (:upper
       (loop for j below n do
	    (loop for i upto j do
		 (setf (vref store
                             ;; #I"i+j*(j+1)/2"
                             (+ i (* 1/2 j (1+ j))))
		       (mref A i j))))))))

(defmethod hegv ((A standard-matrix) (B standard-matrix) &optional (job :N))
  (let ((n (nrows A))
	(element-type (element-type A)))
    (unless (= n (ncols A) (nrows B) (ncols B))
      (error "Arguments are not square matrices of the same size"))
    (unless (equal element-type (element-type B))
      (error "Element types do not match"))
    (let* ((A (copy A))
	   (B (copy B))
	   (itype 1)
	   (etype (cl->lapack-type element-type))
	   (real-p (member etype '(:float :double)))
	   (base-type (if (atom element-type) element-type (second element-type)))
	   (w (zeros n 1 base-type))
	   (v-p (eq job :V))
	   (Z (when real-p (zeros n n element-type)))
	   (n8 (* 8 n))
	   (work (zero-vector n8 element-type))
	   (rwork (unless real-p (zero-vector n8 base-type))))
      (call-lapack-with-error-test
       (lapack fl.lapack::*hegv* etype)
       itype (symbol-name job) "U" n
       (if real-p (symmetric-packed-store A :upper) A) (unless real-p n)
       (if real-p (symmetric-packed-store B :upper) B) (unless real-p n)
       w Z (when real-p n)
       work n8 rwork 0)
      (values w (and v-p (if real-p Z A))))))

(defun test-hegv ()
  (fl.lapack::remove-all-lapack-functions fl.lapack::*hegv*)
  (symmetric-packed-store (eye 3) :upper)
  (let* ((n 200)
	 (A (mrandom n n))
	 (B (eye n)))
    (lapack fl.lapack::*hegv* (cl->lapack-type (element-type B)))
    (time (hegv A B :V))
    nil)
  (hegv (eye 2 2 '(complex double-float)) (eye 2 2 '(complex double-float)))
  (hegv (ones 2 2 '(complex double-float)) (eye 2 2 '(complex double-float)) :V)
  (let ()
    (hegv (ones 2 2 '(complex double-float)) (eye 2 2 '(complex double-float)) :V)
    (hegv (ones 2 2 '(complex double-float)) (eye 2 2 '(complex double-float)) :N)
    (hegv #m((1.0 1.0) (1.0 2.0)) (eye 2) :V)
    (hegv #m((1.0 1.0) (1.0 2.0)) (eye 2) :N))
  (hegv #m((1.0 1.0) (1.0 2.0)) (eye 2) :N)
  (nth-value 1 (hegv #m((1.0 1.0) (1.0 2.0)) (eye 2) :V))
  (hegv #m((0.0 1.0) (1.0 0.0)) (eye 2) :V)
  (hegv (eye 2) (eye 2) :V)

  ;; BUG: no error is reported for this bad input
  (let ((A (eye 2))
	(B #M((1.0 0.0) (0.0 0.0))))
    (hegv A B :V))

  )

#-lispworks
(fl.tests:adjoin-test 'test-hegv)
;;; (test-hegv)
