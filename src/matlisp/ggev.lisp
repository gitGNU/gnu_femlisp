;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; ggev.lisp - Generalized eigenvalue problem
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

(defgeneric ggev (a b &optional job)
  (:documentation
   "Syntax: (GGEV A B [job])

Purpose: Computes the generalized eigenvalues and left/right eigenvectors
of @math{A - s B}.

  1. (GGEV A B :N) => lambda

     Computes the generalized eigenvalues of @math{A - s B}.

  2. (GGEV A B :V) => lambda, V, W

     Computes generalized eigenvalues and eigenvectors of (A - sB).

     @math{ A*V = B*V*diag(lambda), \\ W'*A = diag(lambda)*W'*B}
            
     with V and W orthogonal (unitary).

Remark: The symmetric/hermitian counterpart of this routine is
@function{hegv}."))


(defun combine-real-and-imag-part (rv iv)
  (if (mzerop iv)
      rv
      (lret* ((m (nrows rv))
	      (n (ncols rv))
	      (result (zeros m n `(complex ,(element-type rv)))))
	(dotimes (i m)
	  (dotimes (j n)
	    (setf (mref result i j)
		  (complex (mref rv i j) (mref iv i j))))))))

(defmethod ggev ((A standard-matrix) (B standard-matrix) &optional (job :N))
  (let* ((n (nrows A))
	 (n8 (* 8 n))
	 (element-type (element-type A))
	 (A (copy A))
	 (B (copy B)))
    (unless (= n (ncols A) (nrows B) (ncols B))
      (error "Arguments are not square matrices of the same size"))
    (unless (equal element-type (element-type B))
      (error "Element types do not match"))
    (let* ((alpha (zeros n 1 element-type))
	   (beta (zeros n 1 element-type))
	   (etype (cl->lapack-type element-type))
	   (real-p (member etype '(:float :double)))
	   (v-p (eq job :V))
	   (VL (and v-p (zeros n n element-type)))
	   (VR (and v-p (zeros n n element-type)))
	   (work (zero-vector n8 element-type))
	   (alpha-i (when real-p (zeros n 1 element-type)))
	   (rwork (unless real-p (zero-vector n8 (cdr element-type)))))
      (call-lapack-with-error-test
       (lapack "ggev" etype) (symbol-name job) (symbol-name job)
       n A n B n alpha alpha-i beta
       (or VL A) n (or VR B) n
       work n8 rwork 0)
      (values (if real-p
		  (combine-real-and-imag-part alpha alpha-i)
		  alpha)
	      VL VR))))

(defun test-ggev ()
  (let* ((n 200)
	 (A (mrandom n n))
	 (B (eye n)))
    (time (ggev A B :V))
    nil)
  (ggev (ones 2 2 '(complex double-float)) (eye 2 2 '(complex double-float)))
  (ggev #m((1.0 1.0) (0.0 2.0)) (eye 2) :N)
  (ggev #m((0.0 1.0) (-1.0 0.0)) (eye 2) :V)
  (ggev (eye 2) (eye 2) :V)
  )

(fl.tests:adjoin-test 'test-ggev)
;;; (test-ggev)
