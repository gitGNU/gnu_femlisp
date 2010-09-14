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
	   (base-type (if (atom element-type) element-type (cdr element-type)))
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

  ;; BUG: no error is reported, but the result is bad
  (let ((A #M((6.2172489379008845e-15 4.440892098500626e-15 2.3822801641527197e-22 -1.4210854715202004e-14 -1.4210854715202004e-14 -3.552713678800501e-15)
	      (7.774428950278327e-16 10.08511217749292 1.9497552341246447e-7 -4.369011214675872 -5.623824671657317 0.3783976796049302)
	      (5.955700410381799e-23 1.949755234124645e-7 1.0298405279639302e-14 -8.244351100690419e-8 7.24674978130909e-7 5.842130154715112e-7)
	      (-8.415728594088137e-16 -4.369011214675867 -8.244351100690379e-8 102.1018806730008 131.42632218993467 -8.842988225843465)
	      (-1.9949351081100997e-14 -5.623824671657299 7.246749781309091e-7 131.42632218993467 317.16159112054163 -8.838690372574526)
	      (-9.212035778385909e-15 0.37839767960492665 5.84213015471511e-7 -8.842988225843474 -8.838690372574511 170.65141602931885)))
	(B #M((1.0 -1.6653345369377348e-16 2.897448009914831e-8 4.107825191113079e-15 8.715250743307479e-15 -1.942890293094024e-15)
	      (-2.7755575615628914e-17 0.9999999999999997 1.9341748841415267e-8 -7.53563877964325e-15 -1.1657341758564144e-14 7.549516567451064e-15)
	      (2.897448009914831e-8 1.9341748841415264e-8 1.3023483917789837e-15 -6.8457868359323e-9 5.5569597313562905e-9 3.3136154966616022e-9)
	      (4.163336342344337e-15 -7.632783294297951e-15 -6.845786835932301e-9 1.0000000000003408 4.4794723486063504e-13 -1.307565167252278e-13)
	      (8.659739592076221e-15 -1.1823875212257917e-14 5.556959731356288e-9 4.4775294583132563e-13 1.000000000000597 -1.6767143229401427e-13)
	      (-2.0539125955565396e-15 7.66053886991358e-15 3.3136154966616055e-9 -1.3078427230084344e-13 -1.6765755450620645e-13 1.000000000000044))))
    (hegv A B :V))

  )

#-lispworks
(fl.tests:adjoin-test 'test-hegv)
;;; (test-hegv)
