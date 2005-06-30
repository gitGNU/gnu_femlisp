;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; number-blas.lisp - BLAS operations for numbers
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

;;; Ensure some BLAS operations for numbers

(in-package :fl.matlisp)

(defmethod nrows ((obj number)) 1)
(defmethod ncols ((obj number)) 1)

(defmethod total-entries ((obj number)) 1)

(defmethod copy ((x number)) x)

(defmethod m* ((x number) (y number))
  "Ordinary multiplication on numbers."
  (* x y))
(defmethod m+ ((x number) (y number))
  "Number addition."
  (+ x y))

(defmethod copy! ((x number) (y number))
  "Number copying."
  x)
(defmethod m+! ((x number) (y number))
  "Number addition."
  (+ x y))
(defmethod m*! ((x number) (y number))
  "Number multiplication."
  (* x y))

(defmethod fill-random! ((x number) (s number))
  (random s))

(defmethod scal! (s (x number)) (* s x))

(defmethod axpy! (alpha (x number) (y number))
  (+ (* alpha x) y))

(defmethod l2-norm ((x number)) (abs x))
(defmethod lp-norm ((x number) p)
  (declare (ignore p))
  (abs x))
(defmethod linf-norm ((x number)) (abs x))

(defmethod dot ((x number) (y number)) (* x y))

(defmethod dot-abs ((x number) (y number)) (abs (* x y)))

(defmethod mequalp ((x number) (y number)) (= x y))

(defmethod mzerop ((x number) &optional (threshold 0.0))
  "Tests if the number @arg{x} is lower or equal to @arg{threshold}."
  (<= (abs x) threshold))

(defmethod getrf! ((x number) &optional ipiv)
  (assert (null ipiv))
  (if (zerop x)
      (values 0 nil 0)
      (values x nil t)))

(defmethod getrs! ((x number) (b number) &optional ipiv)
  (assert (null ipiv))
  (/ b x))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Testing
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun test-number-blas ()
  (mzerop (m- (scal 0.5 1.0) (axpy 0.5 1.0 0.0)))
  )

;;; (test-number-blas)
(fl.tests:adjoin-test 'test-number-blas)


