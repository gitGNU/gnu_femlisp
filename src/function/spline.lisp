;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; spline.lisp - Cubic splines (following Stoer/Bulirsch)
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

(in-package :fl.function)

(defun apply-1d-stencil (stencil x)
  (let* ((n (length x))
	 (y (make-array n)))
    (dotimes (i n y)
      (setf (aref y i)
	    (+ (* (aref stencil 0) (aref x (mod (1- i) n)))
	       (* (aref stencil 1) (aref x i))
	       (* (aref stencil 2) (aref x (mod (1+ i) n))))))))
    
(defun M-times (x)
  (apply-1d-stencil #(1.0 4.0 1.0) x))

(defun spline-right-hand-side (x)
  (let ((h (/ 1.0 (length x))))
    (apply-1d-stencil (scal (/ (* h h)) #(6.0 -12.0 6.0)) x)))

(defun solve-moment-system (b &optional (n 50))
  "Solves moment system by Jacobi iteration."
  (loop with x = (scal 0.25 b)
	repeat n do
	(setq x (m+ x (scal 0.25 (m- b (M-times x)))))
	finally (return x)))

#+(or)
(solve-moment-system #(10.0 30.0 12.0) 60)

(defun cubic-spline (y)
  "On a regular partition of the unit interval interpolating values y are
given.  This function returns an interpolating spline."
  (let* ((n (length y)) (h (/ 1.0 n))
	 (d (spline-right-hand-side y))
	 (M (solve-moment-system d))
	 (alpha y)
	 (gamma (scal 0.5 M))
	 (beta (m+ (apply-1d-stencil (scal (/ h) #(0.0 -1.0 1.0)) y)
		   (apply-1d-stencil (scal (/ h 6) #(0.0 -2.0 -1.0)) M)))
	 (delta (apply-1d-stencil (scal (/ 1.0 h 6) #(0.0 -1.0 1.0)) M)))
    (values
     ;; s(f)
     #'(lambda (x)
	 (multiple-value-bind (k xi)
	     (floor (aref x 0) h)
	   (setq k (mod k n))
	   #I"alpha[k]+xi*(beta[k]+xi*(gamma[k]+xi*delta[k]))"))
     ;; s(f)'
     #'(lambda (x)
	 (multiple-value-bind (k xi)
	     (floor (aref x 0) h)
	   (setq k (mod k n))
	   (vector
	    #I"beta[k]+xi*(2.0*gamma[k]+3.0*xi*delta[k])"))))))
     

#+(or)  ; needs package "FL.PLOT"!
(flet ((f (x) #I"sin(2*pi*x)"))
  (let* ((N 3)
	 (h (/ 1.0 N))
	 (plot-h (/ h 10))
	 (ip-coords (map 'vector (curry #'* h) (range< 0 N)))
	 (spline (cubic-spline (vector-map #'f ip-coords))))
    (fl.plot:plot
     (list
      (cons
       "function"
       (loop for x from 0.0 upto 1.0 by plot-h
	     collect (vector x (f x))))
      (cons
       "spline"
       (loop for x from 0.0 upto 1.0 by plot-h
	     collect (vector x (funcall spline (vector x)))))))))

