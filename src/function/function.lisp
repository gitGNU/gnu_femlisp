;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; function.lisp - Functions and functionals
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

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; function class
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass <function> ()
  ((domain-dimension :reader domain-dimension :initarg :domain-dimension)
   (image-dimension :reader image-dimension :initarg :image-dimension))
  (:documentation "The <function> class is an abstract class for a general
function.  This function will usually accept vector arguments, the dimensions of
domain and image are fixed when defining the function. If the function is
differentiable, the gradient matrix can be obtained by evaluating the gradient
slot."))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; functionals
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;(defclass <functional> ()
;  (functional :accessor my-functional :initarg :functional))
;
;(defmethod evaluate ((f' <functional>) f)
;  (funcall (functional f') f))



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Interface
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defgeneric evaluate (f x)
  (:documentation "Generic evaluation of functions on an argument.  Numbers and
arrays are treated as constants.  Special evaluation is defined for multivariate
polynomials on vectors and for <function> objects."))

(defgeneric evaluate-gradient (f x)
  (:documentation "Generic evaluation of gradients of differentiable functions."))

(defgeneric evaluate-k-jet (f k x)
  (:documentation "Generic evaluation of k-jets of C^k-functions."))

(defgeneric smoothness (f)
  (:documentation "Returns smoothness of function. :INFINITY is returned for
infinitely differentiable functions."))

(defgeneric differentiable-p (f &optional k)
  (:documentation "Returns t if f is differentiable or differentiable of the
given degree."))

;;; further generics
;(define-generic differentiate)
;(define-generic gradient)
;(define-generic integrate)

(defmethod evaluate-k-jet (f k x)
  "By default, the k-jet evaluation handles the cases K=0 and K=1 by calls
to EVALUATE and EVALUATE-GRADIENT, respectively."
  (cond ((= k 0) (evaluate f x))
	((= k 1) (evaluate-gradient f x))
	(t (call-next-method))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Self-evaluating objects
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmethod evaluate ((func function) x)
  "An ordinary procedure or generic is funcalled on x."
  (funcall func x))

(defmethod evaluate ((num number) x)
  "Numbers are treated as constant functions."
  (declare (ignore x))
  num)

(defmethod evaluate ((vec array) x)
  "Arrays are treated as constant functions.  This feature should probably not
be used, because an alternative interpretation would be to apply evaluate to all
entries, so misunderstandings may arise."
  (declare (ignore x))
  vec)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; <special-function>
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass <special-function> (<function>)
  ((evaluator :reader evaluator :initarg :evaluator :type (or function null))
   (gradient :reader gradient :initarg :gradient :initform nil :type (or function null))
   (jet :reader jet :initarg :jet :initform nil :type (or function null)))
  (:documentation "A <special-function> provides its own evaluation and gradient
computation."))

(defmethod evaluate ((f <special-function>) x)
  (funcall (evaluator f) x))

(defmethod differentiable-p ((f <special-function>) &optional (k 1))
  (if (= k 1)
      (gradient f)
      (jet f)))

(defmethod evaluate-gradient ((f <special-function>) x)
  (funcall (gradient f) x))

(defmethod evaluate-k-jet ((f <special-function>) k x)
  (funcall (jet f) k x))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; <constant-function>
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass <constant-function> (<function>)
  ((value :reader value :initarg :value :initarg :value))
  (:documentation "For a <constant-function> evaluation and derivative
computation are trivial."))

(defmethod evaluate ((f <constant-function>) x)
  (assert (eq (length x) (domain-dimension f)))
  (value f))

(defmethod differentiable-p ((f <constant-function>) &optional (k 1))
  (declare (ignore k))
  t)

(defmethod smoothness ((f <constant-function>)) :infinity)

(defmethod evaluate-gradient ((f <constant-function>) x)
  (assert (eq (length x) (domain-dimension f)))
  (make-real-matrix (image-dimension f) (domain-dimension f)))

(defun zero-k-jet (f k)
  (make-real-tensor
   (coerce (cons (image-dimension f)
		 (make-list k :initial-element (domain-dimension f)))
	   'fixnum-vec)))

(defmethod evaluate-k-jet ((f <constant-function>) k x)
  (if (<= k 1)
      (call-next-method)
      (zero-k-jet f k)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; <linear-function>
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass <linear-function> (<function>)
  ((A :initform nil :initarg :A)
   (b :initform nil :initarg :b))
  (:documentation "A <linear-function> is determined by a matrix A and a vector b."))

(defmethod initialize-instance :after ((f <linear-function>) &key &allow-other-keys)
  (with-slots (A b domain-dimension image-dimension) f
  (setf image-dimension (nrows A))
  (setf domain-dimension (ncols A))
  (unless b (setf b (make-double-vec (image-dimension f))))))

(defmethod differentiable-p ((f <linear-function>) &optional (k 1))
  (declare (ignore k))
  t)

(defmethod smoothness ((f <linear-function>)) :infinity)

(defmethod evaluate ((f <linear-function>) x)
  (with-slots (A b) f
    (gemm 1.0 A x 1.0 b)))

(defmethod evaluate-gradient ((f <linear-function>) x)
  (or (slot-value f 'A)
      (zeros (image-dimension f) (domain-dimension f))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; transformed functions
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass <linearly-transformed-function> (<function>)
  ((original :reader original :initarg :original)
   (domain-A :reader domain-A :initform nil :initarg :domain-A)
   (domain-b :reader domain-b :initform nil :initarg :domain-b)
   (image-A :reader image-A :initform nil :initarg :image-A)
   (image-b :reader image-b :initform nil :initarg :image-b))
  (:documentation "<linearly-transformed-function> calls the original
function evaluation on linearly transformed coordinates and linearly
transforms also the result."))

(defun transform-function (func &key domain-transform image-transform)
  (let ((domain-A (first domain-transform))
	(domain-b (second domain-transform))
	(image-A (first image-transform))
	(image-b (second image-transform)))
    (when domain-transform
      (assert (= (domain-dimension func) (length domain-b) (nrows domain-A))))
    (when image-transform
      (assert (= (image-dimension func) (ncols image-A))))
    ;; compute resulting values, we always expect pairs.  This could be made
    ;; faster by allowing shifts only.
    (when (typep func '<linearly-transformed-function>)
      (when (domain-A func)
	(setq domain-A (if domain-A
			   (m* (domain-A func) domain-A)
			   (domain-A func)))
	(setq domain-b (if domain-b
			   (gemm 1.0 (domain-A func) domain-b 1.0 (domain-b func))
			   (domain-b func))))
      (when (image-A func)
	(setq image-A (if image-A
			  (m* image-A (image-A func))
			  (image-A func)))
	(setq image-b (if image-b
			  (gemm 1.0 image-A (image-b func) 1.0 image-b)
			  (image-b func)))))
    (make-instance
     '<linearly-transformed-function>
     :domain-dimension (if domain-A (ncols domain-A) (domain-dimension func))
     :image-dimension (if image-A (nrows image-A) (image-dimension func))
     :original (typecase func
		 (<linearly-transformed-function> (original func))
		 (t func))
     :domain-A domain-A :domain-b domain-b
     :image-A image-A :image-b image-b)))

(defun intermediate-coordinates (func pos)
  (if (domain-A func)
      (gemm 1.0 (domain-A func) pos 1.0 (domain-b func))
      pos))

(defmethod evaluate ((func <linearly-transformed-function>) pos)
  (let ((result1 (evaluate (original func)
			   (intermediate-coordinates func pos))))
    (if (image-A func)
	(gemm 1.0 (image-A func) result1 1.0 (image-b func))
	result1)))

(defmethod differentiable-p ((f <linearly-transformed-function>) &optional (k 1))
  (or (< k 1)
      (differentiable-p (original f) k)))

(defmethod evaluate-gradient ((func <linearly-transformed-function>) pos)
  (let ((result1 (evaluate-gradient
		  (original func)
		  (intermediate-coordinates func pos))))
    (let ((result2 (aif (image-A func)
			(m* it result1)
			result1)))
      (aif (domain-A func)
	   (m* result2 it)
	   result2))))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; homotopy
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun homotopy (func1 func2)
  "Returns a function which uses its first coordinate as a homotopy
parameter."
  (let ((domain-dim (domain-dimension func1))
	(image-dim (image-dimension func1)))
    (assert (and (= domain-dim (domain-dimension func2))
		 (= image-dim (image-dimension func2))))
    (make-instance
     '<special-function>
     :domain-dimension (1+ domain-dim) :image-dimension image-dim
     :evaluator
     #'(lambda (x)
	 (let ((s (aref x 0))
	       (xrest (vector-slice x 1 domain-dim)))
	   (axpy! (- 1.0 s) (evaluate func1 xrest)
		  (scal s (evaluate func2 xrest)))))
     :gradient
     (and (differentiable-p func1) (differentiable-p func2)
	  #'(lambda (x)
	      (let ((s (aref x 0))
		    (xrest (vector-slice x 1 domain-dim)))
		(join (ensure-matlisp
		       (m- (evaluate func2 xrest) (evaluate func1 xrest)))
		      (axpy! (- 1.0 s) (evaluate-gradient func1 xrest)
			     (scal s (evaluate-gradient func2 xrest))))))))))
  

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; function composition
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun compose (&rest fns)
  (if (null fns)
      #'identity
      (let ((fn1 (car (last fns)))
            (fns (butlast fns)))
        #'(lambda (&rest args)
            (reduce #'funcall fns 
                    :from-end t
                    :initial-value (apply fn1 args))))))

(defmethod compose-2 (f (g function))
  #'(lambda (&rest args)
      (evaluate f (apply g args))))

(defmethod compose-2 ((func1 <function>) (func2 <function>))
  (make-instance
   '<special-function>
   :domain-dimension (domain-dimension func2)
   :image-dimension (image-dimension func1)
   :evaluator
   #'(lambda (x) (evaluate func1 (evaluate func2 x)))
   :gradient
   (and (differentiable-p func1) (differentiable-p func2)
	#'(lambda (x)
	    (m* (evaluate-gradient func1 (evaluate func2 x))
		(evaluate-gradient func2 x))))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Special functions
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun ellipse-matrix (radius excentricity phi)
  "Returns a matrix A suitable for describing the ellipse as (Ax,x)=1."
  (let ((ev1 (make-real-matrix `(,(cos phi) ,(- (sin phi)))))
	(ev2 (make-real-matrix `(,(sin phi) ,(cos phi)))))
    (gemm-nt! #I"((1+excentricity)/radius)^^2" ev1 ev1
	      #I"((1-excentricity)/radius)^^2" (m* ev2 (transpose ev2)))))

(defun project-to-ellipsoid (midpoint A)
  "Returns a function which projects to the ellipsoid given by
Q(x-midpoint)=1 where Q is the quadratic form associated with the matrix
A."
  (let ((dim (length midpoint)))
    (make-instance
     '<special-function> :domain-dimension dim :image-dimension dim
     :evaluator
     #'(lambda (x)
	 (let ((z (m- x midpoint)))
	   (axpy (/ (sqrt (dot z (m* A z)))) z midpoint)))
     :gradient
     #'(lambda (x)
	 (let* ((z (ensure-matlisp (m- x midpoint)))
		(Az (m* A z))
		(Qz (dot z Az))
		(norm-z (sqrt Qz)))
	   #-(or)(gemm-nt! (/ -1 norm-z Qz) z Az
		     (/ norm-z) (eye dim))
	   #+(or)(m- (scal (/ norm-z) (eye dim))
	       (scal (/ 1 norm-z Qz) (m* z (transpose Az)))))))))

(defun project-to-sphere (midpoint radius)
  "Returns a function which projects to the sphere with given midpoint and
radius."
  (project-to-ellipsoid
   midpoint
   (scal (/ (* radius radius)) (eye (length midpoint)))))

(defun xn-distortion-function (f grad-f dim)
  "Returns a function which distorts the xn-coordinate by a factor f(x').
Also grad-f has to be provided."
  (let ((dim-1 (- dim 1)))
    (make-instance
     '<special-function>
     :domain-dimension dim :image-dimension dim
     :evaluator
     #'(lambda (x)
	 (let ((z (copy x)))
	   (setf (aref z dim-1)
		 (* (aref z dim-1)
		    (funcall f (vector-slice x 0 dim-1))))
	   z))
     :gradient
     (and grad-f
	  #'(lambda (x)
	      (let* ((xn (aref x dim-1))
		     (xrest (vector-slice x 0 dim-1))
		     (value (funcall f xrest))
		     (gradient (funcall grad-f xrest))
		     (result (eye dim)))
		(dotimes (i dim-1)
		  (setf (mref result dim-1 i)
			(* (vref gradient i) xn)))
		(setf (mref result dim-1 dim-1) value)
		result))))))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Function handling
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun numerical-gradient (func &key (shift 1.0d-8))
  "Computes the numerical gradient of func at pos."
  #'(lambda (pos)
      (transpose
       (make-real-matrix
	(loop with dim = (length pos)
	      for k below dim collect
	      (coerce
	       (scal (/ (* 2 shift))
		     (m- (evaluate func (axpy shift (unit-vector dim k) pos))
			 (evaluate func (axpy (- shift) (unit-vector dim k) pos))))
	       'list))))))

(defun interval-method (func a b accuracy)
  "Finds zeros of functions in 1d by the interval method."
  (let ((value_a (evaluate func a))
	(value_b (evaluate func b)))
    (labels ((interval-method (a b)
	       (let* ((c (/ (+ a b) 2))
		      (value_c (evaluate func c)))
		 (cond
		   ((< (abs (- b a)) accuracy) c)
		   ((plusp value_c) (if (= b c) c (interval-method a c)))
		   ((minusp value_c) (if (= a c) c (interval-method c b)))
		   (t c)))))
      (cond
	((zerop value_a) a)
	((zerop value_b) b)
	((> (* value_a value_b) 0) '())
	((plusp value_a) (interval-method b a))
	(t (interval-method a b))))))


;;;; Testing

(defun test-function ()
  (let ((project (project-to-sphere #d(0.5 0.5) 0.5)))
    (evaluate project #d(2.0 0.5))
    (evaluate-gradient project #d(2.0 0.5)))
  (evaluate (project-to-ellipsoid #d(0.0 0.0) (eye 2)) #d(3.0 4.0))
  (evaluate-gradient (project-to-ellipsoid #d(-1.0 0.0) (eye 2)) #d(3.0 4.0))
  (evaluate-gradient (project-to-sphere #d(-1.0 0.0) 1.0) #d(3.0 4.0))
  (let ((distortion
	 (xn-distortion-function #'(lambda (x) #I"1+0.5*sin(x[0])")
				 #'(lambda (x) (vector #I"0.5*cos(x[0])"))
				 2))
	(pos #d(1.7 1.0)))
    (evaluate distortion pos)
    (let ((grad (evaluate-gradient distortion pos))
	  (num-grad (evaluate (numerical-gradient distortion) pos)))
      (assert (< (norm (m- grad num-grad)) 1.0e-4))))
  (interval-method #'(lambda (x) (- (* x x) 2.0)) 0.0 2.0 1e-16)
  (let ((f (make-instance '<linear-function> :A (ones 2) :b #d(1.0 0.0))))
    (evaluate f #d(1.0 1.0)))
  )

;;; (test-function)
(fl.tests:adjoin-test 'test-function)
     