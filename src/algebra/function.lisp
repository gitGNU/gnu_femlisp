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

(in-package :algebra)

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

(defmethod evaluate-k-jet ((f <special-function>) (k fixnum) x)
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

(defmethod smoothness ((f <constant-function>))
  :infinity)

(defmethod evaluate-gradient ((f <constant-function>) x)
  (assert (eq (length x) (domain-dimension f)))
  (make-float-matrix (image-dimension f) (domain-dimension f)))

(defmethod evaluate-k-jet ((f <constant-function>) (k fixnum) x)
  (assert (eq (length x) (domain-dimension f)))
  (cond ((= k 0) (value f))
	((= k 1) (make-float-matrix (image-dimension f) (domain-dimension f)))
	(t (make-real-tensor
	    (list->fixnum-vec
	     (cons (image-dimension f)
		   (make-list k :initial-element (domain-dimension f))))))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; transformed functions
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass <linearly-transformed-function> (<function>)
  ((original :reader original :initarg :original)
   (A :reader trans-A :initarg :A)
   (b :reader trans-b :initarg :b))
  (:documentation "<linearly-transformed-function> calls the original
function evaluation on transformed coordinates."))

(defun transform-function (func A b)
  (assert (= (domain-dimension func) (nrows A)))
  (typecase func
    (<linearly-transformed-function>
     (make-instance
      '<linearly-transformed-function>
      :domain-dimension (ncols A)
      :image-dimension (image-dimension func)
      :original (original func)
      :A (m* (trans-A func) A)
      :b (m+ (m* (trans-A func) b) (trans-b func))))
    (t 
     (make-instance
      '<linearly-transformed-function>
      :domain-dimension (ncols A)
      :image-dimension (image-dimension func)
      :original func :A A :b b))))

(defun intermediate-coordinates (func pos)
  (m+ (m* (trans-A func) pos) (trans-b func)))

(defmethod evaluate ((func <linearly-transformed-function>) pos)
  (evaluate (original func) (intermediate-coordinates func pos)))

(defmethod evaluate-gradient ((func <linearly-transformed-function>) pos)
  (m* (evaluate-gradient (original func) (intermediate-coordinates func pos))
      (trans-A func)))

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
	   (vec+ (scal (- 1.0 s) (evaluate func1 xrest))
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

(defun project-to-sphere (midpoint radius)
  "Returns a function which projects to the sphere with given midpoint and
radius."
  (let ((dim (length midpoint)))
    (make-instance
     '<special-function> :domain-dimension dim :image-dimension dim
     :evaluator
     #'(lambda (x)
	 (let ((z (m- x midpoint)))
	   (axpy (/ radius (norm z)) z midpoint)))
     :gradient
     #'(lambda (x)
	 (let* ((z (ensure-matlisp (m- x midpoint)))
		(norm-z (norm z))
		(unit-z (scal (/ norm-z) z)))
	   (scal (/ radius norm-z)
		 (gemm -1.0 unit-z unit-z 1.0 (eye dim) :nt)))))))

(defun xn-distortion-function (f grad-f dim)
  "Returns a function which distorts the xn-coordinate by a factor f(x').
Also grad-f has to be provided."
  (declare (optimize (debug 3)))
  (let ((dim-1 (- dim 1)))
    (make-instance
     '<special-function>
     :domain-dimension dim :image-dimension dim
     :evaluator
     #'(lambda (x)
	 (let ((z (copy x)))
	   (setf (vec-ref z dim-1)
		 (* (vec-ref z dim-1)
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
		  (setf (matrix-ref result dim-1 i)
			(* (vec-ref gradient i) xn)))
		(setf (matrix-ref result dim-1 dim-1) value)
		result))))))



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; finding roots of functions in 1d
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun interval-method (func a b accuracy)
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


;;;; Testing: (test-function)

(defun test-function ()
  (let ((project (project-to-sphere (double-vec 0.5 0.5) 0.5)))
    (evaluate project (double-vec 1.0 0.5))
    (evaluate-gradient project (double-vec 1.0 0.5)))
  (let ((distortion
	 (xn-distortion-function #'(lambda (x) #I"1+0.5*sin(x[0])")
				 #'(lambda (x) (vector #I"0.5*cos(x[0])"))
				 2)))
  (evaluate-gradient distortion #(1.7 1.0)))
  (interval-method #'(lambda (x) (- (* x x) 2.0)) 0.0 2.0 1e-16))

(tests::adjoin-femlisp-test 'test-function)
     
