;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; vector.lisp - Linear algebra interface for vectors
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

(defclass <vector> ()
  ()
  (:documentation "General vector class."))

(defgeneric vlength (vec)
  (:documentation "Length of vector."))
(defgeneric multiplicity (vec)
  (:documentation "We allow multiple vectors, for solving linear problems
in parallel."))


(defgeneric element-type (vector)
  (:documentation "Type of the elements of the vector/matrix."))

(defgeneric scalar-type (vector)
  (:documentation "Type of the scalars for the vector class."))

(defgeneric total-entries (vector)
  (:documentation "Total number of entries for block vectors."))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Linear algebra interface
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Primary routines
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; Reader and writers
(defgeneric vref (x i)
  (:documentation "Reader for @math{x_i}."))
(defgeneric (setf vref) (value x i)
  (:documentation "Writer for @arg{x_i}."))

;;; Copy
(defgeneric copy (x)
  (:documentation "Returns a deep copy of X.")
  (:method (x)
    "This default method for @function{copy} does a destructive copy to an
empty analog of @arg{x}."
    (copy! x (make-analog x))))


;;; Filling
(defgeneric fill! (x s)
  (:documentation "Fills X with element s."))
(defgeneric fill-random! (x s)
  (:documentation "Fills X with random values (obtained by (random s))."))

;;; Joining
(defgeneric join-horizontal! (x y z)
  (:documentation "Fills Z with the horizontal join of X and Y."))
(defgeneric join-vertical! (x y z)
  (:documentation "Fills Z with the vertical join of X and Y."))

;;; BLAS Level 1
(defgeneric scal! (alpha x)
  (:documentation "X <- \alpha X"))
(defgeneric copy! (x y)
  (:documentation "Y <- X"))
(defgeneric axpy! (alpha x y)
  (:documentation "Y <- alpha*X + Y"))
(defgeneric m+! (x y)
  (:documentation "Y <- X + Y"))
(defgeneric m.*! (x y)
  (:documentation "Y <- X .* Y"))

(defgeneric mzerop (x &optional threshold)
  (:documentation "Returns T if each entry of @arg{x} is smaller or equal
than @arg{threshold}."))
(defgeneric mequalp (x y)
  (:documentation "Returns T if X and Y have equal entries, otherwise NIL."))

(defgeneric dot (x y)
  (:documentation "Returns the dot product of X and Y."))
(defgeneric dot-abs (x y)
  (:documentation "Returns the dot product between |X| and |Y|."))

(defgeneric norm (x &optional p)
  (:documentation "Returns the @arg{p}-norm of @arg{x}."))
(defgeneric l2-norm (x)
  (:documentation "Returns the 2-norm of @arg{x}."))
(defgeneric lp-norm (x p)
  (:documentation "Returns the @arg{p}-norm of @arg{x}."))
(defgeneric linf-norm (x)
  (:documentation "Returns the maximum norm of @arg{x}."))

;;; Iteration

;;; The following is the old interface and will be replaced rather soon by
;;; a macro based approach which allows suitable inlining.

(defgeneric for-each-key (func vec))
(defgeneric for-each-entry (func vec))
(defgeneric for-each-key-and-entry (func vec))
(defgeneric for-each-entry-of-vec1 (func vec1 vec2)
  (:documentation "Calls @arg{func} on all entries of @arg{vec1} with this
entry and the corresponding entry of @arg{vec2}.  Deprecated, use
@function{for-each-key-and-entry} instead.")
  (:method ((fn function) vec1 vec2)
  (for-each-key-and-entry  #'(lambda (key entry)
			       (funcall fn entry (vref vec2 key)))
			   vec1)))

(defgeneric for-each-entry-of-vec2 (func vec1 vec2)
  (:documentation "Calls @arg{func} on all entries of @arg{vec2} with the
corresponding entry of @arg{vec1} and that entry.  Deprecated.")
  (:method ((fn function) vec1 vec2)
  (for-each-key-and-entry  #'(lambda (key entry)
			       (funcall fn (vref vec1 key) entry))
			   vec2)))

(defmacro dovec ((loop-vars vec) &body body)
  "Loops on indices and entries of a vector.  Examples:
@lisp
  (dovec ((key) vec) ...)
  (dovec (entry vec) ...)
  (dovec ((key entry) vec) ...)
@end lisp"
  (let* ((loop-vars (if (consp loop-vars) loop-vars (list nil loop-vars)))
	 (vec-for-each (if (null (car loop-vars))
			   (if (null (cdr loop-vars))
			       (error "no loop variable")
			       'for-each-entry)
			   (if (or (null (cdr loop-vars)) (null (cadr loop-vars)))
			       'for-each-key
			       'for-each-key-and-entry))))
    `(,vec-for-each #'(lambda ,(remove nil loop-vars) ,@body) ,vec)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Derived functionality
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; This is either defined by inlining functions, if we expect those rewritings
;;; to be valid for every vector or matrix class.  Alternatively, it may be
;;; given by a generic function.

(definline x<-0 (x)
  "X <- 0 X.  Uses SCAL!."
  (scal! (coerce 0 (scalar-type x)) x))

(definline scal (alpha x)
  "Returns alpha * X.  Uses SCAL! and COPY."
  (if (numberp x)
      (* alpha x)
      (scal! alpha (copy x))))

(definline normalize! (x &optional (p 2))
  "Scales @arg{x} destructively to have @arg{p}-norm equal to 1."
  (scal! (/ (norm x p)) x))

(definline normalize (x &optional (p 2))
  "Scales @arg{x} to have @arg{p}-norm equal to 1."
  (scal (/ (norm x p)) x))

(defun axpy (alpha x y)
  "Returns alpha X + Y.  Uses AXPY! and COPY."
  (if (numberp x)
      (+ (* alpha x) y)
      (axpy! alpha x (copy y))))

(defgeneric m+ (x y)
  (:documentation "Returns @math{X} + @math{Y}."))
(defmethod m+ (x y)
  "Default method uses M+! and COPY."
  (m+! x (copy y)))

(definline m-! (x y)
  "Y - X -> Y.  Uses AXPY!."
  (axpy! (coerce -1 (scalar-type x)) x y))

(defun m- (x y)
  "Returns X-Y.  Uses AXPY."
  (if (or (numberp x) (numberp y))
      (- x y)
      (axpy (coerce -1 (scalar-type x)) y x)))

(definline m.* (x y)
  "Returns X .* Y.  Uses M.*! and COPY."
  (m.*! x (copy y)))

(defmethod norm (x &optional (p 2))
  (case p
    (2 (l2-norm x))
    (:inf (linf-norm x))
    (t (lp-norm x p))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; some core functionality is recursively defined for block vectors
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmethod fill! (x s)
  "Recursive definition for FILL! usable for sparse block vectors."
  (for-each-entry (rcurry #'fill! s) x)
  x)
(defmethod fill-random! (x s)
  "Recursive definition for FILL-RANDOM! usable for sparse block vectors."
  (for-each-entry (rcurry #'fill-random! s) x)
  x)

(defmethod total-entries (obj)
  (let ((entries 0))
    (for-each-entry
     #'(lambda (entry)
	 (incf entries (total-entries entry)))
     obj)
    entries))

(defmethod mzerop (mat &optional (threshold 0.0))
  (for-each-entry
   #'(lambda (entry)
       (unless (mzerop entry threshold)
	 (return-from mzerop nil)))
   mat)
  t)

(defmethod copy! (x y)
  "Recursive definition for COPY! usable for sparse block vectors."
  (for-each-key #'(lambda (i)
		    (setf (vref y i)
			  (copy! (vref x i) (vref y i))))
		x)
  y)

(defmethod m+! (x y)
  "Recursive definition for M+! usable for sparse block vectors."
  (for-each-key #'(lambda (i)
		    (setf (vref y i)
			  (m+! (vref x i) (vref y i))))
		x)
  y)

(defmethod scal! (s x)
  (for-each-entry #'(lambda (x) (scal! (coerce s (scalar-type x)) x))
		  x)
  x)

(defmethod axpy! (alpha x y)
  "Recursive definition for AXPY! usable for sparse block vectors."
  (for-each-entry-of-vec1
   #'(lambda (x y)
       (axpy! (coerce alpha (scalar-type x)) x y))
   x y)
  y)

(defmethod l2-norm (vec)
  "Recursive definition for the l2-norm."
  (let ((sum 0))
    (for-each-entry #'(lambda (x) (incf sum (expt (l2-norm x) 2)))
		    vec)
    (sqrt sum)))

(defmethod lp-norm (vec (p number))
  "Recursive definition for the lp-norm."
  (let ((sum 0))
    (for-each-entry #'(lambda (x) (incf sum (expt (lp-norm x p) p)))
		    vec)
    (expt sum (/ 1 p))))

(defmethod linf-norm (vec)
  "Recursive definition for the linf-norm."
  (let ((max 0))
    (for-each-entry #'(lambda (x)
			(let ((norm-x (linf-norm x)))
			  (when (> norm-x max) (setq max norm-x))))
		    vec)
    max))

(defmethod dot (x y)
  (let ((sum 0))
    (for-each-entry-of-vec1 #'(lambda (x y) (incf sum (dot x y)))
			    x y)
    sum))

(defmethod dot-abs (x y)
  (let ((sum 0))
    (for-each-entry-of-vec1 #'(lambda (x y) (incf sum (dot-abs x y)))
			    x y)
    sum))

(defmethod mequalp (x y)
  (when (= (total-entries x) (total-entries y))
    (for-each-entry-of-vec1
     #'(lambda (x y)
	 (unless (mequalp x y)
	   (return-from mequalp nil)))
     x y)
    t))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Testing
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun test-vector ()
  )

;;; (test-vector)
(fl.tests:adjoin-test 'test-vector)
