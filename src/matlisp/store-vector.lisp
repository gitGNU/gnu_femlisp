;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; ccs.lisp - Compressed column storage scheme
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

(defclass store-vector (<vector>)
  ((store :reader store :initarg :store :documentation
	  "The vector entries."))
  (:documentation "This mixin yields vector behaviour for a class
containing a store.  The store is a unifom array with elements of a certain
type which can be determined by the funtion @function{element-type}.  It
often is but does not have to be equal to the type of scalars for this
vector which can be obtained by calling the function
@function{scalar-type}."))

(defun store-vector-class-name (type)
  (intern (format nil "~A" (list 'store-vector type))
	  :fl.matlisp))

(eval-when (:compile-toplevel :load-toplevel :execute)
  (defun store-vector (type)
    (assert (subtypep type 'number))
    (let ((class-name (store-vector-class-name type)))
      (or (find-class class-name nil)
	  (prog1
	      (eval `(defclass ,class-name (store-vector) ()))
	    (eval `(defmethod element-type ((vector ,class-name))
		    ',type))
	    (eval `(defmethod scalar-type ((vector ,class-name))
		    ',type)))))))

(defmethod initialize-instance :after ((vec store-vector) &key &allow-other-keys)
  "Coerce the store to the correct type."
  (coerce (store vec) `(simple-array ,(element-type vec) (*))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; BLAS macros for store-vector
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun symbol-store (symbol)
  "Returns a symbol of the form <symbol>-STORE."
  (symconc symbol "-STORE"))

(define-blas-macro 'store-vector
  '(element-copy! (x y) `(setf ,y ,x)))

(define-blas-macro 'store-vector
  '(element-scal! (alpha x) `(setf ,x (* ,alpha ,x))))

(define-blas-macro 'store-vector
  '(element-m+! (x y) `(incf ,y ,x)))

(define-blas-macro 'store-vector
  '(element-m+ (x y) `(+ ,y ,x)))

(define-blas-macro 'store-vector
  '(element-m-! (x y) `(decf ,y ,x)))

(define-blas-macro 'store-vector
  '(element-m.*! (x y) `(setf ,y (* ,y ,x))))

(define-blas-macro 'store-vector
  '(element-m* (x y) `(* ,y ,x)))

(define-blas-macro 'store-vector
  '(element-equal (x y) `(= ,y ,x)))

(define-blas-macro 'store-vector
  '(element-gemm! (alpha x y beta z)
    `(setf ,z (+ (* ,alpha ,x ,y) (* ,beta ,z)))))

(define-blas-macro 'store-vector
    '(with-blas-data (vars &rest body)
      "Gets the stores of all vector variables."
      `(let ,(loop for var in vars collecting
		   `(,(symbol-store var) (slot-value ,var 'store)))
	(declare (ignorable ,@(loop for var in vars collecting (symbol-store var))))
	(declare (type (simple-array element-type (*))
		  ,@(mapcar #'symbol-store vars)))
	,@body)))

(define-blas-macro 'store-vector
    '(for-each-entry (loop-vars &rest body)
      (let ((i (gensym "I")))
	`(symbol-macrolet
	  ,(loop for (xc x) in loop-vars
		 collect `(,xc (aref ,(symbol-store x) ,i)))
	  (dotimes (,i (length ,(symbol-store (cadar loop-vars))))
	    ,@body)))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; BLAS operations
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmacro assert-store-vector-compatibility (x y)
  `(unless (= (length ,(symbol-store x)) (length ,(symbol-store y)))
    (error "Store-vectors are incompatible.")))

;;; Copying

(new-define-blas-template copy! ((x store-vector) (y store-vector))
  (declare (optimize (speed 3) (debug 0) (safety 0)))
  (for-each-entry ((xc x) (yc y)) (setf yc xc))
  y)

;;; BLAS routines involving one store-vector

(new-define-blas-template fill! ((x store-vector) (s number))
  (declare (optimize (speed 3) (debug 0) (safety 0)))
  (for-each-entry ((xc x)) (element-copy! s xc))
  x)

(new-define-blas-template fill-random! ((x store-vector) (s number))
  (declare (optimize (speed 3) (debug 0) (safety 0)))
  (for-each-entry ((xc x)) (element-copy! (random s) xc))
  x)

(new-define-blas-template scal! ((alpha number) (x store-vector))
  (declare (optimize (speed 3) (debug 0) (safety 0)))
  (for-each-entry ((xc x)) (element-scal! alpha xc))
  x)

;;; BLAS routines involving two store-vectors

(new-define-blas-template axpy! ((alpha number) (x store-vector) (y store-vector))
  (declare (optimize speed))
  (assert-store-vector-compatibility x y)
  (for-each-entry ((xc x) (yc y)) (element-m+! (element-m* alpha xc) yc))
  y)

(new-define-blas-template dot ((x store-vector) (y store-vector))
  (declare (optimize speed))
  (let ((sum (coerce 0 'element-type)))
    (declare (type element-type sum))
    (for-each-entry ((xc x) (yc y))
       (element-m+! (element-m* xc yc) sum))
    sum))

(new-define-blas-template mequalp ((x store-vector) (y store-vector))
  "Exact equality test for store-vector."
  (assert-store-vector-compatibility x y)
  (for-each-entry ((xc x) (yc y))
     (unless (element-equal xc yc)
       (return-from mequalp nil)))
  t)

(new-define-blas-template dot-abs ((x store-vector) (y store-vector))
  (declare (optimize speed))
  (let ((sum (coerce 0 'element-type)))
    (declare (type element-type sum))
    (for-each-entry ((xc x) (yc y))
       (element-m+! (abs (element-m* xc yc)) sum))
    sum))

(new-define-blas-template m+! ((x store-vector) (y store-vector))
  (declare (optimize speed))
  (assert-store-vector-compatibility x y)
  (for-each-entry ((xc x) (yc y)) (element-m+! xc yc))
  y)

(new-define-blas-template m.*! ((x store-vector) (y store-vector))
  (declare (optimize speed))
  (assert-store-vector-compatibility x y)
  (for-each-entry ((xc x) (yc y)) (element-m.*! xc yc))
  y)

(defun store-vector-generator (type)
  #'(lambda (n)
      (make-instance (store-vector type)
		     :store (make-array n :element-type type))))

(defun test-store-vector-blas ()
  (dbg-on :blas)
  (let ((x (make-instance (store-vector 'double-float) :store #d(0.0)))
	(y (make-instance (store-vector 'double-float) :store #d(1.0))))
    (describe x)
    (m+! y x)
    (scal! 0.5 x)
    (copy! x y))
  (test-blas 'm+! 10 :generator (store-vector-generator 'single-float))
  (test-blas 'm.*! 10 :generator (store-vector-generator 'double-float))
  (test-blas 'dot 1 :generator (store-vector-generator 'single-float))
  (dbg-off :blas)
  )


