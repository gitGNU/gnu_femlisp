;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; store-vector.lisp - Compressed column storage scheme
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
  ((store :reader store :initarg :store
	  :documentation "The vector entries."))
  (:documentation "This mixin yields vector behaviour for a class
containing a store.  The store is a unifom array with elements of a certain
type which can be determined by the funtion @function{element-type}.  It
often is but does not have to be equal to the type of scalars for this
vector which can be obtained by calling the function
@function{scalar-type}."))

(defclass static-store-vector (store-vector) ()
  (:documentation "Subclass of @class{dynamic-store-vector} for which fast
BLAS operations are possible because the store is a simple and uniform
array."))

(defun number-super-type (types)
  "Returns the union of the number types in @arg{types}."
  (reduce (lambda (t1 t2)
            (if (eql t1 t2)
                t1
                (loop for (set type) in
                     '(((single-float (complex single-float)) (complex single-float))
                       ((double-float (complex double-float)) (complex double-float)))
                     unless (set-exclusive-or set (list t1 t2) :test 'equalp)
                     do (return type)
                     finally (return 'number))))
          types))

(defun uniform-number-type (vec)
  "Tries to find a uniform type for the numbers contained in @arg{vec}."
  (let ((element-type (array-element-type vec)))
    (if (subtypep element-type 'number)
        element-type
        (if (zerop (length vec))
            'number
            (number-super-type (map 'list #'type-of vec))))))

(with-memoization (:type :global :size 4 :id 'store-vector :test 'equal :debug t)
  (defun store-vector (type &key dynamic)
    (memoizing-let ((type type)
                    (dynamic dynamic))
      (let ((class-name
	     (intern (format nil "~A" (list 'store-vector type :dynamic dynamic))
		     "FL.MATLISP")))
	(or (find-class class-name nil)
	    (prog1
		(fl.port:compile-and-eval
		 `(defclass ,class-name
		   (,(if dynamic 'store-vector 'static-store-vector))
		   ()))
	      (fl.port:compile-and-eval
	       `(defmethod element-type ((vector ,class-name)) ',type))
	      (fl.port:compile-and-eval
	       `(defmethod scalar-type ((vector ,class-name)) ',type))))))))

#+(or)  ; here is a performance problem (memoization, shared-initialize, etc)
(time (let ((class (store-vector 'double-float)))
	(loop repeat 100000 do (make-instance class :store nil))))

(defmethod shared-initialize :after ((vec static-store-vector) slot-names
				     &key &allow-other-keys)
  "Coerce the store to a static type."
  (declare (ignore slot-names))
  (when (slot-boundp vec 'store)
    (setf (slot-value vec 'store)
	  (coerce (store vec)
                  `(simple-array ,(upgraded-array-element-type (element-type vec))
                                 (*))))))

(defmethod vref ((vec store-vector) i)
  "Note: access to entries using this generic function is slow.  Therefore,
specialized BLAS routines should be used whenever possible."
    (aref (slot-value vec 'store) i))

(defmethod (setf vref) (value (vec store-vector) i)
  "Writer for vref."
  (setf (aref (slot-value vec 'store) i) value))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; BLAS macros for static-store-vector
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun symbol-store (symbol)
  "Returns a symbol of the form <symbol>-STORE."
  (symconc symbol "-STORE"))

(define-blas-macro 'static-store-vector
  '(element-copy! (x y) `(setf ,y ,x)))

(define-blas-macro 'static-store-vector
  '(element-scal! (alpha x) `(setf ,x (* ,alpha ,x))))

(define-blas-macro 'static-store-vector
  '(element-m+! (x y) `(incf ,y ,x)))

(define-blas-macro 'static-store-vector
  '(element-m+ (x y) `(+ ,y ,x)))

(define-blas-macro 'static-store-vector
  '(element-m-! (x y) `(decf ,y ,x)))

(define-blas-macro 'static-store-vector
  '(element-m.*! (x y) `(setf ,y (* ,y ,x))))

(define-blas-macro 'static-store-vector
  '(element-m* (x y) `(* ,y ,x)))

(define-blas-macro 'static-store-vector
  '(element-equal (x y) `(= ,y ,x)))

(define-blas-macro 'static-store-vector
  '(element-gemm! (alpha x y beta z)
    `(setf ,z (+ (* ,alpha ,x ,y) (* ,beta ,z)))))

(define-blas-macro 'static-store-vector
    '(with-blas-data (vars &rest body)
      "Gets the stores of all vector variables."
      `(let ,(loop for var in vars collecting
		   `(,(symbol-store var) (slot-value ,var 'store)))
	(declare (type (simple-array element-type (*))
		  ,@(mapcar #'symbol-store vars)))
	,@body)))

(define-blas-macro 'static-store-vector
    '(assert-vector-operation-compatibility (x y)
      `(unless (= (length ,(symbol-store x)) (length ,(symbol-store y)))
	 (error "Store-vectors are incompatible."))))

(define-blas-macro 'static-store-vector
  '(vec-for-each-entry (loop-vars &rest body)
    (let ((i (gensym "I")))
      `(with-blas-data ,loop-vars
	,@(when (= (list-length loop-vars) 2)
		`((assert-vector-operation-compatibility ,@loop-vars)))
	(dotimes (,i (length ,(symbol-store (first loop-vars))))
	  (macrolet ((ref (matrix) (list 'aref (symbol-store matrix) ',i)))
	    ,@body))))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; BLAS operations
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;; Copying

(define-blas-template copy! ((x static-store-vector) (y static-store-vector))
  (vec-for-each-entry (x y) (setf (ref y) (ref x)))
  y)

;;; BLAS routines involving one static-store-vector

(define-blas-template fill! ((x static-store-vector) (s number))
  (vec-for-each-entry (x) (element-copy! s (ref x)))
  x)

(define-blas-template fill-random! ((x static-store-vector) (s number))
  (vec-for-each-entry (x) (element-copy! (random s) (ref x)))
  x)

(define-blas-template scal! ((alpha number) (x static-store-vector))
  (vec-for-each-entry (x) (element-scal! alpha (ref x)))
  x)

;;; BLAS routines involving two static-store-vectors

(define-blas-template axpy! ((alpha number) (x static-store-vector) (y static-store-vector))
  (vec-for-each-entry (x y) (element-m+! (element-m* alpha (ref x)) (ref y)))
  y)

(define-blas-template dot ((x static-store-vector) (y static-store-vector))
  (let ((sum (coerce 0 'element-type)))
    (declare (type element-type sum))
    (vec-for-each-entry (x y)
       (element-m+! (element-m* (ref x) (ref y)) sum))
    sum))

(define-blas-template mequalp ((x static-store-vector) (y static-store-vector))
  "Exact equality test for static-store-vector."
  (vec-for-each-entry (x y)
     (unless (element-equal (ref x) (ref y))
       (return-from mequalp nil)))
  t)

(define-blas-template dot-abs ((x static-store-vector) (y static-store-vector))
  (let ((sum (coerce 0 'element-type)))
    (declare (type element-type sum))
    (vec-for-each-entry (x y)
       (element-m+! (abs (element-m* (ref x) (ref y))) sum))
    sum))

(define-blas-template m+! ((x static-store-vector) (y static-store-vector))
  (vec-for-each-entry (x y) (element-m+! (ref x) (ref y)))
  y)

(define-blas-template m.*! ((x static-store-vector) (y static-store-vector))
  (vec-for-each-entry (x y) (element-m.*! (ref x) (ref y)))
  y)

(defun store-vector-generator (type)
  #'(lambda (n)
      (make-instance (store-vector type)
		     :store (zero-vector n type))))

(defun test-store-vector-blas ()
  (dbg-on :blas)
  (dbg-on :compile)

  (let ((x (make-instance (store-vector 'double-float) :store #d(0.0)))
	(y (make-instance (store-vector 'double-float) :store #d(1.0))))
    (copy! x y)
    (describe x)
    (m+! y x)
    (axpy! -1.0 x y)
    ;; (axpy! -1 x y) ; leads to segfault for Allegro CL - not nice
    (scal! 0.5 x))
  (test-blas 'copy! 1 :generator (store-vector-generator '(complex double-float)))
  (test-blas 'm+! 1000 :generator (store-vector-generator 'single-float))
  (test-blas 'm.*! 10 :generator (store-vector-generator 'double-float))
  (test-blas 'dot 1 :generator (store-vector-generator 'single-float))
  (dbg-off)
  )

;;; (test-store-vector-blas)
(fl.tests:adjoin-test 'test-store-vector-blas)
