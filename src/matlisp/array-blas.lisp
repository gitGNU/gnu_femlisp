;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; arrays.lisp - Matlisp for arrays
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

;;; Ensure some BLAS operations also for Lisp vectors and arrays

(in-package :fl.matlisp)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; double-vec
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(deftype double-vec () '(simple-array double-float (*)))

(definline list->double-vec (lst)
  (map 'double-vec #'(lambda (x) (coerce x 'double-float)) lst))

(definline double-vec (&rest comps)
  (list->double-vec comps))

(definline make-double-vec (dim &optional (init 0.0))
  "make-double-vec: double-vec constructor"
  (make-array dim :element-type 'double-float :initial-element (float init 0.0)))

(eval-when (:compile-toplevel :load-toplevel :execute)
  (set-dispatch-macro-character
   #\# #\d  ; dispatch on #d for double-vec
   #'(lambda (stream char n)
       (declare (ignore char n))
       (let ((list (read stream nil (values) t)))
	 `(list->double-vec ',list)))))

(definline unit-vector (dim i)
  (let ((vec (make-double-vec dim)))
    (setf (aref vec i) 1.0)
    vec))

(defgeneric multiplicity (vec)
  (:documentation "We allow multiple vectors, for solving linear problems
in parallel."))

(defmethod multiplicity (vec)
  "If vec is a matrix, the multiplicity is the number of columns by
default."
  (ncols vec))
  
;;; Sequences are considered generally as column vectors

(defmethod vlength ((seq sequence))
  "For sequences, the number of rows is equal to the length."
  (length seq))
(defmethod nrows ((seq sequence))
  "For sequences, the number of rows is equal to the length."
  (length seq))
(defmethod ncols ((seq sequence))
  "For sequences, the number of columns is 1."
  1)

(defmethod element-type ((x vector))
  (array-element-type x))
(defmethod scalar-type ((x vector))
  (array-element-type x))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; conversion of vectors to standard-matrix
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun ensure-matlisp (vec &optional (type :column))
  (etypecase vec
    (vector
     (let ((element-type (array-element-type vec)))
       (assert (subtypep element-type 'number)) ; preliminary
       (make-instance
	(standard-matrix element-type)
	:store vec
	:nrows (if (eq type :column) (length vec) 1)
	:ncols (if (eq type :row) (length vec) 1))))
    (standard-matrix vec)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; vector operations
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmethod vref ((vec vector) index)
  (aref vec index))
(defmethod (setf vref) (val (vec vector) index)
  (setf (aref vec index) val))

(defmethod for-each-key (func (x vector))
  (dotimes (i (length x)) (funcall func i)))

(defmethod for-each-entry (func (seq sequence))
  (map nil func seq))

(defmethod for-each-entry (func (array array))
  "Call @var{func} on each entry of @var{array}."
  (dotimes (i (array-total-size array))
    (funcall func (row-major-aref array i))))

(defmethod for-each-key-and-entry (func (x vector))
  (dotimes (i (length x))
    (funcall func i (aref x i))))

;;; These do not correspond exactly to their name 
(defmethod for-each-entry-of-vec1 (func (x sequence) (y sequence))
  (map nil func x y))
(defmethod for-each-entry-of-vec2 (func (x sequence) (y sequence))
  (map nil func x y))

(defmethod for-each-entry-of-vec1 (func (x array) (y array))
  (array-for-each func x y))
(defmethod for-each-entry-of-vec2 (func (x array) (y array))
  (array-for-each func x y))


;;; matrix-vector multiplication
(defmethod x<-Ay (x (A standard-matrix) y)
  (gemm! 1.0 A (ensure-matlisp y) 0.0 (ensure-matlisp x))
  x)
(defmethod x+=Ay (x (A standard-matrix) y)
  (gemm! 1.0 A (ensure-matlisp y) 1.0 (ensure-matlisp x))
  x)
(defmethod x-=Ay (x (A standard-matrix) y)
  (gemm! -1.0 A (ensure-matlisp y) 1.0 (ensure-matlisp x))
  x)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Matlisp operations for arrays
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; matrix-vector multiplication
(defmethod x<-Ay (x (A standard-matrix) y)
  (gemm! 1.0 A (ensure-matlisp y) 0.0 (ensure-matlisp x))
  x)
(defmethod x+=Ay (x (A standard-matrix) y)
  (gemm! 1.0 A (ensure-matlisp y) 1.0 (ensure-matlisp x))
  x)
(defmethod x-=Ay (x (A standard-matrix) y)
  (gemm! -1.0 A (ensure-matlisp y) 1.0 (ensure-matlisp x))
  x)

(defmethod total-entries ((mat standard-matrix))
   (* (nrows mat) (ncols mat)))

;;; copying of arrays

(defmethod copy ((vec vector)) (map (type-of vec) #'copy vec))
(defmethod copy ((vec list)) (mapcar #'copy vec))

(defmethod copy ((array array) &aux (result (make-analog array)))
  (dotimes (i (array-total-size array) result)
    (setf (row-major-aref result i)
	  (copy (row-major-aref array i)))))

;;; BLAS level 1 for Lisp vectors

(defmacro define-vector-blas-method (symbol args &body body)
  (let ((vec-args (loop for (name arg) in args when (eq arg 'vector) collect name))
	(number-args (loop for (name arg) in args when (eq arg 'number) collect name)))
    (with-gensyms (types functions pos)
      `(let ((,types (vector))
	     (,functions (vector)))
	(defmethod ,symbol ,args
	  (let ((element-type (scalar-type ,(first vec-args))))
	    ;; ensure that function is on position 0
	    (unless (and (plusp (length ,types)) (equalp element-type (aref ,types 0)))
	      (let ((,pos (position element-type ,types)))
		(cond
		  (,pos (rotatef (aref ,types ,pos) (aref ,types 0))
			(rotatef (aref ,functions ,pos) (aref ,functions 0)))
		  (t (setq ,types (concatenate 'vector (list element-type) ,types))
		     (setq ,functions
			   (concatenate
			    'vector
			    (list
			     (let ((*compile-print* nil))
			       (compile
				nil
				(let ((method-source
				       (subst element-type 'element-type
					      '(lambda ,(mapcar #'car args)
						(declare (type (simple-array element-type (*)) ,@vec-args))
						(declare (type element-type ,@number-args))
						(declare (optimize (speed 3) (safety 0)))
						,@body))))
				  (dbg :blas "Generated method: ~%~S~%" method-source)
				  method-source))))
			    ,functions)))))))
	  (funcall (aref ,functions 0) ,@(mapcar #'car args)))))))

(define-vector-blas-method m+! ((x vector) (y vector))
  (dotimes (i (length x) y)
    (incf (aref x i) (aref y i))))

(define-vector-blas-method copy! ((x vector) (y vector))
  (dotimes (i (length x) y)
    (setf (aref y i) (aref x i))))

(define-vector-blas-method m+! ((x vector) (y vector))
  (dotimes (i (length x) y)
    (incf (aref y i) (aref x i))))

#+(or) ; is not generic
(define-vector-blas-method m-! ((x vector) (y vector))
  (dotimes (i (length x) y)
    (decf (aref y i) (aref x i))))

(define-vector-blas-method fill! ((x vector) (s number))
  (dotimes (i (length x) x)
    (setf (aref x i) s)))

(define-vector-blas-method fill-random! ((x vector) (s number))
  (dotimes (k (length x) x)
    (setf (aref x k) (fill-random! x s))))

(define-vector-blas-method scal! ((val number) (vec vector))
  (dotimes (i (length vec) vec)
    (setf (aref vec i) (* (aref vec i) val))))

(defmethod scal! (val (lst list))
  "Call @func{scal!} on each entry of @var{lst}.  Note that the result has
to be freshly consed, because @func{scal!} is a function, not a macro."
  (mapcar #'(lambda (entry) (scal! val entry)) lst))

(defmethod scal! (val (array array))
  "Call @func{scal!} on each entry of @var{array}."
  (call-next-method))

(define-vector-blas-method axpy! ((alpha number) (x vector) (y vector))
  (dotimes (i (length x) y)
    (incf (aref y i) (* alpha (aref x i)))))

(defmethod dot ((x list) (y list))
  (loop for xc in x and yc in y summing (* xc yc)))

(defmacro m-incf (result increment)
  "Adds increment to result which should be a symbol.  If its value is nil
then result is set to increment."
  (with-gensyms (inc)
    `(let ((,inc ,increment))
      (if ,result
	  (if (typep ,inc 'number)
	      (incf ,result ,inc)
	      (m+! ,inc ,result))
	  (setf ,result ,inc)))))

#+(or)
(defmethod gemm! (alpha (A array) (B array) beta (C array) &optional (job :nn))
  "Gemm! using 2d-arrays.  Entries will often be matrices."
  (declare (type (array t (* *)) A B C))
  (let* ((i-dim (array-dimension C 0))
	 (j-dim (array-dimension C 1))
	 (A-transposed-p (member job '(:tn :tt)))
	 (B-transposed-p (member job '(:nt :nt)))
	 (contract-dim (array-dimension A (if A-transposed-p 0 1))))
    (assert (and (= i-dim (array-dimension A (if A-transposed-p 1 0)))
		 (= j-dim (array-dimension B (if B-transposed-p 0 1)))
		 (= contract-dim (array-dimension B (if B-transposed-p 1 0)))))
    (dotimes (i (array-dimension C 0))
      (dotimes (k (array-dimension C 1))
	(dotimes (j contract-dim)
	  (let ((aij (if A-transposed-p (aref A j i) (aref A i j)))
		(bjk (if B-transposed-p (aref B j k) (aref B k j)))
		(cik (aref C i k)))
	    (if (typep cik 'number)
		(setf (aref C i k) (+ (* beta cik) (* alpha (dot aij bjk))))
		(gemm! alpha aij bjk beta cik))))))))

;;; Multiplication of Matlisp matrices and Lisp vectors

(defmethod m*-product-instance ((A standard-matrix) (x vector))
  "Returns a CL vector as result of matrix * CL vector."
  (make-array (nrows A) :element-type (array-element-type x)))

(defmethod gemm-nn! (alpha (A standard-matrix) (x vector) beta (y vector))
  (gemm-nn! alpha A (ensure-matlisp x) beta (ensure-matlisp y))
  y)

(defmethod getrs! ((lr standard-matrix) (b vector) &optional ipiv)
  (getrs! lr (ensure-matlisp b) ipiv)
  b)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Testing
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun test-arrays ()
  (m+! #d(0.0) #d(0.0))
  (m+! #(0.0) #(0.0))
  (let* ((dim 2)
	 (x (make-double-vec dim))
	 (y (make-double-vec dim 1.0)))
    (scalar-type x)
    (x<-0 x)
    (m+! y x)
    (axpy! 2.0 y x)
    (m-! y x)
    (copy #(1.0))
    (format t "X=~A Y=~A" x y)
    (scal 0.5
	  (m- (axpy 0.01 x y)
	      (axpy -0.01 x y)))
    (m* (ensure-matlisp x)
	(transpose (ensure-matlisp y)))
    (gemm 1.0 (eye 2) (ensure-matlisp x) 1.0 (ensure-matlisp y))
    (norm #(1.0 2.0) 1)
  ))

;;; (test-arrays)
(fl.tests:adjoin-test 'test-arrays)


