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

;;; Ensure some Matlisp operations also for Lisp vectors and arrays

(in-package :fl.matlisp)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; vector operations
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; Sequences are considered as column vectors

(defmethod vlength ((seq sequence))
  "For sequences @function{vlength} is equivalent to @function{length}."
  (length seq))
(defmethod nrows ((seq sequence))
  "For sequences @function{nrows} is equivalent to @function{length}."
  (length seq))
(defmethod ncols ((seq sequence))
  "For sequences @function{ncols} returns identically 1."
  1)

(defmethod element-type ((x vector))
  (array-element-type x))
(defmethod scalar-type ((x vector))
  (array-element-type x))

(defmethod vref ((vec vector) index)
  (aref vec index))
(defmethod (setf vref) (val (vec vector) index)
  (setf (aref vec index) val))
(defmethod vref ((vec array) index)
  (row-major-aref vec index))
(defmethod (setf vref) (val (vec array) index)
  (setf (row-major-aref vec index) val))

(defmethod for-each-key (func (x vector))
  (dotimes (i (length x)) (funcall func i)))

(defmethod for-each-entry (func (seq sequence))
  (map nil func seq))

(defmethod for-each-entry (func (array array))
  (dotimes (i (array-total-size array))
    (funcall func (row-major-aref array i))))

(defmethod for-each-key-and-entry (func (x vector))
  (dotimes (i (length x))
    (funcall func i (aref x i))))

(defmethod for-each-key-and-entry (func (array array))
  (dotimes (i (array-total-size array))
    (funcall func i (row-major-aref array i))))

;;; These do not correspond exactly to their name 
(defmethod for-each-entry-of-vec1 (func (x sequence) (y sequence))
  (map nil func x y))
(defmethod for-each-entry-of-vec2 (func (x sequence) (y sequence))
  (map nil func x y))

(defmethod for-each-entry-of-vec1 (func (x array) (y array))
  (array-for-each func x y))
(defmethod for-each-entry-of-vec2 (func (x array) (y array))
  (array-for-each func x y))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Matlisp operations for arrays
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

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
			     (fl.amop:compile-and-eval
			      (subst element-type 'element-type
				     '(lambda ,(mapcar #'car args)
				       (declare (type (simple-array element-type (*)) ,@vec-args))
				       (declare (type element-type ,@number-args))
				       (declare (optimize (speed 3) (safety 0)))
				       ,@body))))
			    ,functions)))))))
	  (funcall (aref ,functions 0) ,@(mapcar #'car args)))))))

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
    (setf (aref x k) (fill-random! (aref x k) s))))

(defmethod scal! (val (lst list))
  "Call @func{scal!} on each entry of @var{lst}.  Note that the result has
to be freshly consed, because @func{scal!} is a function, not a macro."
  (mapcar #'(lambda (entry) (scal! val entry)) lst))

(defmethod scal! (val (array array))
  "Call @func{scal!} on each entry of @var{array}."
  (dotimes (i (array-total-size array) array)
    (setf (row-major-aref array i)
	  (scal! val (row-major-aref array i)))))

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
    (norm #(1.0 2.0) 1)
  ))

;;; (test-arrays)
(fl.tests:adjoin-test 'test-arrays)


