;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; standard-matrix-blas.lisp - some BLAS routines
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

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; BLAS building-blocks for standard-matrix
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun msym-store (matrix-symbol)
  "Returns a symbol of the form <matrix-symbol>-STORE."
  (symconc matrix-symbol "-STORE"))

(defun msym-pos (matrix-symbol)
  "Returns a symbol of the form <matrix-symbol>-POS."
  (symconc matrix-symbol "-POS"))

(defun msym-nrows (matrix-symbol)
  "Returns a symbol of the form <matrix-symbol>-NROWS."
  (symconc matrix-symbol "-NROWS"))

(defun msym-ncols (matrix-symbol)
  "Returns a symbol of the form <matrix-symbol>-NCOLS."
  (symconc matrix-symbol "-NCOLS"))

(defun msym-size (matrix-symbol)
  "Returns a symbol of the form <matrix-symbol>-SIZE."
  (symconc matrix-symbol "-SIZE"))

(define-blas-macro 'standard-matrix
    '(with-blas-data (matrices &rest body)
      "Sets body in an environment with local variables that access store,
nrows and ncols of the given matrices."
      `(let (,@(loop for matrix in matrices nconcing
		    `((,(msym-store matrix) (slot-value ,matrix 'store))
		      (,(msym-nrows matrix) (slot-value ,matrix 'nrows))
		      (,(msym-ncols matrix) (slot-value ,matrix 'ncols))
		      (,(msym-pos matrix) 0))))
	 (declare (ignorable
		   ,@(loop for matrix in matrices nconcing
			  (list (msym-store matrix) (msym-nrows matrix)
				(msym-ncols matrix) (msym-pos matrix)))))
	 (declare
	  (type (simple-array element-type (*)) ,@(mapcar #'msym-store matrices))
	  (type fixnum ,@(loop for matrix in matrices nconcing
			      (list (msym-nrows matrix)
				    (msym-ncols matrix) (msym-pos matrix)))))
	 (macrolet ((ref (matrix) `(aref ,(msym-store matrix) ,(msym-pos matrix))))
	   ,@body))))

#+(or)
(macroexpand-1 '(with-blas-information (x y) double-float
		 (setf (ref x) (ref y))))
  
(define-blas-macro 'standard-matrix
  '(vec-for-each-entry (loop-vars &rest body)
    (let ((i (gensym "I")))
      `(dotimes (,i (length ,(msym-store (first loop-vars))))
	(macrolet ((ref (matrix) (list 'aref (msym-store matrix) ',i)))
	  ,@body)))))

(defmacro vec-for-each-entry (loop-vars &rest body)
  (let ((i (gensym "I")))
    `(macrolet ((ref (matrix) (list 'aref (msym-store matrix) ',i)))
      (dotimes (,i (length ,(msym-store (first loop-vars))))
	,@body))))

(define-blas-macro 'standard-matrix
  '(element-copy! (x y) `(setf ,y ,x)))

(define-blas-macro 'standard-matrix
  '(element-scal! (alpha x) `(setf ,x (* ,alpha ,x))))

(define-blas-macro 'standard-matrix
  '(element-m+! (x y) `(incf ,y ,x)))

(define-blas-macro 'standard-matrix
  '(element-m+ (x y) `(+ ,y ,x)))

(define-blas-macro 'standard-matrix
  '(element-m-! (x y) `(decf ,y ,x)))

(define-blas-macro 'standard-matrix
  '(element-m.*! (x y) `(setf ,y (* ,y ,x))))

(define-blas-macro 'standard-matrix
  '(element-m* (x y) `(* ,y ,x)))

(define-blas-macro 'standard-matrix
  '(element-equal (x y) `(= ,y ,x)))

(define-blas-macro 'standard-matrix
  '(element-gemm! (alpha x y beta z)
    `(setf ,z (+ (* ,alpha ,x ,y) (* ,beta ,z)))))

(defun matrix-loop-expansion (bindings body)
  (let* ((sym1 (car (first bindings)))
	 (type1 (second (first bindings)))
	 (end (gensym "END")))
    `(let ((,end ,(ecase type1
			 (:row-index
			  `(the fixnum
			    (+ ,(msym-pos sym1) ,(msym-nrows sym1))))
			 (:col-index
			  `(the fixnum
			    (* ,(msym-nrows sym1) ,(msym-ncols sym1)))))))
      (do ,(loop
	    for (sym type offset) in bindings
	    collecting
	    `(,(msym-pos sym)
	      ,(if offset
		   `(+ ,(msym-pos sym) ,offset)
		   (msym-pos sym))
	      (+ ,(msym-pos sym)
	       ,(ecase type
		       (:row-index 1)
		       (:col-index (msym-nrows sym))))))
	  ((>= ,(msym-pos sym1) ,end))
	(declare (type fixnum ,@(mapcar #'(lambda (binding)
					    (msym-pos (car binding)))
					bindings)))
	,@body))))

(define-blas-macro 'standard-matrix
  '(matrix-loop (bindings &rest body)
    (matrix-loop-expansion bindings body)))

;;(declaim (inline standard-matrix-compatible-p))
(defmacro standard-matrix-compatible-p (x y)
  `(and (= ,(msym-nrows x) ,(msym-nrows y))
    (= ,(msym-ncols x) ,(msym-ncols y))))

;;;(declaim (inline assert-standard-matrix-compatibility))
(defmacro assert-standard-matrix-compatibility (x y)
  `(unless (standard-matrix-compatible-p ,x ,y)
    (error "Matrices do not have the same format nrows*ncols.")))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Implementation of the BLAS for the standard-matrix class
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;; >>>>

;;; The following routines could be replaced by the BLAS for store-vector.
;;; Unfortunately, this would mean that no matrix compatibility check is
;;; done anymore.  The perfect solution would be to add the check as a
;;; :before method as indicated here.  However, the performance
;;; implications are not clear.  This has to be checked in the future.

(declaim (inline assert-same-size))
(defun assert-same-size (x y)
  (unless (and (= (nrows x) (nrows y)) (= (ncols x) (ncols y))))
    (error "Matrices do not have the same format nrows*ncols."))

#+(or)
(defmethod copy! :before ((x standard-matrix) (y standard-matrix))
  (assert-same-size x y))

;;;; <<<<

(define-blas-template copy! ((x standard-matrix) (y standard-matrix))
  (declare (optimize (speed 3) (debug 0) (safety 0)))
  (assert-standard-matrix-compatibility x y)
  (vec-for-each-entry (x y) (element-copy! (ref x) (ref y)))
  y)

(define-blas-template axpy! ((alpha number) (x standard-matrix) (y standard-matrix))
  (declare (optimize speed))
  (assert-standard-matrix-compatibility x y)
  (vec-for-each-entry (x y) (element-m+! (element-m* alpha (ref x)) (ref y)))
  y)

(define-blas-template m+! ((x standard-matrix) (y standard-matrix))
  (declare (optimize speed))
  (assert-standard-matrix-compatibility x y)
  (vec-for-each-entry (x y) (element-m+! (ref x) (ref y)))
  y)

(define-blas-template mequalp ((x standard-matrix) (y standard-matrix))
  "Exact equality test for standard-matrix."
  (or (and (or (zerop x-nrows) (zerop x-ncols))
	   (or (zerop x-nrows) (zerop x-ncols)))
      (when (and (= x-nrows y-nrows) (= x-ncols y-ncols))
	(vec-for-each-entry (x y)
          (unless (element-equal (ref x) (ref y))
	    (return-from mequalp nil)))
	t)))

;;; Matrix-matrix multiplication

(defun generate-standard-matrix-gemm!-template (job)
  "Generates the GEMM-XX! routine defined by JOB."
  (assert (member job '(:nn :nt :tn :tt)))
  (let ((gemm-job (symconc "GEMM-" (symbol-name job) "!"))
	(x-index-1 :row-index) (x-index-2 :col-index)
	(y-index-1 :row-index) (y-index-2 :col-index)
	(x-length-1 'x-nrows) (x-length-2 'x-ncols)
	(y-length-1 'y-nrows) (y-length-2 'y-ncols))
    (when (member job '(:tn tt))
      (rotatef x-index-1 x-index-2)
      (rotatef x-length-1 x-length-2))
    (when (member job '(:nt tt))
      (rotatef y-index-1 y-index-2)
      (rotatef y-length-1 y-length-2))
    (eval
     `(define-blas-template ,gemm-job
       ((alpha number) (x standard-matrix) (y standard-matrix) (beta number) (z standard-matrix))
       (declare (optimize (speed 3) (space 0) (debug 0) (safety 0)))
       (unless (and (= ,x-length-1 z-nrows)
		    (= ,x-length-2 ,y-length-1)
		    (= ,y-length-2 z-ncols))
	 (error "Size of arguments does not fit for matrix-matrix multiplication."))
       (matrix-loop ((z :col-index) (y ,y-index-2))
	 (matrix-loop ((z :row-index) (x ,x-index-1))
	   (let ((sum (coerce 0 'element-type)))
	     (declare (type element-type sum))
	     (matrix-loop ((x ,x-index-2) (y ,y-index-1))
	       (element-m+! (element-m* (ref x) (ref y)) sum))
	     (setf (ref z) (element-m+ (element-m* alpha sum)
				       (element-m* beta (ref z)))))))
       z))))

;;; activate all of the GEMM-XX! routines
(mapc #'generate-standard-matrix-gemm!-template '(:nn :nt :tn :tt))

(define-blas-template transpose! ((x standard-matrix) (y standard-matrix))
  (declare (optimize (speed 3) (space 0) (debug 0) (safety 0)))
  (unless (and (= x-nrows y-ncols) (= x-ncols y-nrows))
    (error "Size of arguments does not fit for destructive transposition."))
  (matrix-loop ((x :col-index) (y :row-index))
    (matrix-loop ((x :row-index) (y :col-index))
      (setf (ref y) (ref x))))
  y)

(define-blas-template join-horizontal! ((x standard-matrix) (y standard-matrix)
					(z standard-matrix))
  (declare (optimize (speed 3) (space 0) (debug 0) (safety 0)))
  (unless (and (= x-nrows y-nrows z-nrows) (= z-ncols (+ x-ncols y-ncols)))
    (error "Size of arguments does not fit for horizontal join."))
  (matrix-loop ((x :col-index) (z :col-index))
    (matrix-loop ((x :row-index) (z :row-index))
      (setf (ref z) (ref x))))
  (matrix-loop ((y :col-index) (z :col-index (length x-store)))
    (matrix-loop ((y :row-index) (z :row-index))
      (setf (ref z) (ref y))))
  z)
  
(define-blas-template join-vertical! ((x standard-matrix) (y standard-matrix)
				      (z standard-matrix))
  (declare (optimize (speed 3) (space 0) (debug 0) (safety 0)))
  (unless (and (= x-ncols y-ncols z-ncols) (= z-nrows (+ x-nrows y-nrows)))
    (error "Size of arguments does not fit for vertical join."))
  (matrix-loop ((x :col-index) (y :col-index) (z :col-index))
    (matrix-loop ((x :row-index) (z :row-index))
      (setf (ref z) (ref x)))
    (matrix-loop ((y :row-index) (z :row-index x-nrows))
      (setf (ref z) (ref y))))
  z)

(define-blas-template minject ((x standard-matrix) (y standard-matrix)
			       row-off col-off)
  (declare (type fixnum row-off col-off))
  (declare (optimize (speed 3) (space 0) (debug 0) (safety 0)))
  (unless (and (<= 0 row-off) (<= 0 col-off)
	       (<= (+ row-off x-nrows) y-nrows)
	       (<= (+ col-off x-ncols) y-ncols))
    (error "Illegal arguments."))
  (matrix-loop ((x :col-index) (y :col-index (* col-off y-nrows)))
    (matrix-loop ((x :row-index) (y :row-index row-off))
      (setf (ref y) (ref x))))
  y)

(define-blas-template mextract ((x standard-matrix) (y standard-matrix) row-off col-off)
  (declare (type fixnum row-off col-off))
  (declare (optimize (speed 3) (space 0) (debug 0) (safety 0)))
  (unless (and (<= 0 row-off) (<= 0 col-off)
	       (<= (+ row-off y-nrows) x-nrows)
	       (<= (+ col-off y-ncols) x-ncols))
    (error "Illegal arguments."))
  (matrix-loop ((y :col-index) (x :col-index (* col-off x-nrows)))
    (matrix-loop ((y :row-index) (x :row-index row-off))
      (setf (ref y) (ref x))))
  y)

#+(or)
(time
 (let* ((x (make-real-matrix #2a((2.0 1.0) (1.0 2.0))))
	(y (make-real-matrix #2a((2.0 1.0) (1.0 2.0))))
	(z (make-real-matrix 2 2)))
   (dotimes (i 1000000)
     (gemm-nn! 1.0 x y 0.0 z))))

#+(or)
(time 
 (let* ((type 'single-float)
	(n 5)
	(x (make-instance (standard-matrix type) :nrows n :ncols n))
	(y (make-instance (standard-matrix type) :nrows n :ncols n))
	(z (make-instance (standard-matrix type) :nrows n :ncols n)))
   (dotimes (i 100000)
     (gemm-nn! 1.0f0 x y 0.0f0 z))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Testing
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun standard-matrix-generator (type)
  #'(lambda (n)
      (make-instance (standard-matrix type) :nrows n :ncols 1)))

;;;; Testing
(defun test-standard-matrix-blas ()
  (dbg-on :blas)
  (dbg-on :compile)
  (test-blas 'm+! 7992 :generator (standard-matrix-generator '(complex double-float)))
  (scal! 0.5 #m((1.0 2.0)))
  (copy! (eye 2) (eye 2))
  (test-blas 'dot 2 :generator (standard-matrix-generator '(complex double-float)))
  (test-blas 'dot 1 :generator (standard-matrix-generator '(complex double-float)))
  (time (test-blas 'm* 100 :generator 'eye
	     :flop-calculator (lambda (n) (* 2 n n n))))
  (loop for k = 1 then (* 2 k) until (> k 1000) do
	(fl.matlisp::test-blas
	 'fl.matlisp:m+! k :generator 'fl.matlisp:ones
	 :flop-calculator (lambda (n) (* n n))))
  #+(or)
  (loop for k = 1 then (* 2 k) until (> k 1000) do
	(fl.matlisp::test-blas
	 'matlisp:m+! k :generator 'matlisp:ones
	 :flop-calculator (lambda (n) (* n n))))
  (test-blas 'm+! 1 :generator (standard-matrix-generator 'single-float))
  (test-blas 'm.*! 1 :generator (standard-matrix-generator 'double-float))
  (test-blas 'dot 1 :generator (standard-matrix-generator 'single-float))
  (let* ((x #m((2.0 1.0) (-1.0 3.0)))
	 (y #m((2.0 1.0) (1.0 2.0)))
	 (z (make-real-matrix 2 2)))
    (scal 2.0 x)
    (m.* x y)
    (m* x y)
    (m+ x y)
    (axpy -0.5 x y)
    (gemm 1.0 x y 0.0 z)
    (transpose x))
  (transpose #m((1.0 2.0 3.0) (4.0 5.0 6.0)))
  (let ((x (eye 1))
	(y (zeros 1)))
    (copy! x y)
    y)
  (gemm-nt! (/ -27.0) #m((1.5) (0.0)) #m((6.0) (0.0)) (/ 3.0) (eye 2))
  (gemm-nn! (/ -27.0) #m((1.5) (0.0)) #m((6.0 0.0)) (/ 3.0) (eye 2))
  (join #m((1.0 2.0)) #m((3.0 4.0)))
  (join #m((1.0 2.0)) #m((3.0 4.0)) :vertical)
  (join #m((1.0) (2.0)) #m((3.0) (4.0)))
  (join #m((1.0) (2.0)) #m((3.0) (4.0)) :vertical)
  (mequalp #m() #m(()))
  (let* ((x #m((2.0 1.0) (-1.0 2.0)))
	 (y #m(0.0))
	 (z (copy x)))
    (dotimes (i 2)
      (dotimes (j 2)
	(mextract x y i j)
	(minject y x i j)))
    (assert (mequalp x z)))
  (let* ((x #m((2.0 1.0 3.0) (6.0 7.0 8.0) (-1.0 -2.0 -3.0)))
	 (y #m((0.1 0.2) (0.3 0.4)))
	 (z (copy x)))
    (dotimes (i 2)
      (dotimes (j 2)
	(mextract x y i j)
	(minject y x i j)))
    (assert (mequalp x z)))
  (dbg-off :blas)
  (dbg-off :compile)
   )

;;; (time (fl.matlisp::test-standard-matrix-blas))
(fl.tests:adjoin-test 'test-standard-matrix-blas)
