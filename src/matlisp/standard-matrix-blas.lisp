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
;;; BLAS operation counters
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defparameter *blas2-operation-count* nil
  "Counter for BLAS level 2 operations")
(defparameter *blas3-operation-count* nil
  "Counter for BLAS level 3 operations")

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; BLAS building-blocks for standard-matrix
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

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
		     `((,(symbol-store matrix) (slot-value ,matrix 'store))
		       (,(msym-nrows matrix) (slot-value ,matrix 'nrows))
		       (,(msym-ncols matrix) (slot-value ,matrix 'ncols))
		       (,(msym-pos matrix) 0))))
	(declare (ignorable
		  ,@(loop for matrix in matrices nconcing
			  (list (symbol-store matrix) (msym-nrows matrix)
				(msym-ncols matrix) (msym-pos matrix)))))
	(declare
	 (type (simple-array element-type (*)) ,@(mapcar #'symbol-store matrices))
	 (type fixnum ,@(loop for matrix in matrices nconcing
			      (list (msym-nrows matrix)
				    (msym-ncols matrix) (msym-pos matrix)))))
	(macrolet ((ref (matrix) `(aref ,(symbol-store matrix) ,(msym-pos matrix))))
	  (let ()
	    ,@body)))))

(defun matrix-loop-expansion (bindings body)
  (let* ((sym1 (car (first bindings)))
	 (type1 (second (first bindings)))
	 ;;(off1 (third (first bindings)))
	 (end1 (fourth (first bindings)))
	 (end (gensym "END")))
    `(let ((,end ,(ecase type1
			 (:row-index
			  `(the fixnum
                                (+ ,(msym-pos sym1)
                                   ,(or end1 (msym-nrows sym1)))))
			 (:col-index
			  `(the fixnum
                                (* ,(msym-nrows sym1)
                                   ,(or end1 (msym-ncols sym1))))))))
      (do ,(loop
	    for (sym type offset) in bindings
	    collecting
	    `(,(msym-pos sym)
	      ,(if offset
		   `(the fixnum
                         (+ ,(msym-pos sym)
                            ,(ecase type
                               (:row-index offset)
                               (:col-index `(the fixnum
                                                 (* ,(msym-nrows sym) ,offset))))))
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

(define-blas-macro 'standard-matrix
    '(assert-vector-operation-compatibility (x y)
      `(unless (and (= ,(msym-nrows x) ,(msym-nrows y))
		(= ,(msym-ncols x) ,(msym-ncols y)))
	(error "Matrices do not have the same format nrows*ncols."))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Implementation of the BLAS for the standard-matrix class
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;; The following routines could be replaced by the BLAS for store-vector.
;;; Unfortunately, this would mean that no matrix compatibility check is
;;; done anymore.  The perfect solution would be to add the check as a
;;; :before method as indicated below.  However, the performance
;;; implications are unclear (e.g., Allegro performance drops by a factor
;;; 20 for copy! on size 1 matrices).  This has to be checked in the
;;; future.

(declaim (inline assert-same-size))
(defun assert-same-size (x y)
  (unless (and (= (nrows x) (nrows y)) (= (ncols x) (ncols y)))
    (error "Matrices do not have the same format nrows*ncols.")))

#+(or)
(defmethod copy! :before ((x standard-matrix) (y standard-matrix))
  (assert-same-size x y))

(defmacro define-blas2-template (name args &body body)
  `(define-blas-template ,name ,args
     #+(or)
     (when *blas2-operation-count*
      (incf *blas2-operation-count* (* (nrows x) (ncols x))))
    ,@body))

(define-blas2-template copy! ((x standard-matrix) (y standard-matrix))
  (vec-for-each-entry (x y) (setf (ref y) (ref x)))
  y)

(define-blas2-template axpy! ((alpha number) (x standard-matrix) (y standard-matrix))
  (vec-for-each-entry (x y) (element-m+! (element-m* alpha (ref x)) (ref y)))
  y)

(define-blas2-template dot ((x standard-matrix) (y standard-matrix))
  (let ((sum (coerce 0 'element-type)))
    (declare (type element-type sum))
    (vec-for-each-entry (x y)
       (element-m+! (element-m* (ref x) (ref y)) sum))
    sum))

(define-blas2-template mequalp ((x standard-matrix) (y standard-matrix))
  (vec-for-each-entry (x y)
     (unless (element-equal (ref x) (ref y))
       (return-from mequalp nil)))
  t)

(define-blas2-template dot-abs ((x standard-matrix) (y standard-matrix))
  (let ((sum (coerce 0 'element-type)))
    (declare (type element-type sum))
    (vec-for-each-entry (x y)
       (element-m+! (abs (element-m* (ref x) (ref y))) sum))
    sum))

(define-blas2-template m+! ((x standard-matrix) (y standard-matrix))
  (vec-for-each-entry (x y) (element-m+! (ref x) (ref y)))
  y)

(define-blas2-template m.*! ((x standard-matrix) (y standard-matrix))
  (vec-for-each-entry (x y) (element-m.*! (ref x) (ref y)))
  y)

;;; Matrix-matrix multiplication

(defmacro conditional-compile (condition form1 form2)
  (aif (eval condition) form1 form2))

(defmacro if-lapack-function ((function type) form1 form2)
  (let ((lapack-func (aand (cl->lapack-type type nil)
                           (nth-value 1 (lapack function it)))))
    `(let ((it ',lapack-func))
       ,(if lapack-func form1 form2))))

(defun generate-standard-matrix-gemm!-template (job)
  "Generates the GEMM-XX! routine defined by JOB."
  (assert (member job '(:nn :nt :tn :tt)))
  (let ((gemm-job (symconc "GEMM-" (symbol-name job) "!"))
	(x-index-1 :row-index) (x-index-2 :col-index)
	(y-index-1 :row-index) (y-index-2 :col-index)
	(x-length-1 'x-nrows) (x-length-2 'x-ncols)
	(y-length-1 'y-nrows) (y-length-2 'y-ncols))
    (when (member job '(:tn :tt))
      (rotatef x-index-1 x-index-2)
      (rotatef x-length-1 x-length-2))
    (when (member job '(:nt :tt))
      (rotatef y-index-1 y-index-2)
      (rotatef y-length-1 y-length-2))
    (eval
     `(define-blas-template ,gemm-job
       ((alpha number) (x standard-matrix) (y standard-matrix) (beta number) (z standard-matrix))
       (with-blas-data (x y z)
	 (unless (and (= ,x-length-1 z-nrows)
                      (= ,x-length-2 ,y-length-1)
                      (= ,y-length-2 z-ncols))
	   (error "Size of arguments does not fit for matrix-matrix multiplication."))
         (when *blas3-operation-count*
           (incf *blas3-operation-count* (* z-nrows z-ncols ,y-length-1)))
         (when (and (plusp z-nrows) (plusp z-ncols))
           (conditional-compile
            (cl->lapack-type 'element-type nil)
            (call-lapack (load-time-value (lapack "gemm" 'element-type))
                         ,(ecase job ((:nn :nt) "N") ((:tn :tt) "T"))
                         ,(ecase job ((:nn :tn) "N") ((:nt :tt) "T"))
                         z-nrows z-ncols ,x-length-2
                         alpha x-store x-nrows y-store y-nrows
                         beta z-store z-nrows)
            (matrix-loop ((z :col-index) (y ,y-index-2))
                         (matrix-loop ((z :row-index) (x ,x-index-1))
                                      (let ((sum (coerce 0 'element-type)))
                                        (declare (type element-type sum))
                                        (matrix-loop ((x ,x-index-2) (y ,y-index-1))
                                                     (element-m+! (element-m* (ref x) (ref y)) sum))
                                        (setf (ref z) (element-m+ (element-m* alpha sum)
                                                                  (element-m* beta (ref z))))))))))
       z))))

;; (gemm-nn! 1.0 (ones 2) (ones 2) 1.0 (zeros 2))
;;; activate all of the GEMM-XX! routines
(mapc #'generate-standard-matrix-gemm!-template '(:nn :nt :tn :tt))

(define-blas-template transpose! ((x standard-matrix) (y standard-matrix))
  (with-blas-data (x y)
    (unless (and (= x-nrows y-ncols) (= x-ncols y-nrows))
      (error "Size of arguments does not fit for destructive transposition."))
    (matrix-loop ((x :col-index) (y :row-index))
		 (matrix-loop ((x :row-index) (y :col-index))
			      (setf (ref y) (ref x)))))
  y)

(define-blas-template minject! ((x standard-matrix) (y standard-matrix)
				row-off col-off)
  (declare (type fixnum row-off col-off))
  (with-blas-data (x y)
    (unless (and (<= 0 row-off) (<= 0 col-off)
		 (<= (+ row-off x-nrows) y-nrows)
		 (<= (+ col-off x-ncols) y-ncols))
      (error "Illegal arguments."))
    (matrix-loop ((x :col-index) (y :col-index col-off))
		 (matrix-loop ((x :row-index) (y :row-index row-off))
			      (setf (ref y) (ref x)))))
  y)

(defmacro define-matrix-matrix-operation (name operation)
  `(define-blas-template ,name ((x standard-matrix) (y standard-matrix)
                                y-row-off y-col-off
                                x-row-off x-col-off x-row-end x-col-end)
     (declare (type fixnum y-row-off y-col-off x-row-off x-col-off x-row-end x-col-end))
     (with-blas-data (x y)
       (let ((x-m (- x-row-end x-row-off))
             (x-n (- x-col-end x-col-off)))
         (declare (type fixnum x-m x-n))
         (unless (and (<= 0 y-row-off) (<= 0 y-col-off)
                      (< y-row-off y-nrows) (< y-col-off y-ncols)
                      (<= 0 x-row-off x-row-end) (<= 0 x-col-off x-col-end)
                      (<= x-row-end x-nrows) (<= x-col-end x-ncols)
                      (<= (+ y-row-off x-m) y-nrows)
                      (<= (+ y-col-off x-n) y-ncols))
           (error "Illegal arguments."))
         #+(or)
         (when *blas2-operation-count*
           (incf *blas2-operation-count* (* x-m x-n)))
         (matrix-loop ((x :col-index x-col-off x-col-end)
                       (y :col-index y-col-off))
                      (matrix-loop ((x :row-index x-row-off x-row-end)
                                    (y :row-index y-row-off))
                                   ,operation)))
       y)))

(define-matrix-matrix-operation extended-minject! (setf (ref y) (ref x)))
(define-matrix-matrix-operation matop-x->y! (setf (ref y) (ref x)))
(define-matrix-matrix-operation matop-x<-y! (setf (ref x) (ref y)))
(define-matrix-matrix-operation matop-y+=x! (incf (ref y) (ref x)))
(define-matrix-matrix-operation matop-y-=x! (decf (ref y) (ref x)))

(define-blas-template mextract! ((x standard-matrix) (y standard-matrix) row-off col-off)
  (declare (type fixnum row-off col-off))
  (with-blas-data (x y)
    (unless (and (<= 0 row-off) (<= 0 col-off)
		 (<= (+ row-off x-nrows) y-nrows)
		 (<= (+ col-off x-ncols) y-ncols))
      (error "Illegal arguments."))
    (matrix-loop ((x :col-index) (y :col-index col-off))
		 (matrix-loop ((x :row-index) (y :row-index row-off))
			      (setf (ref x) (ref y)))))
  x)

(define-blas-template extended-mclear! ((x standard-matrix)
                                        x-row-off x-col-off x-row-end x-col-end)
  (declare (type fixnum x-row-off x-col-off x-row-end x-col-end))
  (with-blas-data (x)
    (unless (and (<= 0 x-row-off x-row-end) (<= 0 x-col-off x-col-end)
                 (<= x-row-end x-nrows) (<= x-col-end x-ncols))
      (error "Illegal arguments."))
    (let ((zero (coerce 0 'element-type)))
      (matrix-loop ((x :col-index x-col-off x-col-end))
                   (matrix-loop ((x :row-index x-row-off x-row-end))
                                (setf (ref x) zero))))
    x))

(defmethod vector-slice ((mat standard-matrix) offset size)
  (let ((result (make-instance (standard-matrix (element-type mat))
			       :nrows size :ncols (ncols mat))))
    (mextract! result mat offset 0)))

(defmethod matrix-slice ((mat standard-matrix) &key (from-row 0) (from-col 0)
			 (nrows (- (nrows mat) from-row))
			 (ncols (- (ncols mat) from-col)))
  (let ((result (make-instance (standard-matrix (element-type mat))
			       :nrows nrows :ncols ncols)))
    (mextract! result mat from-row from-col)))

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
  (copy! #m(1.0) #m(0.0))
  (copy! (funcall (standard-matrix-generator 'single-float) 1)
	 (funcall (standard-matrix-generator 'single-float) 1))
  (dbg-off :blas)
  (dbg-off :compile)
  (test-blas 'copy! 1 :generator (standard-matrix-generator '(complex double-float)))
  (test-blas 'm+! 100 :generator (standard-matrix-generator 'double-float))
  ;; the following shows very slow performance
  (test-blas 'm+! 7992 :generator (standard-matrix-generator '(complex double-float)))
  (scal! 0.5 #m((1.0 2.0)))
  (copy! (eye 2) (eye 2))
  (test-blas 'dot 100 :generator (standard-matrix-generator '(complex double-float)))
  (test-blas 'dot 1 :generator (standard-matrix-generator '(complex double-float)))
  (time (test-blas 'm* 100 :generator 'eye
                           :flop-calculator (lambda (n) (* 2 n n n))))

  ;; performance measurement of BLAS routine GEMM!
  (loop for k from 4 upto 11
        collect
        ;; measure/calculate performance as GFLOPs
        (let* ((n (expt 2 k))
               (repetitions (floor (/ 1d10 (expt n 3))))
               (A (ones n)) (B (ones n)) (C (zeros n)))
          (print n) (print repetitions)
          (let ((time (fl.utilities::measure-time
                       (lambda () (gemm! 1.0 A B 1.0 C))
                       repetitions)))
            (/ (* 2 repetitions (expt n 3)) time 1e9))))
  
  (loop for k = 1 then (* 2 k) until (> k 2000) do
	(fl.matlisp::test-blas
	 'fl.matlisp:m+! k :generator 'fl.matlisp:ones
	 :flop-calculator (lambda (n) (* n n))))
  (test-blas 'm+! 1 :generator (standard-matrix-generator 'single-float))
  (test-blas 'm.*! 1 :generator (standard-matrix-generator 'double-float))
  (test-blas 'dot 100 :generator (standard-matrix-generator 'double-float))
  (test-blas 'dot 1 :generator (standard-matrix-generator 'single-float))
  (let* ((x #m((2.0 1.0) (-1.0 3.0)))
	 (y #m((2.0 1.0) (1.0 2.0)))
	 (z (make-real-matrix 2 2)))
    (scal 2.0 x)
    (m.* x y)
    (m* x y)
    (m+ x y)
1    (axpy -0.5 x y)
    (gemm 1.0 x y 0.0 z)
    (transpose x))
  (transpose #m((1.0 2.0 3.0) (4.0 5.0 6.0)))
  (let ((x (eye 1))
	(y (zeros 1)))
    (copy! x y)
    y)
  (gemm! (/ -27.0) #m((1.5) (0.0)) #m((6.0) (0.0)) (/ 3.0) (eye 2) :nt)
  (gemm! (/ -27.0) #m((1.5) (0.0)) #m((6.0 0.0)) (/ 3.0) (eye 2))
  (join :horizontal #m((1.0 2.0)) #m((3.0 4.0)))
  (join :vertical #m((1.0 2.0)) #m((3.0 4.0)))
  (join :horizontal #m((1.0) (2.0)) #m((3.0) (4.0)))
  (join :vertical #m((1.0) (2.0)) #m((3.0) (4.0)))
  (let* ((x #m((2.0 1.0) (-1.0 2.0)))
	 (y #m(0.0))
	 (z (copy x)))
    (dotimes (i 2)
      (dotimes (j 2)
	(mextract! y x i j)
	(minject! y x i j)))
    (assert (mequalp x z)))
  (let* ((x #m((2.0 1.0 3.0) (6.0 7.0 8.0) (-1.0 -2.0 -3.0)))
	 (y #m((0.1 0.2) (0.3 0.4)))
	 (z (copy x)))
    (dotimes (i 2)
      (dotimes (j 2)
	(mextract! y x i j)
	(minject! y x i j)))
    (assert (mequalp x z)))
  
  (time (let* ((x (laplace-full-matrix 2))
	       (y (copy x)) (z (zeros 2)))
	  (dotimes (i 1000000)
	    (gemm! 1.0 x y 0.0 z))))
  
  (lret ((*blas3-operation-count* 0))
    (time (let* ((n 5) (type 'single-float)
                 (x (zeros n n type))
                 (y (ones n n type))
                 (z (ones n n type)))
            (dotimes (i 100000)
              (gemm! 1.0f0 x y 0.0f0 z)))))
    

  ;;; something special - used in matheum web pages
  (let ((A (make-instance (standard-matrix 'integer) :content #2a((2 3) (4 5)))))
    (m* A A))
  )

;;; (fl.matlisp::test-standard-matrix-blas)
(fl.tests:adjoin-test 'test-standard-matrix-blas)
