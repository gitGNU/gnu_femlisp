;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; vector-blas.lisp - vector BLAS routines
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

(defmethod find-applicable-macros (symbol classes)
  "Finds applicable BLAS macros defined on SYMBOL for the given classes."
  (let ((definitions (get symbol 'BLAS)))
    (loop for entry in (get symbol 'BLAS)
	  for specifiers = (car entry)
	  when (every #'subtypep classes specifiers)
	  collect entry)))

(defun more-specific-p (spec1 spec2)
  "Returns T if and only if SPEC1 is really more specific than SPEC2."
  (mapcar #'(lambda (class1 class2)
	      (unless (eq class1 class2)
		(cond
		  ((subtypep class1 class2) (return-from more-specific-p T))
		  ((subtypep class2 class1) (return-from more-specific-p NIL)))))
	  spec1 spec2)
  nil)

(defmethod find-most-specific-macro (symbol classes)
  "Finds a most specific BLAS macro for the given classes."
  (let* ((applicable (find-applicable-macros symbol classes))
	 (best (car applicable)))
    (loop for other in (cdr applicable)
	  when (more-specific-p (car other) (car best))
	  (setq best other))
    best))
;;; Problems : not compatible with CLOS priority

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; BLAS macros
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun clear-blas-macros (class)
  (setf (get class 'BLAS-MACROS) ()))

(defun define-blas-macro (class macro-name macro-definition)
  (setf (geta (get class 'BLAS-MACROS) macro-name)
	macro-definition))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; BLAS building-blocks
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun sym-store (matrix-symbol)
  "Returns a symbol of the form <matrix-symbol>-STORE."
  (symconc matrix-symbol "-STORE"))

(defun sym-size (matrix-symbol)
  "Returns a symbol of the form <matrix-symbol>-SIZE."
  (symconc matrix-symbol "-SIZE"))

(defmethod vector-store ((vec vector))
  `(store (x) x))

(defmethod vector-store ((vec <store-vector>))
  `(vector-store (x) `(slot-value ,x 'store)))

(defmethod with-vector-blas-data (vec)
  "Default method valid for vectors with a store."
  `(with-vector-blas-information (vectors &rest body)
    (let* ,(loop for vector in vectors nconcing
		 `((,(sym-store vector) (vector-store ,vector))
		   (,(sym-size vector) (length ,(sym-store vector)))))
      (declare (ignorable
		,@(loop for vector in vectors nconcing
			(list (sym-store vector) (sym-size vector)))))
      ,@body)))

(defmethod vec-for-each-entry (vec)
  '(vec-for-each-entry (loop-vars &rest body)
    (let ((i (gensym "I")))
      `(dotimes (,i ,(sym-size vec))
	(macrolet ((ref (vector) (list 'aref (sym-store vector) ',i)))
	  ,@body)))))

(defmethod element-copy! (vec)
  '(element-copy! (x y) `(setf ,y ,x)))

(defmethod element-scal! (vec)
  '(element-scal! (alpha x) `(setf ,x (* ,alpha ,x))))

(defmethod element-m+! (vec)
  '(element-m+! (x y) `(incf ,y ,x)))

(defmethod element-m-! (vec)
  '(element-m-! (x y) `(decf ,y ,x)))

(defmethod element-m.*! (vec)
  '(element-m.*! (x y) `(setf y (* ,y ,x))))

(defmethod element-m* (vec)
  '(element-m* (x y) `(* ,y ,x)))

(defmethod element-equal (vec)
  '(element-equal (x y) `(= ,y ,x)))

(defmethod element-gemm! (vec)
  '(element-gemm! (alpha x y beta z)
    `(setf ,z (+ (* ,alpha ,x ,y) (* ,beta ,z)))))

(defmethod with-vector-blas-info (vec)
  `(with-vector-blas-info (vectors &body)
  (macrolet (,(vector-store vec)
		,(vector-size vec)
		,(vec-for-each-entry vec)
		,(element-copy! vec)
		,(element-scal! vec)
		,(element-m+! vec)
		,(element-m-! vec)
		,(element-m.*! vec)
		,(element-m* vec)
		,(element-equal vec)
		,(element-gemm! vec)
		,(with-vector-blas-data (vec)))
      (with-vector-blas-data
    
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

(define-blas-macro 'standard-matrix 'matrix-loop
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
;;; define-blas-template
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun specialized-method-code (name template-args actual-args body)
  (flet ((matrix-p (template-arg)
	   (and (consp template-arg)
		(eq 'standard-matrix (second template-arg))))
	 (number-p (template-arg)
	   (and (consp template-arg)
		(eq 'number (second template-arg)))))
    (let* ((pos (or (position-if #'(lambda (x)
				     (member x '(&optional &key &allow-other-keys)))
				 template-args)
		  (length template-args)))
	   (primary-args (subseq template-args 0 pos))
	   (rest-args (subseq template-args pos)))
      (let (matrices number-args numbers specialized-class element-type)
	(loop for template-arg in primary-args
	      and actual-arg in actual-args
	      for arg = (if (consp template-arg)
			    (car template-arg)
			    template-arg)
	      do
	      (cond
		((matrix-p template-arg)
		 (cond (specialized-class
			(unless (eq specialized-class (class-of actual-arg))
			  (break)
			  (error "Template depends on different classes.")))
		       (t
			(setq specialized-class (class-of actual-arg))
			(setq element-type (element-type actual-arg))))
		 (push arg matrices))
		((number-p template-arg)
		 (push arg number-args)
		 (push actual-arg numbers))
		;; otherwise: do nothing
		))
	(dolist (number numbers)
	  (unless (subtypep (type-of number) element-type)
	    (error "Type of number does not fit with element-type.")))
	`(defmethod
	  ,name
	  (,@(loop for arg in primary-args collect
		   (if (matrix-p arg)
		       `(,(car arg) ,(class-name specialized-class))
		       arg))
	   ,@rest-args)
	  (declare (type ,element-type ,@number-args))
	  (with-blas-information (,@matrices)
	    ,(subst element-type 'element-type
		    `(macrolet ,(mapcar #'cdr (get 'standard-matrix 'BLAS-MACROS))
		      (let ()
			,@body)))))))))

(eval-when (:compile-toplevel :load-toplevel :execute)
  (defun remove-subclass-methods (gf template-args)
    "Removes all methods dispatching on subclasses of the template
arguments."
    (loop for method in (copy-seq (pcl:generic-function-methods gf))
	  when (every #'subtypep
		      (pcl::method-specializers method)
		      (mapcar #'(lambda (arg)
				  (if (consp arg)
				      (second arg)
				      T))
			      template-args))
	  do (remove-method gf method)))

  (defun dispatcher-code (name template-args body)
    "Generates a method which generates code for a type-specialized method."
    (whereas ((gf (and (fboundp name) (symbol-function name))))
      (remove-subclass-methods gf template-args))
    (let ((actual-args (gensym "ACTUAL-ARGS")))
      `(defmethod ,name ,template-args
	,@(when (stringp (car body)) (list (car body)))
	(let ((,actual-args
	       (list ,@(loop for arg in template-args
			     unless (member arg '(&optional))
			     collect (if (consp arg) (car arg) arg)
			     do (assert (not (member arg '(&key &allow-other-keys))))))))
	  ;; define specialized method
	  (let ((method-source
		 (specialized-method-code
		  ',name ',template-args ,actual-args
		  ',(if (stringp (car body)) (cdr body) body))))
	    (dbg-when :blas "Generated method: ~%~A~%" method-source)
	    (eval method-source))
	  ;; retry call
	  (apply ',name ,actual-args))))))

(defmacro define-blas-template (name args &body body)
  (dispatcher-code name args body))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Implementation of the BLAS for the standard-matrix class
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; Copying

(define-blas-template copy! ((x standard-matrix) (y standard-matrix))
  (declare (optimize (speed 3) (debug 0) (safety 0)))
  (assert-standard-matrix-compatibility x y)
  (vec-for-each-entry (x y) (element-copy! (ref x) (ref y)))
  y)

#+(or)
(copy! (eye 2) (eye 2))

;;; BLAS routines involving one matrix

(define-blas-template fill! ((x standard-matrix) (s number))
  (declare (optimize (speed 3) (debug 0) (safety 0)))
  (vec-for-each-entry (x) (element-copy! s (ref x)))
  x)

(define-blas-template fill-random! ((x standard-matrix) (s number))
  (declare (optimize (speed 3) (debug 0) (safety 0)))
  (vec-for-each-entry (x) (element-copy! (random s) (ref x)))
  x)

(define-blas-template scal! ((alpha number) (x standard-matrix))
  (declare (optimize (speed 3) (debug 0) (safety 0)))
  (vec-for-each-entry (x) (element-scal! alpha (ref x)))
  x)

;;; BLAS routines involving two matrices

(define-blas-template axpy! ((alpha number) (x standard-matrix) (y standard-matrix))
  (declare (optimize speed))
  (assert-standard-matrix-compatibility x y)
  (vec-for-each-entry (x y) (element-m+! (element-m* alpha (ref x)) (ref y)))
  y)

(define-blas-template dot ((x standard-matrix) (y standard-matrix))
  "Dot product for standard-matrix.  This is not a perfectly fast routine
due to the condition in the innermost and due to consing when returning the
result.  At bottlenecks, the use of gemm! should be preferred."
  (declare (optimize speed))
  (let ((sum (coerce 0 'element-type)))
    (declare (type element-type sum))
    (vec-for-each-entry (x y)
       (element-m+! (element-m* (ref x) (ref y)) sum))
    sum))

(define-blas-template mequalp ((x standard-matrix) (y standard-matrix))
  "Dot product for standard-matrix.  This is not a perfectly fast routine
due to the condition in the innermost and due to consing when returning the
result.  At bottlenecks, the use of gemm! should be preferred."
  (or (and (or (zerop x-nrows) (zerop x-ncols))
	   (or (zerop x-nrows) (zerop x-ncols)))
      (when (and (= x-nrows y-nrows) (= x-ncols y-ncols))
	(vec-for-each-entry (x y)
          (unless (element-equal (ref x) (ref y))
	    (return-from mequalp nil)))
	t)))

(define-blas-template dot-abs ((x standard-matrix) (y standard-matrix))
  (declare (optimize speed))
  (let ((sum (coerce 0 'element-type)))
    (declare (type element-type sum))
    (vec-for-each-entry (x y)
       (element-m+! (abs (element-m* (ref x) (ref y))) sum))
    sum))

(define-blas-template m+! ((x standard-matrix) (y standard-matrix))
  (declare (optimize speed))
  (assert-standard-matrix-compatibility x y)
  (vec-for-each-entry (x y) (element-m+! (ref x) (ref y)))
  y)

#+(or) (m+! #m((1.0)) #m((1.0)))

(define-blas-template m.*! ((x standard-matrix) (y standard-matrix))
  (declare (optimize speed))
  (assert-standard-matrix-compatibility x y)
  (vec-for-each-entry (x y) (element-m.*! (ref x) (ref y)))
  y)

;;; Matrix-matrix multiplication
(defun generate-standard-matrix-gemm!-template (job)
  "Generates the GEMM-XX! routine defined by JOB."
  (declare (optimize speed))
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
	   (let ((sum (element-scal! beta (ref z))))
	     (declare (type element-type sum))
	     (matrix-loop ((x ,x-index-2) (y ,y-index-1))
	       (element-m+! (element-m* alpha (element-m* (ref x) (ref y))) sum))
	     (setf (ref z) sum))))
       z))))
;;; activate all of the GEMM-XX! routines
(mapc #'generate-standard-matrix-gemm!-template '(:nn :nt :tn :tt))

#+(or)
(let* ((x (make-real-matrix #2a((2.0 1.0) (1.0 2.0))))
       (y (make-real-matrix #2a((2.0 1.0) (1.0 2.0))))
       (z (make-real-matrix 2 2)))
     (gemm-nn! 1.0 x y 0.0 z))

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

(defun typed-vector-generator (type)
  #'(lambda (n)
      (make-instance (standard-matrix type) :nrows n :ncols 1)))

(defun test-blas-1 (fsym size &optional
		    (genvec (typed-vector-generator 'double-float)))
  (let ((x (funcall genvec size))
	(y (funcall genvec size))
	(fn (symbol-function fsym)))
    (format
     t "~A-~D: ~$ MFLOPS~%" fsym size
     (loop with after = 0
	   for before = (get-internal-run-time) then after
	   and count of-type fixnum = 1 then (* count 2)
	   do
	   (loop repeat count do (funcall fn x y))
	   (setq after (get-internal-run-time))
	   (when (> (/ (- after before) internal-time-units-per-second)
		    fl.utilities::*mflop-delta*)
	       (return (/ (* 2 size count internal-time-units-per-second)
			  (* 1e6 (- after before)))))))))

;;;; Testing

(defun test-standard-matrix-blas ()
  (scal! 0.5 #m((1.0 2.0)))
  (copy! (eye 2) (eye 2))
  (test-blas-1 'dot 2 (typed-vector-generator '(complex double-float)))
  (test-blas-1 'dot 1 (typed-vector-generator '(complex double-float)))
  (test-blas-1 'm+! 20 (typed-vector-generator 'double-float))
  (test-blas-1 'm+! 1 (typed-vector-generator 'single-float))
  (test-blas-1 'm.*! 1 (typed-vector-generator 'double-float))
  (test-blas-1 'dot 1 (typed-vector-generator 'single-float))
  (let* ((x #m((2.0 1.0) (-1.0 2.0)))
	 (y #m((2.0 1.0) (1.0 2.0)))
	 (z (make-real-matrix 2 2)))
    (scal 2.0 x)
    (m.* x y)
    (m* x y)
    (m+ x y)
    (axpy -0.5 x y)
    (gemm 1.0 x y 0.0 z)
    (transpose x))
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
   )

;;; (test-standard-matrix-blas)
(fl.tests:adjoin-test 'test-standard-matrix-blas)