;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; matlisp.lisp - Matlisp extensions
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
;;; matlisp corrections
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmethod print-element ((matrix real-matrix) element stream)
  (format stream "~13,9,,,,,'Eg" element))

(setq matlisp::*print-matrix* 5)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; matlisp enhancements
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; from matlisp
(defun make-column-vector (seq)
  (let* ((n (length seq))
	 (store (make-array n :element-type 'real-matrix-element-type)))
    (declare (type fixnum n))
    (dotimes (k n)
      (declare (type fixnum k))
      (setf (aref store k) (coerce (elt seq k) 'real-matrix-element-type)))
    (make-instance 'real-matrix :nrows n :ncols 1 :store store)))

(defun make-row-vector (seq)
  (let* ((n (length seq))
	 (store (make-array n :element-type 'real-matrix-element-type)))
    (declare (type fixnum n))
    (dotimes (k n)
      (declare (type fixnum k))
      (setf (aref store k) (coerce (elt seq k) 'real-matrix-element-type)))
    (make-instance 'real-matrix :nrows 1 :ncols n :store store)))

(definline double-vec->matlisp-column-vector (dv)
  (assert (typep dv 'double-vec))
  (make-instance 'real-matrix :nrows (length dv) :ncols 1 :store dv))

(definline ensure-matlisp (vec &optional (type :column))
  (etypecase vec
      (vector (ecase type
		(:column (make-column-vector vec))
		(:row (make-row-vector vec))))
      (standard-matrix vec)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; vector operations for matlisp matrices
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(definline fortran-matrix-indexing (row col nrows ncols)
  (declare (optimize (speed 3) (safety 0)))
  (declare (type positive-fixnum row col nrows ncols))
  (declare (ignore ncols))
  (the fixnum (+ row (the fixnum (* col nrows)))))

(with-fast-clos ()
(defmethod** mat-ref ((mat real-matrix) i j)
  (declare (optimize (speed 3) (safety 0)))
  (declare (type real-matrix mat))
  (declare (type positive-fixnum i j))
  (let ((store (matlisp::store mat))
	(n (matlisp:nrows mat))
	(m (matlisp:ncols mat)))
    (declare (type double-vec store) (type positive-fixnum n m))
    (aref store (fortran-matrix-indexing i j n m))))
  
(defmethod* (setf mat-ref) ((val double-float) (mat real-matrix) i j)
  (declare (optimize (speed 3) (safety 0)))
  (let ((store (matlisp::store mat))
	(n (matlisp:nrows mat))
	(m (matlisp:ncols mat)))
    (declare (type double-vec store) (type positive-fixnum n m))
    (setf (aref store (fortran-matrix-indexing i j n m))
	  val)))

(defmethod** vec-ref ((mat real-matrix) i)
  "For vector-like access to matlisp matrices."
  (declare (optimize (speed 3) (safety 0)))
  (declare (type positive-fixnum i))
  (let ((store (matlisp::store mat)))
    (declare (type double-vec store))
    (aref store i)))

(defmethod** (setf vec-ref) (val (mat real-matrix) i)
  "For vector-like access to matlisp matrices."
  (declare (optimize (speed 3) (safety 0)))
  (declare (type positive-fixnum i))
  (let ((store (matlisp::store mat)))
    (declare (type double-vec store))
    (setf (aref store i) val)))
)  ; with-fast-clos


(defmethod multiplicity ((mat real-matrix))
  (ncols mat))

(defmethod make-analog ((mat real-matrix))
  (make-real-matrix (nrows mat) (ncols mat)))

(defmethod make-analog ((arr array))
  (make-array (array-dimensions arr) :element-type (array-element-type arr)))

;;; define a row for the abstract looping interface
(defmethod** matrix-row ((mat real-matrix) i)
  (cons mat i))

(defmethod** for-each-row-key ((fn function) (mat real-matrix))
  "Loop through row keys."
  (dotimes (i (nrows mat))
    (funcall fn i)))
(defmethod** for-each-row ((fn function) (mat real-matrix))
  (for-each-row-key #'(lambda (i) (matrix-row mat i)) mat))
(defmethod** for-each-key-in-row ((fn function) (mat real-matrix) row-key)
  (declare (ignore row-key))
  (dotimes (i (ncols mat))
    (funcall fn i)))
;;; better: (for-each-key-in-row ((fn function) (row cons)) ...)
(defmethod** for-each-col-key ((fn function) (mat real-matrix))
  " Loop through column keys."
  (dotimes (i (ncols mat))
    (funcall fn i)))

(defmethod for-each-key-and-entry ((func function) (mat real-matrix))
  (dotimes (i (nrows mat))
    (dotimes (j (ncols mat))
      (funcall func i j (mat-ref mat i j)))))
    
;;; Define vector blas commands for matlisp matrices.  The technique is
;;; working simply on the store.
(defun matlisp-blas-caller (proc) (list proc))
(defun matlisp-vector-transformer (obj) (list 'matlisp::store obj))
(initialize-vector-blas-methods
 'standard-matrix 'double-float #'matlisp-blas-caller
 :vector-transformer #'matlisp-vector-transformer)

;;; matrix-vector multiplication
(defmethod x<-Ay (x (A standard-matrix) y)
  (gemm! 1.0d0 A (ensure-matlisp y) 0.0d0 (ensure-matlisp x))
  x)
(defmethod x+=Ay (x (A standard-matrix) y)
  (gemm! 1.0d0 A (ensure-matlisp y) 1.0d0 (ensure-matlisp x))
  x)
(defmethod x-=Ay (x (A standard-matrix) y)
  (gemm! -1.0d0 A (ensure-matlisp y) 1.0d0 (ensure-matlisp x))
  x)
(defmethod x+=s*Ay ((s double-float) x (A standard-matrix) y)
  (gemm! s A (ensure-matlisp y) 1.0d0 (ensure-matlisp x))
  x)

(defmethod x<-random ((x standard-matrix) s)
  (x<-random (matlisp::store x) s))

(progn   ;alternative with explicit multiplication
  (definline mat-row-vec-loop (A nrows ncols i y)
    (declare (type double-vec y A))
    (declare (type fixnum nrows ncols i))
    (assert (= (length y) ncols))
    (let ((sum 0.0d0))
      (declare (optimize (speed 3) (safety 0)))
      (declare (type double-float sum))
      (dotimes (j ncols)
	(declare (fixnum j))
	(incf sum (* (aref A (fortran-matrix-indexing i j nrows ncols)) (aref y j))))
      sum))
  
  (defmethod x<-Ay ((x array) (A real-matrix) (y array))
    (dotimes (i (nrows A) x)
      (setf (vec-ref x i) (mat-row-vec-loop (matlisp::store A) (nrows A) (ncols A) i y))))
  (defmethod x+=Ay ((x array) (A real-matrix) (y array))
    (dotimes (i (nrows A) x)
      (incf (vec-ref x i) (mat-row-vec-loop (matlisp::store A) (nrows A) (ncols A) i y))))
  (defmethod x-=Ay ((x array) (A real-matrix) (y array))
    (dotimes (i (nrows A) x)
      (decf (vec-ref x i) (mat-row-vec-loop (matlisp::store A) (nrows A) (ncols A) i y))))
  )

(defmethod total-entries ((mat standard-matrix))
   (* (nrows mat) (ncols mat)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Matlisp operations for arrays
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmethod axpy (alpha (x number) (y number))
  (+ y (* alpha x)))

(defmethod axpy (alpha x y)
  (let ((result (copy y)))
    (x+=s*y result alpha x)))

(defmethod axpy! (alpha x y)
  (x+=s*y y alpha x))

(defmethod dot ((vec1 array) (vec2 array) &optional conjugate-p)
  (let ((sum 0.0d0))
    (array-for-each
     #'(lambda (entry1 entry2)
	 (incf sum (dot entry1 entry2 conjugate-p)))
     vec1 vec2)
    sum))

(defmethod dot-abs ((mat1 real-matrix) (mat2 real-matrix) &optional conjugate-p)
  (declare (ignore conjugate-p))
  (let ((sum 0.0d0))
    (dotimes (i (nrows mat1))
      (dotimes (j (ncols mat1))
	(incf sum (abs (* (mat-ref mat1 i j) (mat-ref mat2 i j))))))
    sum))

(defmethod m+ ((x number) (y array))
  (+ x (aref y 0)))

(defmethod m+ ((x array) (y number))
  (+ (aref x 0) y))

(defmethod m+ ((vec1 array) (vec2 array))
  (map (type-of vec1) #'+ vec1 vec2))

(defmethod m+! ((vec1 array) (vec2 array))
  (dotimes (i (length vec1))
    (incf (aref vec2 i) (aref vec1 i))))

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

(defmethod m- ((x number) (y number)) (- x y))

(defmethod m- ((vec1 array) (vec2 array))
  (if (vectorp vec1)
      (map (type-of vec1) #'- vec1 vec2)
      (let ((result (make-analog vec1)))
	(assert (= (array-rank result) 2))
	(dotimes (i (array-dimension result 0))
	  (dotimes (j (array-dimension result 1))
	    (setf (aref result i j) (m- (aref vec1 i j) (aref vec2 i j)))))
	result)))

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

(defmethod gemm! (alpha A y beta x &optional (job :nn))
  (gemm! alpha (ensure-matlisp A) (ensure-matlisp y) beta (ensure-matlisp x) job))

(defmethod gemm (alpha A y beta x &optional (job :nn))
  (gemm alpha (ensure-matlisp A) (ensure-matlisp y) beta (ensure-matlisp x) job))

(defmethod m* ((x standard-matrix) (y array))
  (matlisp::store (m* x (make-column-vector y))))

(defun vector->column-array (x)
  (let ((result (make-array (list (length x) 1))))
    (dotimes (i (length x) result)
      (setf (aref result i 0) (aref x i)))))
  
(defmethod m* ((A array) (B array))
  (let ((A (if (= 1 (array-rank A)) (vector->column-array A) A))
	(B (if (= 1 (array-rank B)) (vector->column-array B) B)))
    (assert (= (array-dimension A 1) (array-dimension B 0)))
    (let ((C (make-array (list (array-dimension A 0) (array-dimension B 1)))))
      (gemm! 1.0 A B 0.0 C)
      C)))

(defmethod m*! ((x standard-matrix) y)
  (m*! x (ensure-matlisp y))
  y)

(defmethod m*! ((x standard-matrix) (y array))
  (m*! x 
       (make-instance 'real-matrix :nrows (length y) :ncols 1 :store y))
  y)


(definline m/-check (a b)
  (if b
      (typecase b
	(number t)
	(standard-matrix (unless (= (nrows b) (nrows a))
			   (error "dimensions of A,B given to M/ do not match")))
	(vector (unless (= (length b) (nrows a))
		  (error "dimensions of A,B given to M/ do not match")))
	(t (error "argument B given to M/ is not a matrix, vector or number")))
      (unless (square-matrix-p a)
	(error "argument A given to M/ is not a square matrix"))))

(definline m/-solve (a b)
  (multiple-value-bind (x ipiv f info)
      (gesv a b)
    (declare (ignore ipiv f))
    (if (numberp info)		     
	(error "argument A given to M/ is singular to working machine precision")
	x)))
  
(defmethod m/ ((x number) &optional y)
  (if y (/ y x) (/ x)))

(defmethod m/ :before ((a standard-matrix) &optional b)
  (m/-check a b))

(defmethod m/ ((a standard-matrix) &optional b)
  (if b
      (typecase b
	(number (m./ a b))
	(vector (matlisp::store (m/-solve a (double-vec->matlisp-column-vector b))))
	(standard-matrix (m/-solve a b)))
      (m/-solve a (eye (nrows a)))))

(defmethod m/! :before ((a standard-matrix) &optional b)
  (m/-check a b))

(definline m/!-solve (a b)
  (multiple-value-bind (x ipiv f info)
      (gesv! a b)
    (declare (ignore ipiv f))
    (if (numberp info)		     
	(error "argument A given to M/ is singular to working machine precision")
	x)))

(defmethod m/! ((a standard-matrix) &optional b)
  (if b
      (typecase b
	(number (m./! a b))
	(vector (matlisp::store (m/!-solve a (double-vec->matlisp-column-vector b))))
	(standard-matrix (m/!-solve a b)))
      (m/!-solve (copy a) (eye (nrows a)))))

;;; scal for numbers does exist already in matlisp

(defmethod norm ((vec array) &optional (p 2))
  (case p
    ((:inf) (loop for x across vec maximize (norm x p)))
    (otherwise (expt (loop for x across vec sum (expt (norm x p) p)) (/ 1 p)))))

;;(defmethod ncols ((x number)) 1)
;;(defmethod nrows ((x number)) 1)

;;; recursive definition of scaling for block vectors 

(defmethod scal! (s x)
  (for-each-entry (curry #'scal! s) x))

(defmethod scal (s x)
  (let ((result (copy x)))
    (scal! s result)
    result))

(defmethod scal ((val number) (arr array))
  (if (vectorp arr)
      (map (type-of arr) #'(lambda (x) (* val x)) arr)
      (call-next-method)))
  
(defmethod scal! ((val number) (arr array))
  (if (vectorp arr)
      (dotimes (i (length arr) arr)
	(setf (aref arr i) (* (aref arr i) val)))
      (call-next-method)))

(defmethod matlisp::matrix-ref-1d ((vec array) i) (aref vec i))
(defmethod matlisp::matrix-ref-2d ((vec array) i j) (aref vec i j))

(defmethod (setf matlisp::matrix-ref-1d) (value (vec array) i)
  (setf (aref vec i) value))
(defmethod (setf matlisp::matrix-ref-2d) (value (vec array) i j)
  (setf (aref vec i j) value))

(defun matlisp::unit-vector (dim i)
  (let ((vec (make-float-matrix dim 1)))
    (setf (vec-ref vec i) 1.0d0)
    vec))

(defmethod print-cell ((matrix real-matrix)
			  cell
			  stream)
  (format stream "~18,10,,,'*,,'EE" cell))

;;; copying of arrays

(defmethod copy ((mat array))
  (cond
    ((vectorp mat) (copy-seq mat))
    ((= (array-rank mat) 2)
     (let ((result (make-analog mat)))
	 (dotimes (i (array-dimension mat 0))
	   (dotimes (j (array-dimension mat 1))
	     (setf (aref result i j) (copy (aref mat i j)))))
	 result))
    (t (error "only rank2 operation implemented here"))))

;;; new matlisp commands

(defun det-from-lr (lr pivot)
  "This routine computes the determinant using a given LR decomposition."
  (loop with prod of-type double-float = 1.0d0
	for i below (number-of-rows lr) do
	(setq prod (* prod (mat-ref lr i i)))
	(unless (= (aref pivot i) (1+ i)) (setq prod (- prod)))
	finally (return prod)))

(defun det (mat)
  "Example of use: (det [[-1 2]' [2 3]'])  -> -7.0"
  (declare (values double-float))
  (if (or (zerop (nrows mat)) (zerop (ncols mat)))
      1.0d0  ; yields correct value for volume in the vertex case
      (multiple-value-bind (lr pivot)
	  (getrf! (copy mat))
	(det-from-lr lr pivot))))

(defun area-of-span (mat)
  "Computes the area/volume spanned by k vectors in Rn given as columns of
the argument mat.  For k=n, this is (abs (det mat))."
  (sqrt
   (loop for set in (k-subsets (range< 0 (nrows mat)) (ncols mat))
	 for det = (det (submatrix mat :row-indices set))
	 summing (* det det))))

(defgeneric getrs! (ldu ipiv result &key trans))
(defgeneric getrs (ldu ipiv result &key trans))

(defmethod mzerop ((x number) &optional (threshold 0.0))
  (<= (abs x) threshold))

(defmethod mzerop ((mat standard-matrix) &optional (threshold 0.0))
  (every #'(lambda (x) (<= (abs x) threshold))
	 (matlisp::store mat)))

(defmethod midentity-p ((x number) &optional (threshold 0.0))
  (<= (abs (- x 1.0)) threshold))

(defmethod midentity-p (mat &optional (threshold 0.0))
  (for-each-key-and-entry
   #'(lambda (rk ck entry)
       (unless (if (eql rk ck)
		   (midentity-p entry threshold)
		   (mzerop entry threshold))
	 (return-from midentity-p (values nil rk ck entry))))
   mat)
  t)

(defmethod submatrix ((mat real-matrix) &key row-indices col-indices)
  (unless row-indices (setq row-indices (range< 0 (nrows mat))))
  (unless col-indices (setq col-indices (range< 0 (ncols mat))))
  (let ((result (make-real-matrix (length row-indices) (length col-indices))))
    (loop for i from 0 and row-ind in row-indices do
      (loop for j from 0 and col-ind in col-indices do
	    (setf (mat-ref result i j) (mat-ref mat row-ind col-ind))))
    result))

(defmethod matrix-slice ((mat real-matrix) &key (from-row 0) (from-col 0) nrows ncols)
  (unless nrows (setq nrows (- (nrows mat) from-row)))
  (unless ncols (setq ncols (- (ncols mat) from-col)))
  (submatrix mat :row-indices (make-set from-row nrows)
	     :col-indices (make-set from-col ncols)))

(defmethod make-row-vector-for ((mat real-matrix) &optional (multiplicity 1))
  (make-float-matrix (nrows mat) multiplicity))

(defmethod make-column-vector-for ((mat real-matrix) &optional (multiplicity 1))
  (make-float-matrix (ncols mat) multiplicity))

(defmethod x<-id ((mat standard-matrix))
  (dotimes (i (nrows mat))
    (dotimes (j (ncols mat))
      (setf (mat-ref mat i j) 0.0d0))
    (setf (mat-ref mat i i) 1.0d0)))

(defmethod clear-row ((mat standard-matrix) (row fixnum)
		      &optional row2)
  (declare (ignore row2))
  (dotimes (k (ncols mat))
    (setf (mat-ref mat row k) 0.0d0)))

(defmethod clear-column ((mat standard-matrix) (col fixnum)
			 &optional col2)
  (declare (ignore col2))
  (dotimes (k (ncols mat))
    (setf (mat-ref mat k col) 0.0d0)))

(defun zero-row? (mat i &optional ignore-diag)
  (dotimes (j (ncols mat) t)
    (unless (or (and ignore-diag (= i j))
		(zerop (mat-ref mat i j)))
      (return nil))))

(defun zero-column? (mat j &optional ignore-diag)
  (dotimes (i (nrows mat) t)
    (unless (or (and ignore-diag (= i j))
		(zerop (mat-ref mat i j)))
      (return nil))))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Testing
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun test-matlisp-vector-combination ()
  ;;(c::info function c::info (quote algebra::vec-ref))
  (let ((x (make-double-vec 2))
	(y (make-double-vec 2 1.0d0))
	(A [[1 2]' [3 4]'])
	(B (make-real-matrix 2 2)))
    (x<-0 x)
    (x+=y x y)
    (x+=s*y x 2.0d0 y)
    (x-=y x y)
    (x-=s*y x 3.0d0 y)
    (assert (= (aref x 0) -1.0d0))
    (x<-0 B)
    (x+=y B A)
    (x+=s*y B 2.0d0 A)
    (x-=y B A)
    (x-=s*y B 3.0d0 A)
    (assert (= (mat-ref B 0 0) -1.0d0))
  ))

;;; (test-matlisp-vector-combination)
(tests::adjoin-femlisp-test 'test-matlisp-vector-combination)

