;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; matrix.lisp - Matrix class and interface
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

(defclass <matrix> (<vector>)
  ()
  (:documentation "General matrix class."))

(defgeneric nrows (mat)
  (:documentation "Number of matrix rows."))
(defgeneric ncols (mat)
  (:documentation "Number of matrix columns."))

(defgeneric mref (A i j)
  (:documentation "Returns the matrix element @code{A[i,j]}.")
  (:method ((x vector) i (j (eql 0)))
    (vref x i))
  (:method ((x <vector>) i (j (eql 0)))
    (vref x i)))

(defgeneric (setf mref) (value A i j)
  (:documentation "Writes the matrix element @code{A[i,j]}.")
  (:method (value (x vector) i (j (eql 0)))
    (setf (vref x i) value))
  (:method (value (x <vector>) i (j (eql 0)))
    (setf (vref x i) value)))

;;; For these routines, inlining may help a lot.

(defgeneric for-each-row-key (func mat)
  (:documentation "Loop through row keys."))

(defgeneric for-each-col-key (func mat)
  (:documentation "Loop through column keys."))

(defgeneric for-each-key-in-row (func mat row-key)
  (:documentation "Loop through col-keys in row."))

(defgeneric for-each-key-in-col (func mat col-key)
  (:documentation "Loop through row-keys in column col."))

(defgeneric for-each-entry-in-row (func mat row-key)
  (:documentation "Loop through col-keys in row."))

(defgeneric for-each-entry-in-col (func mat col-key)
  (:documentation "Loop through entries in column col."))

(defgeneric for-each-key-and-entry-in-row (func mat row-key)
  (:documentation "Loop through col-keys and entries in row."))

(defgeneric for-each-key-and-entry-in-col (func mat col-key)
  (:documentation "Loop through row-keys and entries in col."))

(defmacro dorows ((key mat) &body body)
  "Syntax: @slisp{(dorows (key mat) ...)}."
  `(for-each-row-key (lambda (,key) ,@body) ,mat))

(defmacro docols ((key mat) &body body)
  "Syntax: @slisp{(docols (key mat) ...)}."
  `(for-each-col-key (lambda (,key) ,@body) ,mat))

(defmacro dorow ((loop-vars mat row) &body body)
  (multiple-value-bind (key entry func)
      (cond ((atom loop-vars) (values loop-vars nil 'for-each-key-in-row))
            ((single? loop-vars) (values nil (first loop-vars) 'for-each-entry-in-row))
            (t (values (first loop-vars) (second loop-vars) 'for-each-key-and-entry-in-row)))
    `(,func (lambda (,@(when key (list key))
                     ,@(when entry (list entry)))
              ,@body) ,mat ,row)))

(defmacro docol ((loop-vars mat col) &body body)
  (multiple-value-bind (key entry func)
      (cond ((atom loop-vars) (values loop-vars nil 'for-each-key-in-col))
            ((single? loop-vars) (values nil (first loop-vars) 'for-each-entry-in-col))
            (t (values (first loop-vars) (second loop-vars) 'for-each-key-and-entry-in-col)))
    `(,func (lambda (,@(when key (list key))
                     ,@(when entry (list entry)))
              ,@body) ,mat ,col)))

;;;; derived functionality

(defgeneric row-keys (mat)
  (:documentation "All row keys for a matrix.")
  (:method ((mat <matrix>))
    (mapper-collect #'for-each-row-key mat)))

(defgeneric col-keys (mat)
  (:documentation "All column keys for a matrix.")
  (:method ((mat <matrix>))
    (mapper-collect #'for-each-col-key mat)))

(defgeneric keys-of-row (mat key)
  (:documentation "All column keys in the given row for a matrix.")
  (:method ((mat <matrix>) key)
    (mapper-collect #'for-each-key-in-row mat key)))

(defgeneric keys-of-column (mat key)
  (:documentation "All row keys in the given column for a matrix.")
  (:method ((mat <matrix>) key)
    (mapper-collect #'for-each-key-in-col mat key)))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; specializations of vector methods
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmethod vref ((mat <matrix>) (indices list))
  "Vector referencing on matrices is done by default by matrix referencing
a list of indices."
  (apply #'mref mat indices))

(defmethod (setf vref) (value (mat <matrix>) (indices list))
  "Vector referencing on matrices is done by default by matrix referencing
a list of indices."
  (setf (apply #'mref mat indices) value))

(defmethod for-each-entry-and-key (fn (mat <matrix>))
  (for-each-row-key
   (lambda (rk)
     (for-each-key-and-entry-in-row
      (lambda (ck entry)
        (funcall fn entry rk ck))
      mat rk))
   mat))

(defmethod for-each-entry-and-vector-index (func (mat <matrix>))
  (for-each-entry-and-key
   (lambda (entry i j) (funcall func entry (list i j)))
   mat))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; matrix tests
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defgeneric msquare-p (mat)
  (:documentation "Returns T, iff @arg{mat} is square.")
  (:method (mat)
    (= (nrows mat) (ncols mat))))

(defgeneric msymmetric-p (mat)
  (:documentation "Returns T, if @arg{mat} is symmetric.")
  (:method (mat)
    (mequalp mat (transpose mat))))

(defgeneric midentity-p (number &optional threshold)
  (:documentation "Returns T, if @arg{mat} is the identity, i.e. if the
  elementwise difference to the identity is not larger than
  @arg{threshold}."))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; matrix printing
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defgeneric show (matrix &key &allow-other-keys)
  (:documentation "Shows the contents of @arg{matrix} in a readable form.")
  (:method :around (matrix &key &allow-other-keys)
           (call-next-method)
           matrix)
  (:method (matrix &rest args)
    "The default method describes its argument."
    (declare (ignore args))
    (describe matrix)))

(defgeneric display (matrix &rest args)
  (:documentation "Formats the contents of @arg{matrix} in rectangular
  form.")
  (:method (matrix &rest args)
    "The default method describes its argument."
    (declare (ignore args))
    (describe matrix)))
  

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; special matrix methods
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; BLAS Level 3
(defgeneric gemm-nn! (a x y b z)
  (:documentation "General matrix-matrix multiplication:
@math{Z <- alpha * X * Y + beta * Z}"))
(defgeneric gemm-nt! (a x y b z)
  (:documentation "General matrix-matrix multiplication:
@math{Z <- alpha * X * Y' + beta * Z}"))
(defgeneric gemm-tn! (a x y b z)
  (:documentation "General matrix-matrix multiplication:
@math{Z <- alpha * X' * Y + beta * Z}"))
(defgeneric gemm-tt! (a x y b z)
  (:documentation "General matrix-matrix multiplication:
@math{Z <- alpha * X' * Y' + beta * Z}"))

(defmacro define-default-gemm! (name At yt)
  (multiple-value-bind (i j)
      (if At (values 'j 'i) (values 'i 'j))
    `(defmethod ,name (alpha A y beta z)
      (unless (= beta 1) (scal! beta z))
      (dovec ((entry i j) A z)
	(dotimes (l (multiplicity z))
	  (setf (mref z ,i l)
		(,name alpha entry (mref y ,@(if yt `(l ,j) `(,j l)))
		       1.0 (mref z ,i l))))))))

(define-default-gemm! gemm-nn! nil nil)
(define-default-gemm! gemm-nt! nil t)
(define-default-gemm! gemm-tn! t nil)
(define-default-gemm! gemm-tt! t t)

;;; LAPACK
(defgeneric getrf! (A &optional ipiv)
  (:documentation "Computes the PA=LU decomposition of @arg{A} which is
stored again in @arg{A}.  @arg{ipiv} can be a pre-allocated vector which
the routine fills with the indices for column pivoting, or NIL which
implies that the routine allocates such a vector itself.  If @arg{ipiv} is
@symbol{:none}, no pivoting is done.  Returns @arg{A} as the first value,
the pivot vector as a second value, and a boolean as the third value
indicating that the decomposition succeeded."))

(defgeneric getrs! (LU b &optional ipiv)
  (:documentation "Solves the PA=LU decomposition specified by @arg{LU} and
@arg{ipiv} for the rhs @arg{b}.  The result is stored in @arg{b}."))

;;; Special
(defgeneric transpose! (x y)
  (:documentation "Sets Y to the transpose of X."))
(defgeneric matrix-transpose-instance (x)
  (:documentation "Returns a zero matrix for storing the transpose of X."))

;;; Matrix manipulation

(defgeneric minject! (x y row-offset col-offset)
  (:documentation "Inject matrix X in matrix Y at the position given by
ROW-OFFSET and COL-OFFSET."))

(defgeneric mextract! (x y row-offset col-offset)
  (:documentation "Extract matrix X out of matrix Y from the position given
by ROW-OFFSET and COL-OFFSET."))

;;; derived functionality
(defgeneric vector-slice (x offset size)
  (:documentation "Extract a subvector of size @arg{size} out of @arg{x}
starting from position @arg{offset}."))

(defgeneric matrix-slice (x &key from-row from-col nrows ncols)
  (:documentation "Extract a submatrix of size @arg{nrows} @math{times}
@arg{ncols} out of @arg{x} starting from position
@arg{from-row}/@arg{from-col}."))

(defgeneric join-horizontal! (result &rest matrices)
  (:documentation "Joins @arg{matrices} horizontally into result.")
  (:method (result &rest matrices)
    (loop with n = (nrows result) and k = 0
       for mat in matrices do
	 (assert (= n (nrows mat)))
	 (minject! mat result 0 k)
	 (incf k (ncols mat))
       finally (assert (= k (ncols result))))
    result))

(defgeneric join-vertical! (result &rest matrices)
  (:documentation "Joins @arg{matrices} vertically into result.")
  (:method (result &rest matrices)
    (loop with n = (ncols result) and k = 0
       for mat in matrices do
	 (assert (= n (ncols mat)))
	 (minject! mat result k 0)
	 (incf k (nrows mat))
       finally (assert (= k (nrows result))))
    result))

(defgeneric join-instance (orientation matrix &rest matrices)
  (:documentation "Compute an instance for storing the join of
  @arg{orientation} applied to matrix and matrices."))

(defun join (orientation &rest matrices)
  "Joins @arg{matrices} either horizontally or vertically depending on
@arg{orientation}.  Due to the call to @function{zeros} this is not yet a
generic function."
  (unless matrices (error "No arguments to join"))
  (apply (ecase orientation
	   (:horizontal #'join-horizontal!)
	   (:vertical #'join-vertical!))
	 (apply #'join-instance orientation matrices)
	 matrices))

;;; Matrix-vector routines

(defgeneric make-domain-vector-for (mat &optional multiplicity))
(defgeneric make-image-vector-for (mat &optional multiplicity))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; General methods
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmethod nrows ((mat <matrix>))
  "Default and very inefficient method for computing the number of rows of a
matrix."
  (lret ((n 0))
    (for-each-row-key (_ (incf n)) mat)))

(defmethod ncols ((mat <matrix>))
  "Default and very inefficient method for computing the number of columns of a
matrix."
  (lret ((n 0))
    (for-each-col-key (_ (incf n)) mat)))

(defmethod multiplicity ((vec <matrix>))
  "If @arg{vec} should be a matrix, the multiplicity is the number of
columns by default."
  (ncols vec))
  
(defmethod matrix-transpose-instance (x)
  (make-instance (class-of x) :nrows (ncols x) :ncols (nrows x)))

(defgeneric transpose (x)
  (:documentation "Transpose the matrix @arg{x}.")
  (:method (x)
    "The default method casts performs a combination of
@function{transpose!} and @function{matrix-transpose-instance}."
    (transpose! x (matrix-transpose-instance x))))

(defun gemm! (alpha x y beta z &optional (job :nn))
  "Dispatches on the optional job argument (member :nn :tn :nt :tt) and
calls the corresponding generic function, e.g. GEMM-NN!."
  (ecase job
    (:nn (gemm-nn! alpha x y beta z))
    (:nt (gemm-nt! alpha x y beta z))
    (:tn (gemm-tn! alpha x y beta z))
    (:tt (gemm-tt! alpha x y beta z))))

(defun gemm (alpha x y beta z  &optional (job :nn))
  "Rewriting of GEMM in terms of GEMM!."
  (gemm! alpha x y beta (copy z) job))

(defgeneric m*-product-instance (x y)
  (:documentation "Allocates an instance for the product of X and Y.")
  (:method (x y)
    (make-instance (class-of y) :nrows (nrows x) :ncols (ncols y))))

(defgeneric m* (x y)
  (:documentation "Multiply X by Y.")
  (:method (x y)
    "By default M* is rewritten in terms of GEMM!."
    (gemm-nn! (coerce 1 (scalar-type x)) x y
	      (coerce 0 (scalar-type x)) (m*-product-instance x y))))
  
(defgeneric m*-tn-product-instance (x y)
  (:documentation "Allocates an instance for the product of X^t and Y.")
  (:method (x y)
    (make-instance (class-of y) :nrows (ncols x) :ncols (ncols y))))

(defgeneric m*-tn (x y)
  (:documentation "Multiply X^t by Y.")
  (:method (x y)
    "By default, M*-TN is rewritten in terms of GEMM!."
    (gemm-tn! (coerce 1 (scalar-type x)) x y
	      (coerce 0 (scalar-type x)) (m*-tn-product-instance x y))))

(defun getrf (x &optional ipiv)
  "Rewriting for GETRF  in terms of GETRF!."
  (getrf! (copy x) ipiv))

(defun getrs (lu b &optional ipiv)
  "Rewriting for GETRS in terms of GETRS!."
  (getrs! lu (copy b) ipiv))

(defmethod gesv! (A b)
  "Solves a linear system A X = B for X via GETRF and GETRS!."
  (multiple-value-bind (LR ipiv info)
      (getrf A)
    (if (numberp info)
	(error "argument A given to GESV! is singular to working machine precision")
	(getrs! LR b ipiv))))

(defun gesv (A b)
  "Rewriting for GESV in terms of GESV!."
  (gesv! A (copy b)))

(defmethod mat-diff (mat1 mat2)
  (dovec ((entry i j) mat1)
    (whereas ((entry2 (mref mat2 i j)))
      (unless (mzerop (m- entry entry2))
	(format t "(~A,~A) : ~A <--> ~A~%" i j entry entry2)))))
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Submatrices
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass <submatrix> (<matrix>)
  ((matrix :initarg :matrix :documentation "The \"supermatrix\".")
   (row-keys :reader row-keys :initarg :row-keys :type vector
	     :documentation "The row indices of the submatrix.")
   (col-keys :reader col-keys :initarg :col-keys :type vector
	     :documentation "The column indices of the submatrix."))
  (:documentation "Describes an ordered submatrix of a matrix.  Only a
restricted set of operations is allowed for these matrices and element
access is slow.  They are indexed with ordinary integers."))

(defmethod mref ((submat <submatrix>) i j)
  (mref (slot-value submat 'matrix)
	(aref (row-keys submat) i)
	(aref (col-keys submat) j)))

(defmethod (setf mref) (value (submat <submatrix>) i j)
  (setf (mref (slot-value submat 'matrix)
	      (aref (row-keys submat) i)
	      (aref (col-keys submat) j))
	value))

(defmethod vref ((submat <submatrix>) k)
  "A @class{<submatrix>} can be accessed as one-dimensional vector."
  (multiple-value-bind (i j) (floor k (nrows submat))
    (mref submat i j)))

(defmethod (setf vref) (value (submat <submatrix>) k)
  (multiple-value-bind (i j) (floor k (nrows submat))
    (setf (mref submat i j) value)))

(defmethod for-each-entry-and-key (func (mat <submatrix>))
  (dotimes (i (nrows mat))
    (dotimes (j (ncols mat))
      (funcall func (mref mat i j) i j))))

(defmethod nrows ((mat <submatrix>)) (length (row-keys mat)))

(defmethod ncols ((mat <submatrix>)) (length (col-keys mat)))


