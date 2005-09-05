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
  (:documentation "Number of matrix rows."))

(defgeneric mref (A i j)
  (:documentation "Returns the matrix element @code{A[i,j]}."))
(defgeneric (setf mref) (value A i j)
  (:documentation "Writes the matrix element @code{A[i,j]}."))

(defmethod vref ((mat <matrix>) (index-pair cons))
  "Vector referencing on matrices is done by default by matrix referencing
a pair of index-pair.  Some matrices may allow for other vector indexing
schemes."
  (mref mat (car index-pair) (cdr index-pair)))

(defmethod (setf vref) (value (mat <matrix>) (index-pair cons))
  "Vector referencing on matrices is done by default by matrix referencing
a pair of index-pair."
  (setf (mref mat (car index-pair) (cdr index-pair)) value))


;;; BLAS Level 3
(defgeneric gemm-nn! (a x y b z)
  (:documentation "General matrix-matrix multiplication: Z <- alpha * X * Y
+ beta * Z"))
(defgeneric gemm-nt! (a x y b z)
  (:documentation "General matrix-matrix multiplication: Z <- alpha * X * Y'
+ beta * Z"))
(defgeneric gemm-tn! (a x y b z)
  (:documentation "General matrix-matrix multiplication: Z <- alpha * X' * Y
+ beta * Z"))
(defgeneric gemm-tt! (a x y b z)
  (:documentation "General matrix-matrix multiplication: Z <- alpha * X' * Y'
+ beta * Z"))

(defgeneric m*-product-instance (x y)
  (:documentation "Returns a zero matrix for storing the product of X and Y."))

;;; LAPACK
(defgeneric getrf! (x &optional ipiv)
  (:documentation "Computes the PLU decomposition of X (overwriting X).
Returns X and as a second value the permuations vector."))
(defgeneric getrs! (x b &optional ipiv)
  (:documentation "Solves the given PLU decomposition for the rhs b while
overwriting b."))

;;; Special
(defgeneric transpose! (x y)
  (:documentation "Sets Y to the transpose of X."))
(defgeneric matrix-transpose-instance (x)
  (:documentation "Returns a zero matrix for storing the transpose of X."))

;;; Matrix manipulation

(defgeneric join (x y &optional orientation)
  (:documentation "Joins X and Y horizontally or vertically depending on the
value of orientation."))

(defgeneric matrix-slice (matrix &key from-row from-col nrows ncols)
  (:documentation "Return the specified submatrix of MATRIX."))

(defgeneric vector-slice (vector offset size)
  (:documentation "Returns the specified subvector."))

(defgeneric minject (x y row-offset col-offset)
  (:documentation "Inject matrix X in matrix Y at the position given by
ROW-OFFSET and COL-OFFSET."))

(defgeneric mextract (x y row-offset col-offset)
  (:documentation "Extract matrix Y out of matrix X from the position given
by ROW-OFFSET and COL-OFFSET."))

(defgeneric vector-slice (x offset size)
  (:documentation "Extract a subvector of size @arg{size} out of @arg{x}
starting from position @arg{offset}."))

(defgeneric matrix-slice (x &key from-row from-col nrows ncols)
  (:documentation "Extract a submatrix of size @arg{nrows} @math{times}
@arg{ncols} out of @arg{x} starting from position
@arg{from-row}/@arg{from-col}."))

;;; Matrix-vector routines

(defgeneric make-domain-vector-for (mat &optional multiplicity))
(defgeneric make-image-vector-for (mat &optional multiplicity))

;;; Old iteration interface.  For these routines, inlining may help a lot.

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

(defmacro dorows ((loop-vars mat) &body body)
  (let* ((loop-vars (if (consp loop-vars) loop-vars (list loop-vars)))
	 (mat-for-each-row
	  (if (null (car loop-vars))
	      (if (null (cdr loop-vars))
		  (error "no loop variable")
		  'for-each-row)
	      (if (or (null (cdr loop-vars)) (null (cadr loop-vars)))
		  'for-each-row-key
		  'for-each-key-and-row))))
    `(,mat-for-each-row #'(lambda ,(remove nil loop-vars) ,@body) ,mat)))

(defmacro dorow ((loop-vars row) &body body)
  (let* ((loop-vars (if (consp loop-vars) loop-vars (list nil loop-vars)))
	 (row-for-each (if (null (car loop-vars))
			   (if (null (cdr loop-vars))
			       (error "no loop variable")
			       'for-each-entry-in-row)
			   (if (or (null (cdr loop-vars)) (null (cadr loop-vars)))
			       'for-each-key-in-row
			       'for-each-key-and-entry-in-row))))
    `(,row-for-each #'(lambda ,(remove nil loop-vars) ,@body) ,row)))

;;; vector-matrix operations
(defgeneric x<-Ay (x A y))
(defgeneric x+=Ay (x A y))
(defgeneric x-=Ay (x A y))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; General methods
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmethod multiplicity (vec)
  "If @arg{vec} is a matrix, the multiplicity is the number of columns by
default."
  (ncols vec))
  
(defmethod matrix-transpose-instance (x)
  (make-instance (class-of x) :nrows (ncols x) :ncols (nrows x)))

(defun transpose (x)
  "Rewriting for TRANSPOSE."
  (transpose! x (matrix-transpose-instance x)))

(defun gemm! (alpha x y beta z &optional (job :nn))
  "Dispatches on the optional job argument (member :nn :tn :nt :tt) and
calls the corresponding generic function, e.g. GEMM-NN!."
  (ecase job
    (:nn (gemm-nn! alpha x y beta z))
    (:nt (gemm-nt! alpha x y beta z))
    (:tn (gemm-tn! alpha x y beta z))
    (:tt (gemm-tt! alpha x y beta z))))

(definline gemm (alpha x y beta z  &optional (job :nn))
  "Rewriting of GEMM in terms of GEMM!."
  (gemm! alpha x y beta (copy z) job))

(defmethod m*-product-instance (x y)
  "Default method allocates an instance of class of Y."
  (make-instance (class-of y) :nrows (nrows x) :ncols (ncols y)))

(defmethod m*-product-instance-tn (x y)
  "Default method allocates an instance of class of Y."
  (make-instance (class-of y) :nrows (ncols x) :ncols (ncols y)))

(defmethod m* (x y)
  "By default M* is rewritten in terms of GEMM!."
  (gemm-nn! (coerce 1 (scalar-type x)) x y
	    (coerce 0 (scalar-type x)) (m*-product-instance x y)))
  
(defun m*-tn (x y)
  "By default M*-TN is rewritten in terms of GEMM-TN!."
  (gemm-tn! (coerce 1 (scalar-type x)) x y
	    (coerce 0 (scalar-type x)) (m*-product-instance-tn x y)))
  
(definline getrf (x &optional ipiv)
  "Rewriting for GETRF  in terms of GETRF!."
  (getrf! (copy x) ipiv))

(definline getrs (lu b &optional ipiv)
  "Rewriting for GETRS in terms of GETRS!."
  (getrs! lu (copy b) ipiv))

(defmethod gesv! (A b)
  "Solves a linear system A X = B for X via GETRF and GETRS!."
  (multiple-value-bind (LR ipiv info)
      (getrf A)
    (if (numberp info)
	(error "argument A given to GESV! is singular to working machine precision")
	(getrs! LR b ipiv))))

(definline gesv (x b)
  "Rewriting for GESV in terms of GESV!."
  (gesv! x (copy b)))

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

(defmethod for-each-key ((fn function) (mat <submatrix>))
  (dotimes (i (* (nrows mat) (ncols mat)))
    (funcall fn i)))

(defmethod for-each-key-and-entry ((fn function) (mat <submatrix>))
  (dotimes (i (nrows mat))
    (dotimes (j (ncols mat))
      (funcall fn (+ (* j (nrows mat)) i) (mref mat i j)))))

(defmethod for-each-entry ((fn function) (mat <matrix>))
  (dotimes (i (nrows mat))
    (dotimes (j (ncols mat))
      (funcall fn (mref mat i j)))))

(defmethod nrows ((mat <submatrix>)) (length (row-keys mat)))

(defmethod ncols ((mat <submatrix>)) (length (col-keys mat)))

(defmethod nr-of-entries ((submat <submatrix>))
  (* (nrows submat) (ncols submat)))



