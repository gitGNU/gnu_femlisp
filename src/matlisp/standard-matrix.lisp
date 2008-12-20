;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; standard-matrix.lisp - Matlisp alternative
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

;;;; This is the implementation of the standard-matrix class and its
;;;; methods which replaces Matlisp within Femlisp.

(defclass standard-matrix (<matrix>)
  ((nrows :reader nrows :initarg :nrows
	  :documentation "Number of rows in the matrix")
   (ncols :reader ncols :initarg :ncols
	  :documentation "Number of columns in the matrix"))
  (:documentation "Mixin for dense matrices."))

;;; if standard-matrix is considered as multiple vectors
(defmethod vlength ((mat standard-matrix))
  (nrows mat))
(defmethod multiplicity ((mat standard-matrix))
  (ncols mat))
(defmethod total-entries ((mat standard-matrix))
   (* (nrows mat) (ncols mat)))

(defmethod element-type (matrix)
  "Default method returns T."
  (declare (ignore matrix))
  T)

(defmethod scalar-type (matrix)
  "Default method returns NUMBER."
  (declare (ignore matrix))
  'NUMBER)

(declaim (inline indexing))
(defun indexing (i j nrows ncols)
  "Computes column-major indexing for compatibility with Matlisp/Fortran."
  (declare (ignore ncols))
  (declare (type positive-fixnum i j nrows ncols))
  (the fixnum (+ i (the fixnum (* nrows j)))))

(defmethod initialize-instance ((matrix standard-matrix) &key content &allow-other-keys)
  "Handles the content parameter and sets the store to a suitable vector.
If content is a 2d array, the dimensions can be deduced."
  (call-next-method)
  (assert (typep matrix 'store-vector))
  (typecase content
    (sequence
     ;; put content in a normal form
     (setq content
	   (map 'list #'(lambda (row)
			  (etypecase row
			    (number (list row))
			    (list row)
			    (vector (coerce row 'list))))
		content))
     (unless (slot-boundp matrix 'nrows)
       (setf (slot-value matrix 'nrows) (length content)))
     (assert (same-p (map 'list #'length content)))
     (unless (slot-boundp matrix 'ncols)
       (setf (slot-value matrix 'ncols) (length (car content)))))
    (array
     (unless (slot-boundp matrix 'nrows)
       (setf (slot-value matrix 'nrows) (array-dimension content 0)))
     (unless (slot-boundp matrix 'ncols)
       (setf (slot-value matrix 'ncols) (array-dimension content 1)))))
  (with-slots (nrows ncols store)
      matrix
    (let ((n (* nrows ncols))
	  (type (element-type matrix)))
      (if (slot-boundp matrix 'store)
	  (assert (and (= n (length store))
		       (subtypep (upgraded-array-element-type type)
				 (array-element-type store))))
	  (setq store (if (subtypep type 'number)
			  (zero-vector n type)
			  (make-array n :element-type type))))
      (when content
	(etypecase content
	  (sequence
	   (loop for i from 0 and row in content do
		 (loop for j from 0 and entry in row do
		       (setf (aref store (indexing i j nrows ncols)) entry))))
	  (array
	   (assert (equal (array-dimensions content) (list nrows ncols)))
	   (dotimes (i nrows)
	     (dotimes (j ncols)
	       (setf (aref store (indexing i j nrows ncols)) (aref content i j))))))))))

(with-memoization (:type :local :size 2 :id 'standard-matrix)
  (defun standard-matrix (type)
    "Defines the programmatic class @class{standard-matrix} for element type
@arg{type} as extensions of the programmatic class @class{store-vector}."
    (memoizing-let ((type type))
      (fl.amop:find-programmatic-class
       (list 'standard-matrix (store-vector type))
       (intern (format nil "~A" (list 'STANDARD-MATRIX type)) "FL.MATLISP")))))

(inlining
 (defun standard-matrix-p (obj)
   "Tests if @arg{obj} is a @class{standard-matrix}."
   (typep obj 'standard-matrix)))

(defun make-real-matrix (&rest args)
  "Generates a real matrix as specified by its arguments.  If two arguments
are provided, they should be numbers which are interpreted as rows and
columns.  If only one argument is provided, it should be either a number
meaning the rows and columns of a square matrix or a nested list or vector
structure defining the contents matrix."
  (quickly
    (cond
      ((single? args)
       (let ((arg (car args)))
	 (cond ((numberp arg) (zeros arg arg))
	       (t (make-instance (standard-matrix 'double-float)
				 :content arg)))))
      ((single? (cdr args))
       (zeros (first args) (second args)))
      (t (error "Too many arguments.")))))

(eval-when (:compile-toplevel :load-toplevel :execute)
  (set-dispatch-macro-character
   #\# #\m  ; dispatch on #m for real matrices
   #'(lambda (stream char n)
       (declare (ignore char n))
       (let ((list (read stream nil nil t)))
	 ;; TODO: try without quoting
	 `(make-real-matrix ',list)))))

(defun make-real-vector (dim &optional (value 0.0))
  "Generates a real matrix of dimension @arg{dim} x 1."
  (make-instance (standard-matrix 'double-float)
		 :nrows dim :ncols 1
		 :content (make-double-vec dim value)))
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; conversion of vectors to standard-matrix
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun uniform-number-type (vec)
  "Tries to find a uniform type for the numbers contained in @arg{vec}."
  (lret ((element-type (array-element-type vec)))
    (unless (subtypep element-type 'number)
      (when (plusp (length vec))
	(setq element-type nil)
	(loop for x across vec
	      unless (subtypep (type-of x) element-type)
	      do (setq element-type
		       (upgraded-array-element-type
			`(or ,element-type ,(type-of x)))))))))

(defgeneric ensure-matlisp (obj &optional type)
  (:documentation "Tries to coerce @arg{obj} into Matlisp format.")
  (:method ((obj standard-matrix) &optional type)
    (declare (ignore type))
    obj)
  (:method ((obj number) &optional type)
    (declare (ignore type))
    (ensure-matlisp (vector obj)))
  (:method ((obj vector) &optional (type :column))
    (let ((element-type (uniform-number-type obj)))
      (unless (subtypep (array-element-type obj) element-type)
	(setq obj (coerce obj `(simple-array ,element-type (*)))))
      #+(or ecl gcl)
      (progn
	(when (eq element-type 'long-float)
	  (setq element-type 'double-float))
	(when (eq element-type 'short-float)
	  (setq element-type 'single-float)))
      (?2
       (make-instance
	(standard-matrix element-type)
	:store obj 
	:nrows (if (eq type :column) (length obj) 1)
	:ncols (if (eq type :row) (length obj) 1))
       ;; this branch avoids calling make-instance with keyword arguments
       ;; on a class unknown at compile time, because this is a
       ;; performance problem for several CL implementations
       (lret ((result (make-instance (standard-matrix element-type))))
	 (with-slots (store nrows ncols) result
	   (setf store obj
		 nrows (if (eq type :column) (length obj) 1)
		 ncols (if (eq type :row) (length obj) 1)))))))
  (:method ((obj list) &optional (type :column))
    (ensure-matlisp (coerce obj 'vector) type))
    )

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Access to the entries
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(declaim (inline standard-matrix-indexing))
(defun standard-matrix-indexing (i j nrows)
  (declare (type (and fixnum (integer 0)) i j))
  (the fixnum (+ i (the fixnum (* j nrows)))))

(defmethod mref ((matrix standard-matrix) i j)
  "We choose Fortran-like column-major indexing to make the use of Matlisp
routines easier.  Note: access to matrix entries using this generic
function is slow.  Therefore, specialized BLAS routines should be used
whenever possible."
  (declare (type (and fixnum (integer 0)) i j))
  (with-slots (store nrows ncols) matrix
    (assert (and (<= 0 i) (<= 0 j) (< i nrows) (< j ncols)))
    (aref store (standard-matrix-indexing i j nrows))))


(defmethod (setf mref) (value (matrix standard-matrix) i j)
  "Writer for mref."
  (declare (type (and fixnum (integer 0)) i j))
  (with-slots (store nrows ncols) matrix
    (assert (and (<= 0 i) (<= 0 j) (< i nrows) (< j ncols)))
    (setf (aref store (standard-matrix-indexing i j nrows))
	  value)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Printing matrices
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defvar *print-matrix* 5
  "Maximum number of columns and/or rows to print. NIL: no elements, T: all
elements.")

(defvar *print-matrix-pretty* nil
  "T means that a newline is printed after each row.")

(defvar *print-matrix-element-format* nil
  "Format of matrix element field to be printed.  A useful format is
  \"~10,2,2E\" for debugging purposes.")

(defgeneric print-element (matrix element stream)
  (:documentation "Prints the element @arg{element} of @arg{matrix} to the
stream @arg{stream}."))

(defmethod print-element ((matrix standard-matrix) element stream)
  (if *print-matrix-element-format*
      (format stream *print-matrix-element-format* element)
      (format stream "~S" element)))

(defun print-matrix (matrix stream)
  (if *print-matrix*
      (let ((*print-circle* nil))
	(if (and (slot-boundp matrix 'nrows)
		 (slot-boundp matrix 'ncols)
		 (slot-boundp matrix 'store))
	    (with-slots (nrows ncols) matrix
	      (let ((m nrows) (n ncols))
		(when (numberp *print-matrix*)
		  (setq m (min m *print-matrix*))
		  (setq n (min n *print-matrix*)))
		(princ "#M(" stream)
		(dotimes (i m)  
		  (princ "(" stream)
		  (dotimes (j n)
		    (print-element matrix (mref matrix i j) stream)
		    (unless (= j (1- n)) (format stream " ")))
		  (if (> ncols n)
		      (princ " ...)" stream)
		      (princ ")" stream))
		  (when *print-matrix-pretty*
		    (terpri stream))
		  (unless (= i (1- m)) (format stream " ")))
		(if (> nrows m)
		    (princ " ...)" stream)
		    (princ ")" stream))))
	    (princ "#M(??)" stream)))
      (print-unreadable-object (matrix stream :type t :identity t))))

(defmethod print-object ((matrix standard-matrix) stream)
  (print-matrix matrix stream))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Copy
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Special matrices
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; Constructors for special standard matrices

(defun zeros (n &optional (m n) (type 'double-float))
  "Returns nxn or (if m is provided) nxm zeros.  The value is freshly
allocated."
  (?2 (make-instance (standard-matrix type) :nrows n :ncols m)
      ;; this branch avoids calling make-instance with keyword arguments on
      ;; a class unknown at compile time, because this is a performance
      ;; problem for several CL implementations
      (lret ((result (make-instance (standard-matrix type))))
	(with-slots (store nrows ncols) result
	  (setf nrows n ncols m
		store (make-array (* n m) :element-type type
				  :initial-element (coerce 0 type)))))))

(defun eye (n &optional (m n) (type 'double-float))
  "Returns the nxn identity matrix.  The value is freshly allocated."
  (let ((result (zeros n m type))
	(one (coerce 1 type)))
    (dotimes (i (min n m) result)
      (setf (mref result i i) one))))

(defun ones (n &optional (m n) (type 'double-float))
  "Returns nxn or (if m is provided) nxm ones.  The value is freshly
allocated."
  (let ((vec (make-array (* n m) :element-type type
			 :initial-element (coerce 1 type))))
    (?2 (make-instance (standard-matrix type) :nrows n :ncols m :store vec)
	;; this branch avoids calling make-instance with keyword arguments
	;; on a class unknown at compile time, because this is a
	;; performance problem for several CL implementations
	(lret ((result (make-instance (standard-matrix type))))
	  (with-slots (store nrows ncols) result
	    (setf nrows n ncols m store vec))))))

(defun column (&rest values)
  (lret ((result (zeros (length values) 1)))
    (loop for value in values and i from 0 do
	  (setf (vref result i) value))))

(defun row (&rest values)
  (lret ((result (zeros 1 (length values))))
    (loop for value in values and i from 0 do
	  (setf (vref result i) value))))

(defun mrandom (n &optional m (type 'double-float) (range 1.0))
  "Returns a random nxn or (if m is provided) nxm matrix.  The value is
freshly allocated."
  (unless m (setq m n))
  (fill-random! (zeros n m) (coerce range type)))

(defun diag (vec)
  "Returns a diagonal matrix with diagonal entries from vec."
  (let* ((n (vlength vec))
	 (result (zeros n n (element-type vec))))
    (dotimes (i n result)
      (setf (mref result i i) (vref vec i)))))

(defmethod diagonal ((A standard-matrix))
  "Returns the diagonal of @arg{A} as a vector."
  (lret* ((n (min (nrows A) (ncols A)))
	  (result (zero-vector n (element-type A))))
    (dotimes (i n)
      (setf (aref result i) (mref A i i)))))

(defun stencil-matrix (stencil sizes)
  "@arg{stencil} should be an array of rank dim, each dimension being 3,
@arg{sizes} is either an integer or an integer vector of length dim.  Each
size denotes the number of interior mesh-points in the respective
dimension."
  (declare (optimize debug))
  (let ((dim (array-rank stencil)))
    (when (numberp sizes)
      (setq sizes (make-list dim :initial-element sizes)))
    (let* ((n (reduce #'* sizes))
	   (offsets (loop for k from 0
		       and size in sizes
		       for offset = 1 then (* offset size)
		       collect offset))
	   (result (zeros n)))
      (flet ((ind->off (indices)
	       (loop for k in indices and off in offsets
		  summing (* k off))))
	(dotuple (i sizes)
	  (let ((ind_i (ind->off i)))
	    (dotuple (j (make-list dim :initial-element 3))
	      (let ((ij (mapcar (lambda (x y) (+ x y -1)) i j)))
		(when (every (lambda (x s) (< -1 x s)) ij sizes)
		  (let ((ind_ij (ind->off ij)))
		    (setf (mref result ind_i ind_ij)
			  (apply #'aref stencil j)))))))))
      result)))

(defun laplace-stencil (dim)
  "Returns the @{2*dim+1}-point stencil for the Laplace matrix on a uniform
mesh as an array of rank @arg{dim}."
  (lret ((result (make-array (make-list dim :initial-element 3)
			     :initial-element 0.0)))
    (setf (apply #'aref result (make-list dim :initial-element 1))
	  (* 2.0 dim))
    (flet ((ind (k dir)
	     (append (make-list k :initial-element 1)
		     (list (ecase dir (+ 2) (- 0)))
		     (make-list (- dim k 1) :initial-element 1))))
      (loop for k below dim do
	   (setf (apply #'aref result (ind k '+)) -1.0
		 (apply #'aref result (ind k '-)) -1.0)))))

(defun laplace-full-matrix (n &optional (dim 1))
  "Generates the matrix for a @arg{dim}-dimensional Laplace problem
discretized with the @math{2*@arg{dim}+1}-point stencil on a structured
mesh with Dirichlet boundary conditions."
  (stencil-matrix (laplace-stencil dim)
		  (make-list dim :initial-element n)))

(defmethod make-analog ((x standard-matrix))
  "Constructs a zero matrix with the same size as X."
  (zeros (nrows x) (ncols x) (element-type x)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Join of standard matrices
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmethod join-instance (orientation (matrix standard-matrix) &rest matrices)
  (let (etype m n)
    (loop for mat in (cons matrix matrices) do
	 (cond ((null etype)
		(setq  etype (element-type mat)
		       m (nrows mat)
		       n (ncols mat)))
	       (t
		(assert (equalp etype (element-type mat)))
		(ecase orientation
		  (:horizontal (assert (= m (nrows mat)))
			       (incf n (ncols mat)))
		  (:vertical (assert (= n (ncols mat)))
			     (incf m (nrows mat)))))))
    (zeros m n etype)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; LAPACK routines
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmethod fl.lapack::lapack-convert ((mat standard-matrix))
  (fl.port::vector-sap (store mat)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Iteration
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; this old interface will probably soon be replaced by BLAS macros.

(defmethod for-each-entry (fn (mat standard-matrix))
  (dotimes (i (nrows mat))
    (dotimes (j (ncols mat))
      (funcall fn (mref mat i j)))))
(defmethod for-each-row-key (fn (mat standard-matrix))
  "Loop through row keys."
  (dotimes (i (nrows mat))
    (funcall fn i)))
(defmethod for-each-key-in-row (fn (mat standard-matrix) row-key)
  "Loop through column keys in row."
  (declare (ignore row-key))
  (dotimes (i (ncols mat))
    (funcall fn i)))
(defmethod for-each-col-key (fn (mat standard-matrix))
  "Loop through column keys."
  (dotimes (i (ncols mat))
    (funcall fn i)))
(defmethod for-each-entry-and-key (fn (mat standard-matrix))
  (dotimes (i (nrows mat))
    (dotimes (j (ncols mat))
      (funcall fn (mref mat i j) i j))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Testing
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(eval-when (:load-toplevel :execute)
  ;; to have the #M reader macro working in its modified form this use is
  ;; wrapped in an EVAL-WHEN
  (defun test-standard-matrix ()
    ;;(setq *print-matrix* nil)
    (join-instance :horizontal #m(1.0) #m(1.0))
    (make-instance (standard-matrix 'integer) :content #2a((2 3) (4 5)))
    (describe (standard-matrix 'double-float))
    (standard-matrix 'single-float)
    (describe (eye 1))
    (eye 2)
    (zeros 1)
    (row (cos 1.0))
    (make-real-matrix #((2.0 3.0) (4.0 5.0)))
    (copy (eye 5))
    (norm (row 1.0 2.0) 1)
    )
  )

;;; (fl.matlisp::test-standard-matrix)
(fl.tests:adjoin-test 'test-standard-matrix)