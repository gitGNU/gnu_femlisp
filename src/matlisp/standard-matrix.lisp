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
	  (setq store (zero-vector n type)))
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

(defun standard-matrix (type)
  "Defines the programmatic class @class{standard-matrix} for element type
@arg{type} as extensions of the programmatic class @class{store-vector}."
  (assert (subtypep type 'number))
  (fl.amop:find-programmatic-class
   (list 'standard-matrix (store-vector type))
   (intern (format nil "~A" (list 'STANDARD-MATRIX type)) "FL.MATLISP")))

(defun make-real-matrix (&rest args)
  "Generates a real matrix as specified by its arguments.  If two arguments
are provided, they should be numbers which are interpreted as rows and
columns.  If only one argument is provided, it should be either a number
meaning the rows and columns of a square matrix or a nested list or vector
structure defining the contents matrix."
  (let ((class (standard-matrix 'double-float)))
    (cond
      ((= (length args) 1)
       (let ((arg (car args)))
	 (cond ((typep arg class) (copy arg))
	       ((numberp arg) (make-instance class :nrows arg :ncols arg))
	       (t (make-instance class :content arg)))))
      ((= (length args) 2)
       (make-instance class :nrows (first args) :ncols (second args)))
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

(defun complex-to-real-vector (cvec)
  (declare (type (simple-array (complex double-float) *) cvec))
  (declare (optimize speed (safety 0)))
  (let* ((n (length cvec))
	 (n2 (* 2 n)))
    (let ((dvec (make-double-vec n2)))
      (dotimes (i n)
	(let ((c (aref cvec i))
	      (k (* 2 i)))
	  (setf (aref dvec k) (realpart c))
	  (setf (aref dvec (1+ k)) (imagpart c))))
      dvec)))

(defun real-to-complex-vector (dvec)
  (declare (type double-vec dvec))
  (declare (optimize speed (safety 0)))
  (let ((n2 (length dvec)))
    (assert (evenp n2))
    (let ((n  (/ n2 2)))
      (declare (type fixnum n))
      (let ((cvec (make-array n :element-type '(complex double-float))))
	(dotimes (i n)
	  (setf (aref cvec i)
		(let ((k (* 2 i)))
		  (complex (aref dvec k) (aref dvec (1+ k))))))
	cvec))))

(defun extend-matlisp-function (func)
  "If @package{MATLISP} is available, and the argument @arg{func} is a
generic function in this package, this function is extended to be
applicable to matrices in @arg{FL.MATLISP}.  This is done by defining a
method for @function{no-applicable-method} which converts the arguments,
calls @arg{func} and reconverts the returned values.  If @package{MATLISP}
is not available, NIL is returned."
  (whereas ((matlisp-package (find-package "MATLISP")))
    (when (and (typep func 'generic-function)
	       (eq (symbol-package (fl.amop::generic-function-name func))
		   matlisp-package))
      (destructuring-bind (standard-matrix real-matrix complex-matrix
					   nrows ncols store)
	  (mapcar (rcurry #'intern matlisp-package)
		  '("STANDARD-MATRIX" "REAL-MATRIX" "COMPLEX-MATRIX"
		    "NROWS" "NCOLS" "STORE"))
	(defmethod no-applicable-method ((gf (eql func)) &rest args)
	  (let (alist)
	    ;; setup translation table with arguments
	    (loop for obj in args do
		 (when (and (typep obj 'fl.matlisp:standard-matrix)
			    (not (assoc obj alist)))
		   (push
		    (cons obj
			  (if (eq (element-type obj) 'double-float)
			      (make-instance real-matrix
					     :nrows (nrows obj) :ncols (ncols obj) :store (store obj))
			      (make-instance complex-matrix
					     :nrows (nrows obj) :ncols (ncols obj) :store
					     (complex-to-real-vector (store obj)))))
		    alist)))
	    (unless alist
	      (error "No matching method for the generic function ~S, when called with arguments ~S."
		     gf args))
	    (let ((return-values
		   (multiple-value-list
		    (apply gf (loop for obj in args collect
				   (or (geta alist obj) obj))))))
	      ;; augment translation table with returned values
	      (loop for obj in return-values do
		   (when (and (typep obj standard-matrix)
			      (not (rassoc obj alist)))
		     (push (cons (make-instance (fl.matlisp:standard-matrix
						 (cond ((typep obj real-matrix) 'double-float)
						       ((typep obj complex-matrix) '(complex double-float))
						       (t (error "Unknown Matlisp matrix."))))
						:nrows (funcall nrows obj) :ncols (funcall ncols obj)
						:store (let ((store (funcall store obj)))
							 (if (typep obj real-matrix)
							     store
							     (real-to-complex-vector store))))
				 obj)
			   alist)))
	      (apply #'values (loop for obj in return-values collect
				   (or (car (rassoc obj alist)) obj))))))))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; conversion of vectors to standard-matrix
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun ensure-matlisp (vec &optional (type :column))
  (etypecase vec
    (vector
     (let ((element-type (array-element-type vec)))
       (assert (subtypep element-type 'number)) ; preliminary
       #+ecl
       (progn
	 (when (eq element-type 'long-float)
	   (setq element-type 'double-float))
	 (when (eq element-type 'short-float)
	   (setq element-type 'single-float)))
       (make-instance
	(standard-matrix element-type)
	:store vec
	:nrows (if (eq type :column) (length vec) 1)
	:ncols (if (eq type :row) (length vec) 1))))
    (number
     (let ((result (make-instance (standard-matrix (type-of vec)) :nrows 1 :ncols 1)))
       (setf (vref result 0) vec)
       result))
    (standard-matrix vec)
    ))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Access to the entries
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmethod mref ((matrix standard-matrix) i j)
  "We choose Fortran-like column-major indexing to make the use of Matlisp
routines easier.  Note: access to matrix entries using this generic
function is slow.  Therefore, specialized BLAS routines should be used
whenever possible."
  (declare (type (and fixnum (integer 0)) i j))
  (with-slots (store nrows ncols) matrix
    (assert (and (<= 0 i) (<= 0 j) (< i nrows) (< j ncols)))
    (aref store (the fixnum (+ i (the fixnum (* j nrows)))))))


(defmethod (setf mref) (value (matrix standard-matrix) i j)
  "Writer for mref."
  (declare (type (and fixnum (integer 0)) i j))
  (with-slots (store nrows ncols) matrix
    (assert (and (<= 0 i) (<= 0 j) (< i nrows) (< j ncols)))
    (setf (aref store (the fixnum (+ i (the fixnum (* j (nrows matrix))))))
	  value)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Printing
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defvar *print-matrix* 5
  "Maximum number of columns and/or rows to print. NIL: no elements, T: all
elements.")

(defgeneric print-element (matrix element stream)
  (:documentation "Prints a matrix element to stream."))

(defmethod print-element ((matrix standard-matrix) element stream)
  (format stream "~a" element))

(defun print-matrix (matrix stream)
  (if *print-matrix*
      (let ((*print-circle* nil))
	  (cond
	    ((and (slot-boundp matrix 'nrows)
		  (slot-boundp matrix 'ncols)
		  (slot-boundp matrix 'store))
	     (with-slots (nrows ncols store)
		 matrix
	       (loop with m = nrows and n = ncols
		  initially
		  (princ "#M(" stream)
		  (when (numberp *print-matrix*)
		    (setq m (min nrows *print-matrix*))
		    (setq n (min ncols *print-matrix*)))
		  for i below m do
		  (princ "(" stream)
		  (loop for j below n do
		       (print-element matrix (mref matrix i j) stream)
		       (unless (= j (1- n)) (format stream " "))
		       finally
		       (when (> ncols n)
			 (format stream " ...")))
		  (princ ")" stream)
		  (unless (= i (1- m)) (format stream " "))
		  finally
		  (when (> nrows m)
		    (format stream " ..."))
		  (princ ")" stream))))
	    (t (princ "#M(??)" stream))))
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

(defmethod eye (n &optional m (type 'double-float))
  "Returns the nxn identity matrix.  The value is freshly allocated."
  (let ((result (make-instance (standard-matrix type) :nrows n :ncols (or m n)))
	(one (coerce 1 type)))
    (dotimes (i (if m (min n m) n) result)
      (setf (mref result i i) one))))

(defmethod zeros (n &optional m (type 'double-float))
  "Returns nxn or (if m is provided) nxm zeros.  The value is freshly
allocated."
  (make-instance (standard-matrix type) :nrows n :ncols (or m n)))

(defmethod ones (n &optional m (type 'double-float))
  "Returns nxn or (if m is provided) nxm ones.  The value is freshly
allocated."
  (unless m (setq m n))
  (make-instance (standard-matrix type) :nrows n :ncols m
		 :store (make-array (* n m) :element-type type
				    :initial-element (coerce 1 type))))

(defmethod mrandom ((n integer) &optional m (type 'double-float) (range 1.0))
  "Returns a random nxn or (if m is provided) nxm matrix.  The value is
freshly allocated."
  (unless m (setq m n))
  (fill-random! (zeros n m) (coerce range type)))

(defmethod diag (vec)
  "Returns a diagonal matrix with diagonal entries from vec."
  (let* ((n (vlength vec))
	 (result (make-instance (standard-matrix (element-type vec))
				:nrows n :ncols n)))
    (dotimes (i n result)
      (setf (mref result i i) (vref vec i)))))

(defun laplace-full-matrix (n)
  "Generates the matrix for a 1-dimensional Laplace problem discretized
with the 3-point stencil on a structured mesh."
  (let ((result (make-instance (standard-matrix 'double-float) :nrows n :ncols n)))
    (dotimes (i n)
      (setf (mref result i i) 2.0)
      (when (> i 0) (setf (mref result i (1- i)) -1.0))
      (when (< i (1- n)) (setf (mref result i (1+ i)) -1.0)))
    result))

(defmethod make-analog ((x standard-matrix))
  "Constructs a zero matrix with the same size as X."
  (zeros (nrows x) (ncols x) (element-type x)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Matrix combination
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmethod join ((x standard-matrix) (y standard-matrix) &optional (orientation :horizontal))
  "Joins standard matrices X and Y either horizontally or vertically
depending on the value of ORIENTATION."
  (flet ((layout (m n)
	   (make-instance (standard-matrix (element-type x))
			  :nrows m :ncols n)))
    (ecase orientation
      (:horizontal
       (join-horizontal! x y (layout (nrows x) (+ (ncols x) (ncols y)))))
      (:vertical
       (join-vertical! x y (layout (+ (nrows x) (nrows y)) (ncols x)))))))

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
(defmethod for-each-key-and-entry (fn (mat standard-matrix))
  (dotimes (i (nrows mat))
    (dotimes (j (ncols mat))
      (funcall fn i j (mref mat i j)))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Testing
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(eval-when (:load-toplevel :execute)
  ;; to have the #M reader macro working in its modified form this use is
  ;; wrapped in an EVAL-WHEN
  (defun test-standard-matrix ()
    ;;(setq *print-matrix* nil)
    (make-instance (standard-matrix 'integer) :content #2a((2 3) (4 5)))
    (describe (standard-matrix 'double-float))
    (standard-matrix 'single-float)
    (describe (eye 1))
    (eye 2)
    (zeros 1)
    (make-real-matrix `((,(cos 1.0))))
    (make-real-matrix #((2.0 3.0) (4.0 5.0)))
    (copy (eye 5))
    (norm #m((1.0 2.0)) 1)
    )
  )

;;; (fl.matlisp::test-standard-matrix)
(fl.tests:adjoin-test 'test-standard-matrix)