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

(defclass standard-matrix (<matrix> <store-vector>)
  ((nrows
    :initarg :nrows :accessor nrows :type fixnum
    :documentation "Number of rows in the matrix")
   (ncols
    :initarg :ncols :accessor ncols :type fixnum
    :documentation "Number of columns in the matrix")
   (store
    :reader store :initarg :store :type (simple-array * (*))
    :documentation "The actual storage for the matrix."))
  (:documentation "Basic matrix class."))

;;; if standard-matrix is considered as multiple vectors
(defmethod vlength ((mat standard-matrix))
  (nrows mat))
(defmethod multiplicity ((mat standard-matrix))
  (ncols mat))

(defmethod element-type (matrix)
  "Default method returns T."
  T)

(defmethod scalar-type (matrix)
  "Default method returns NUMBER."
  'NUMBER)

(declaim (inline indexing))
(defun indexing (i j nrows ncols)
  "Computes column-major indexing for compatibility with Matlisp/Fortran."
  (declare (ignore ncols))
  (declare (type positive-fixnum i j nrows ncols))
  (the fixnum (+ i (the fixnum (* nrows j)))))

(defmethod initialize-instance :after ((matrix standard-matrix) &key content &allow-other-keys)
  "Handles the content parameter and sets the store to a suitable vector.
If content is a 2d array, the dimensions can be deduced."
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
		       (eq (upgraded-array-element-type type)
			   (array-element-type store))))
	  (setq store (make-array n :element-type type)))
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

;;; (make-instance 'standard-matrix :content #2a((2 3) (4 5)))

(defun standard-matrix (type)
  (assert (subtypep type 'number))
  (let ((class-name (intern (format nil "~A" (list 'standard-matrix type))
			    :fl.matlisp)))
    (or (find-class class-name nil)
	(prog1
	    (eval `(defclass ,class-name (standard-matrix)
		    ((store :type (simple-array ,type (*))))))
	  (eval `(defmethod element-type ((matrix ,class-name))
		  ',type))
	  (eval `(defmethod scalar-type ((matrix ,class-name))
		  ',type))))))

(defun make-real-matrix (&rest args)
  "Generates a matrix with double-float entries."
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
       (let ((list (read stream nil (values) t)))
	 `(make-real-matrix ',list)))))

(defun make-real-vector (dim &optional (value 0.0))
  "Generates a dimx1matrix with double-float entries."
  (make-instance (standard-matrix 'double-float)
		 :nrows dim :ncols 1
		 :content (make-double-vec dim value)))


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

(defmethod vref ((matrix standard-matrix) i)
  "Note: access to matrix entries using this generic function is slow.
Therefore, specialized BLAS routines should be used whenever possible."
    (aref (slot-value matrix 'store) i))

(defmethod (setf vref) (value (matrix standard-matrix) i)
  "Writer for vref."
  (setf (aref (slot-value matrix 'store) i) value))

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
	(with-slots (nrows ncols store)
	    matrix
	  (cond
	    ((= ncols 1)
	     (format stream "#M(~{~A~^ ~})" (coerce store 'list)))
	    ((= nrows 1)
	     (format stream "#M((~{~A~^ ~}))" (coerce store 'list)))
	    (t (loop with m = nrows and n = ncols
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
		     (princ ")" stream))))))
      (print-unreadable-object (matrix stream :type t :identity t))))

(defmethod print-object ((matrix standard-matrix) stream)
  (print-matrix matrix stream))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Copy
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmethod copy ((matrix standard-matrix))
  "Copy constructor for standard-matrix."
  (make-instance
   (class-of matrix)
   :nrows (slot-value matrix 'nrows)
   :ncols (slot-value matrix 'ncols)
   :content (slot-value matrix 'store)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Special matrices
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; Constructors for special standard matrices

(defmethod eye (n &optional m (type 'double-float))
  "Returns the nxn identity matrix.  The value is freshly allocated."
  (let ((result (make-instance (standard-matrix type) :nrows n :ncols n))
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

(defmethod diag (vec)
  "Returns a diagonal matrix with diagonal entries from vec."
  (let* ((n (vlength vec))
	 (result (make-instance (standard-matrix (element-type vec))
				:nrows n :ncols n)))
    (dotimes (i n result)
      (setf (mref result i i) (vref vec i)))))

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

(defmethod for-each-entry ((func function) (mat standard-matrix))
  (dotimes (i (nrows mat))
    (dotimes (j (ncols mat))
      (funcall func (mref mat i j)))))
(defmethod for-each-row-key ((fn function) (mat standard-matrix))
  "Loop through row keys."
  (dotimes (i (nrows mat))
    (funcall fn i)))
(defmethod for-each-key-in-row ((fn function) (mat standard-matrix) row-key)
  "Loop through column keys in row."
  (declare (ignore row-key))
  (dotimes (i (ncols mat))
    (funcall fn i)))
(defmethod for-each-col-key ((fn function) (mat standard-matrix))
  "Loop through column keys."
  (dotimes (i (ncols mat))
    (funcall fn i)))
(defmethod for-each-key-and-entry ((func function) (mat standard-matrix))
  (dotimes (i (nrows mat))
    (dotimes (j (ncols mat))
      (funcall func i j (mref mat i j)))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Testing
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun test-standard-matrix ()
  (describe (standard-matrix 'double-float))
  (standard-matrix 'single-float)
  (describe (eye 1))
  (eye 2)
  (zeros 1)
  (make-real-matrix `((,(cos 1.0))))
  (make-real-matrix #((2.0 3.0) (4.0 5.0)))
  (make-instance 'standard-matrix :content #((2.0 3.0) (4.0 5.0)))
  (copy (eye 5))
  (norm #m((1.0 2.0)) 1))

;;; (test-standard-matrix)
