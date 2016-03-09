;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; compressed.lisp - Compressed sparse storage scheme
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;
;;; Copyright (C) 2006 Nicolas Neuss, University of Karlsruhe.
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

(defclass compressed-pattern ()
  ((sizes :initarg :sizes :accessor sizes
	  :documentation "Vector of matrix sizes, at the moment only length
2 is allowed here.  The first is the dimension which is not compressed, the
second is the dimension which gets compressed.")
   (orientation :accessor orientation :initform :column :type (member :row :column) :initarg :orientation
		:documentation "Denotes if rows or columns are compressed.")
   (starts :accessor starts :initarg :starts
	   :documentation "Vector with start indices of compressed columns/rows.")
   (indices :accessor indices :initarg :indices
	    :documentation "Vector with compressed row/column indices.")
   (offsets :accessor offsets :initarg :offsets :initform nil
	    :documentation "Vector of offsets.  This is only non-nil, if
the pattern supports identification."))
  (:documentation "A compressed sparse pattern.  Note: we use int32 vectors
for @slot{starts} and @slot{indices}, such that they do not have to be
copied for a call to the alien sparse solvers."))

(defmethod initialize-instance :after
    ((pattern compressed-pattern) &key ((:pattern pattern-list)) &allow-other-keys)
  "This is a more elaborate compressed-pattern constructor.  A sparse
matrix of the form @math{ | * 0 0 0 a | | 0 a 0 0 0 | } can be described by
its sizes as a vector \(dimension of non-compressed, dimension/compressed)
together with the pattern '( ((* . 0) (a . 4)) ((a . 1)) ). Here, * means a
non-identified value.  Other symbols can be used to identify entries."
  ;; handle pattern parameter
  (when pattern-list
    (when (some (lambda (name)
		  (and (slot-boundp pattern name)
		       (slot-value pattern name)))
		'(starts indices offsets))
      (error "Either pattern or starts/indices/offsets should be given."))
    ;; initialize slots from pattern parameter
    (with-slots (sizes starts indices offsets) pattern
      (unless (= (length pattern-list) (aref sizes 0))
	(error "First entry of sizes and the length of the pattern do not match."))
      (multiple-value-setq (starts indices offsets)
	(loop with table = (make-hash-table)
	      and offset = 0
	      for row in pattern-list
	      do (assert (reduce #'< row :initial-value -1 :key #'cdr)) 
	      appending (mapcar #'cdr row) into indices
	      summing (length row) into nr-entries
	      appending
	      (mapcar (lambda (entry)
			(cond ((eq (car entry) '*)
			       (prog1 offset (incf offset)))
			      ((gethash (car entry) table))
			      (t (prog1
				     (setf (gethash (car entry) table) offset)
				   (incf offset)))))
		      row)
	      into offsets
	      collect nr-entries into starts
	      finally (return (values (cons 0 starts) indices offsets))))))
  ;; ensure that all slots are of the correct type
  (with-slots (sizes starts indices offsets) pattern
    (assert (= (length starts) (1+ (aref sizes 0))))
    (assert (= (length indices) (number-of-nonzero-entries pattern)))
    (setq starts (coerce starts 'int-vec))
    (setq indices (coerce indices 'int-vec))
    (setq sizes (coerce sizes 'fixnum-vec))
    (when offsets 
      (assert (= (length indices) (length offsets)))
      (setq offsets (coerce offsets 'int-vec)))))

(defgeneric transposed-pattern (pattern)
  (:documentation "Transpose a sparse matrix pattern.")
  (:method ((pattern compressed-pattern))
      (with-slots (sizes starts indices orientation) pattern
        (assert (= 2 (length sizes)))
        (make-instance 'compressed-pattern :sizes sizes
                                           :starts starts :indices indices
                                           :orientation (ecase orientation
                                                          (:row :column)
                                                          (:column :row))))))

(defgeneric number-of-nonzero-entries (pattern)
  (:documentation "Number of nonzero entries of a sparse matrix pattern.")
  (:method ((pattern compressed-pattern))
      (with-slots (starts) pattern
        (elt starts (1- (length starts))))))

(with-memoization (:id 'full-compressed-pattern)
  (defun full-compressed-pattern (nrows ncols &optional (orientation :column))
    "Returns a full compressed pattern."
    (when (eq orientation :column)
      (rotatef nrows ncols))
    (memoizing-let ((nrows nrows) (ncols ncols) (orientation orientation))
      (make-instance
       'compressed-pattern
       :sizes (vector nrows ncols)
       :starts (coerce (loop for i from 0 upto nrows
                          collect (* i ncols))
                       'int-vec)
       :indices (coerce (loop for i from 0 below (* nrows ncols)
                           collect (mod i ncols))
                        'int-vec)
       :orientation orientation))))

;;; for backward compatibility
(defun full-ccs-pattern (nrows ncols)
  (full-compressed-pattern nrows ncols :column))
(defun full-crs-pattern (nrows ncols)
  (full-compressed-pattern nrows ncols :row))

(defmethod in-pattern-p ((pattern compressed-pattern) &rest rest)
  (with-slots (starts indices orientation) pattern
    (destructuring-bind (i j)
        (ecase orientation
          (:row rest)
          (:column (reverse rest)))
      (loop for k from (aref starts i) below (aref starts (1+ i))
           thereis (= j (aref indices k))))))

(defclass compressed-matrix (<matrix>)
  ((pattern :reader pattern :initarg :pattern :type compressed-pattern
	  :documentation "A compressed pattern."))
  (:documentation "A compressed sparse matrix.  This is an abstract class
which is made concrete by mixing it with a store-vector containing the
entries."))

(defmethod initialize-instance :after ((cm compressed-matrix) &key &allow-other-keys)
  (assert (typep cm 'store-vector))
  (with-slots (store pattern) cm
    (if (slot-boundp cm 'store)
	(assert (= (length (store cm))
		   (number-of-nonzero-entries (pattern cm))))
	(setf store (zero-vector (number-of-nonzero-entries pattern) (element-type cm))))))

(inlining
 (defun find-compressed-offset (cm i j)
   (with-slots (pattern) cm
     (with-slots (starts indices orientation) pattern
       (when (eq orientation :row) (rotatef i j))
       (position i indices :start (aref starts j) :end (aref starts (1+ j)))))))

(defmethod mref ((cm compressed-matrix) i j)
  "This method is at the moment relatively inefficient because it ignores
any ordering."
  (let ((offset (find-compressed-offset cm i j)))
    (if offset
	(aref (store cm) offset)
	(coerce 0 (element-type cm)))))

(defmethod nrows ((cm compressed-matrix))
  (with-slots (sizes orientation) (pattern cm)
    (aref sizes (ecase orientation (:row 0) (:column 1)))))

(defmethod ncols ((cm compressed-matrix))
  (with-slots (sizes orientation) (pattern cm)
    (aref sizes (ecase orientation (:row 1) (:column 0)))))

(defmethod in-pattern-p ((cm compressed-matrix) &rest indices)
  (apply #'in-pattern-p (pattern cm) indices))

(defmethod make-domain-vector-for ((cm compressed-matrix) &optional (multiplicity 1))
  (make-instance (standard-matrix (element-type cm))
		 :nrows (ncols cm) :ncols multiplicity))

(defmethod make-image-vector-for ((cm compressed-matrix) &optional (multiplicity 1))
  (make-instance (standard-matrix (element-type cm))
		 :nrows (nrows cm) :ncols multiplicity))

(with-memoization (:type :global :size 4 :id 'compressed-matrix)
  (defun compressed-matrix (type)
    "Construct a compressed sparse matrix with entries of @arg{type}."
    (memoizing-let ((type type))
      (assert (subtypep type 'number))
      (fl.amop:find-programmatic-class
       (list 'compressed-matrix (store-vector type))))))

(defun make-full-compressed-matrix (nrows ncols &optional (orientation :column))
  (make-instance (compressed-matrix 'double-float)
		 :pattern (full-compressed-pattern nrows ncols orientation)))

(defun make-full-crs-matrix (nrows ncols)
  (make-full-compressed-matrix nrows ncols :row))

(defmethod transpose ((cm compressed-matrix))
  "A compressed matrix can be transposed easily by transposing its
pattern."
  (with-slots (pattern store) cm
    (make-instance (class-of cm) :pattern (transposed-pattern pattern)
		   :store store)))

(defmethod for-each-entry-and-key ((fn function) (x compressed-matrix))
  (let* ((store (store x))
	 (pattern (pattern x))
	 (sizes (sizes pattern))
	 (starts (starts pattern))
	 (indices (indices pattern))
	 (offsets (offsets pattern))
	 (orientation (orientation pattern)))
    (declare (type (simple-array * (*)) store))
    (declare (type (member :row :column) orientation))
    (declare (type int-vec starts indices))
    (declare (type (or null int-vec) offsets))
    (dotimes (i (aref sizes 0))
      (loop for k of-type int from (aref starts i) below (aref starts (1+ i)) do
	   (let* ((j (aref indices k))
		  (l (if offsets (aref offsets k) k)))
	     (declare (type int j l))
	     (if (eq orientation :row)
		 (funcall fn (aref store l) i j)
		 (funcall fn (aref store l) j i)))))))

(defgeneric compressed->matlisp (cm)
  (:documentation "Converts a compressed matrix into matlisp format.")
  (:method ((cm compressed-matrix))
      (lret ((result (zeros (nrows cm) (ncols cm) (element-type cm))))
        (for-each-entry-and-key
         (lambda (value i j)
           (setf (mref result i j) value))
         cm))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; GEMM!
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmethod m*-product-instance ((A compressed-matrix) (x vector))
  (zero-vector (nrows A) (array-element-type x)))

(inlining
 (defun compressed-gemm! (alpha x y beta z x-transposed y-transposed)
   (unless (= beta (coerce 1 (type-of beta))) (scal! beta z))
   (let ((ny (nrows y))
	 (nz (nrows z))
	 (multiplicity (ncols y))
	 (storey (store y))
	 (storez (store z)))
     #+sbcl (declare (sb-ext:muffle-conditions sb-ext:code-deletion-note))
     (labels ((zref (j l)
		(aref storez (+ j (* l nz))))
	      ((setf zref) (value j l)
		(setf (aref storez (+ j (* l nz))) value))
	      (yref (j l)
		(aref storey (if y-transposed
				 (+ l (* j ny))
				 (+ j (* l ny)))))
	      (base-op (entry i j)
		(dotimes (l multiplicity)
		  (setf (zref i l)
			(+ (* alpha entry (yref j l))
			   (zref i l))))))
       (for-each-entry-and-key 
	(lambda (xc i j)
	  (if x-transposed
	      (base-op xc j i)
	      (base-op xc i j)))
	x)))
   z))

(defmethod gemm-nn! (alpha (x compressed-matrix) (y standard-matrix) beta (z standard-matrix))
  (compressed-gemm! alpha x y beta z nil nil))
(defmethod gemm-nt! (alpha (x compressed-matrix) (y standard-matrix) beta (z standard-matrix))
  (compressed-gemm! alpha x y beta z nil t))
(defmethod gemm-tn! (alpha (x compressed-matrix) (y standard-matrix) beta (z standard-matrix))
  (compressed-gemm! alpha x y beta z t nil))
(defmethod gemm-tt! (alpha (x compressed-matrix) (y standard-matrix) beta (z standard-matrix))
  (compressed-gemm! alpha x y beta z t t))

#+(or)
(let* ((n 1000)
       (A (time (fl.matlisp::five-point-stencil-matrix n n)))
       (b (ones (* n n) 1))
       (x (ones (* n n) 1)))
  (time (gemm-nn! 1.0 A b 1.0 x))
  )

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; GESV!
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defparameter *default-cm-solver*
  (if fl.start::*superlu-library*
      'fl.alien::superlu
      (when fl.start::*umfpack-library*
	'fl.alien::umfpack))
  "Default solver for the CM format.  At the moment this can be UMFPACK or
SuperLU.")

(defmethod gesv! ((mat compressed-matrix) (vec standard-matrix))
  "Solve the system by calling an external sparse solver."
  (if *default-cm-solver*
      (with-slots (sizes starts indices orientation)
	  (pattern mat)
	(let ((nrows (aref sizes 1))
	      (ncols (aref sizes 0)))
	  (assert (= nrows (nrows vec)))
	  (unless (zerop
                   (funcall *default-cm-solver*
                            nrows ncols (number-of-nonzero-entries (pattern mat))
                            starts indices (store mat)
                            (ncols vec) (store vec) (store vec)
                            (ecase orientation (:column 0) (:row 1))))
            (error "External direct solver did not succeed."))
	  vec))
      (gesv! (compressed->matlisp mat) vec)))
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Test direct solvers on the CM scheme
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun five-point-stencil-matrix (nx ny &optional (orientation :column))
  "Generate a CM matrix for the five-point stencil on a grid with the
given number of active grid points in each dimension."
  (declare (type (integer 0 10000) nx ny))
  (let* ((nrows (* nx ny))
	 (ncols nrows)
	 (nnz (- (* 5 nrows) (* 2 (+ nx ny))))
	 (starts (zero-vector (1+ nrows) '(signed-byte 32)))
	 (indices (zero-vector nnz '(signed-byte 32)))
	 (store (zero-vector nnz 'double-float))
	 (pos 0))
    (declare (type (integer 0 100000000) pos))
    (dotimes (j ny)
      (declare (type (integer 0 10000) j))
      (dotimes (i nx)
	(let ((k (+ i (* j nx))))
	  (declare (optimize speed))
	  (setf (aref starts k) pos)
	  (flet ((connect (l value)
		   (setf (aref store pos) value)
		   (setf (aref indices pos) l)
		   (incf pos)))
	    (when (plusp j) (connect (- k nx) -1.0))
	    (when (plusp i) (connect (- k 1) -1.0))
	    (connect k 4.0)
	    (when (< i (1- nx)) (connect (+ k 1) -1.0))
	    (when (< j (1- ny)) (connect (+ k nx) -1.0))))))
    (assert (= pos nnz))
    (setf (aref starts nrows) pos)
    ;; return matrix
    (make-instance
     (compressed-matrix 'double-float)
     :pattern (make-instance 'compressed-pattern :sizes (vector nrows ncols)
			     :starts starts :indices indices
			     :orientation orientation)
     :store store)))


(defun direct-solver-performance-test (solver n)
  (when (evenp n) (setf n (1- n)))
  (let* ((cm (five-point-stencil-matrix n n))
	 (rhs (ones (nrows cm) 1)))
    (time
     (let ((*default-cm-solver* solver))
       (progn (gesv! cm rhs)
	      (/ (vref rhs (floor (nrows cm) 2))
		 (expt (1+ n) 2)))))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Testing
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun test-compressed ()
  (when fl.start::*superlu-library*
    (direct-solver-performance-test 'fl.alien::superlu 50))
  (when fl.start::*umfpack-library*
    (direct-solver-performance-test 'fl.alien::umfpack 100))
  (let ((*print-matrix* t))
    (princ (compressed->matlisp (five-point-stencil-matrix 4 4))))
  (let* ((n 300)
	 (mat (five-point-stencil-matrix n n)))
    (time
     (for-each-entry-and-key (lambda (x i j)
			       (declare (ignore x i j))
			       nil)
			     mat)))
  (make-instance (store-vector 'single-float)
		 :store (zero-vector 1 'single-float))
  (let* ((pattern (make-instance
		   'compressed-pattern :sizes #(1 1)
		   :starts (int-vec 0 1)
		   :indices (int-vec 0)))
	 (cm (make-instance
	       (compressed-matrix 'double-float) :pattern pattern :store #d(2.0)))
	 (rhs #m(1.0)))
    (mref cm 0 0)
    (gesv! cm rhs)
    (transpose cm))
  (transpose (five-point-stencil-matrix 2 2))
  (make-full-crs-matrix 3 2)
  (make-instance 'compressed-pattern
		 :sizes #(2 2) :orientation :row
		 :pattern '( ((a . 0))  ((a . 1)) ))
  )

;;; (test-compressed)
(fl.tests:adjoin-test 'test-compressed)


