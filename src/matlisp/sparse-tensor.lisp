;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; sparse-tensor.lisp - sparse tensors of arbitrary rank
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
;;; Class definition and basic methods
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass <sparse-tensor> (tensor)
  ((rank :reader rank :initarg :rank :type fixnum
	 :documentation "Tensor rank.")
   (indices :accessor indices :initarg :indices :type vector
	    :initform (make-array 0 :element-type 'fixnum :adjustable t)
	    :documentation "The (nonzero) indices of the first slot.")
   (entries :accessor entries :initarg :entries :type vector
	    :initform (make-array 0 :adjustable t)
	    :documentation "The (nonzero) entries of the first slot."))
  (:documentation "A general sparse tensor class which is implemented as a
sparse vector containing full-or sparse tensor entries."))

(defmethod multiplicity ((tensor <sparse-tensor>))
  (assert (= (rank tensor) 1))
  1)

(defmethod initialize-instance :after ((tensor <sparse-tensor>) &key initial-contents &allow-other-keys)
  (loop for (entry . indices) in initial-contents do
	(setf (apply #'tensor-ref tensor indices) entry)))

(defun sparse-tensor (contents)
  "Constructor for @class{sparse-tensor}."
  (make-instance '<sparse-tensor>
		 :rank (1- (length (first contents)))
		 :initial-contents contents))

(defmethod make-analog ((tensor <sparse-tensor>))
  (with-slots (rank indices entries) tensor
    (make-instance
     '<sparse-tensor> :rank rank :indices indices
     :entries (map 'vector (lambda (entry) (make-analog entry))
		   entries))))

;;; sparse tensors of rank 2 behave as matrices

(defun diagonal-sparse-tensor (values &optional ncomps)
  "Constructs a sparse tensor of rank 2 where @arg{values} are the diagonal
entries.  If @arg{ncomps} is given then the tensor dimension is nxn with
each diagonal entry being @arg{values}."
  (lret ((result (make-instance '<sparse-tensor> :rank 2)))
    (dotimes (i (or ncomps (length values)))
      (setf (mref result i i)
	    (if ncomps values (aref values i))))))

(defun find-first-position>= (indices index)
  "Finds the position for index, which is the first position which is >=
index in the sorted list."
  (declare (type fixnum index)
	   (values fixnum))
  (labels ((find-position (a b)
	     (declare (type fixnum a b))
	     (let ((c (ceiling (+ a b) 2)))
	       (if (= c b)
		   b			; return position, even if index was not found
		   (let ((index2 (aref indices c)))
		     (cond ((> index2 index) (find-position a c))
			   ((= index2 index) c)
			   (t (find-position c b))))))))
    (let* ((length (length indices))
	   (length-1 (1- length)))
      (if (zerop length)
	  0
	  (let ((index-0 (aref indices 0)))
	    (if (<= index index-0)
		0
		(let ((index-e (aref indices length-1)))
		  (cond ((= index index-e) length-1)
			((> index index-e) length)
			(t (find-position 0 length-1))))))))))

(defgeneric in-pattern-p (tensor &rest indices)
  (:documentation "Returns T, if the indices are in the nonzero pattern.")
  (:method (object &rest indices)
    "The default method returns T iff no indices are given."
    (declare (ignore object))
    (null indices))
  (:method ((tensor array) &rest indices)
    (every #'< indices (array-dimensions tensor)))
  (:method ((vec <vector>) &rest indices)
    (destructuring-bind (index) indices
      (and (plusp index) (< index (nrows vec)))))
  (:method ((mat <matrix>) &rest indices)
    (destructuring-bind (i j) indices
      (and (plusp i) (plusp j)
	   (< i (nrows mat)) (< j (ncols mat)))))
  (:method ((tensor full-tensor) &rest indices)
    (every #'< indices (dimensions tensor)))
  )

(defmethod in-pattern-p ((tensor <sparse-tensor>) &rest indices)
  (if (null indices)
      t					; indices are in this pattern
      (whereas ((pos (position (first indices) (indices tensor))))
	(let ((rest (rest indices)))
	  (if (null rest)
	      pos
	      (apply #'in-pattern-p (aref (entries tensor) pos) rest))))))

(definline get-place-in-tensor (tensor &rest search-indices)
  "Returns a getter/setter pair for the specified location."
  (labels
      ((get-place (tensor search-indices)
	 (let* ((indices (indices tensor))
		(entries (entries tensor))
		(index (first search-indices))
		(pos (find-first-position>= indices index))
		(length (length indices)))
	   (unless (and (< pos length)
			(= (aref indices pos) index))
	     ;; generate new (sparse!) place
	     (setq indices (adjust-array indices (1+ length) :initial-element -1))
	     (setq entries (adjust-array entries (1+ length) :initial-element nil))
	     (loop for i from length above pos do
		   (setf (aref indices i) (aref indices (1- i)))
		   (setf (aref entries i) (aref entries (1- i))))
	     (setf (aref indices pos) index)
	     (if (cdr search-indices)
		 (setf (aref entries pos)
		       (make-instance '<sparse-tensor> :rank (1- (rank tensor))))))
	   (if (null (cdr search-indices))
	       ;; getter/setter pair
	       (values #'(lambda () (aref entries pos))
		       #'(lambda (value) (setf (aref entries pos) value)))
	       ;; continue search
	       (get-place (aref entries pos) (cdr search-indices))))))
    (if (null search-indices)
	;; getter for tensor
	(values #'(lambda () tensor)
		#'(lambda (value)
		    (declare (ignore value))
		    (error "No setter possible.")))
	(get-place tensor search-indices))))

(defmethod tensor-ref ((tensor <sparse-tensor>) &rest indices)
  (funcall (apply #'get-place-in-tensor tensor indices)))
(defmethod (setf tensor-ref) (value (tensor <sparse-tensor>) &rest indices)
  (funcall (nth-value 1 (apply #'get-place-in-tensor tensor indices))
	   value))

;;; for rank 2 tensors we allow mref access
(defmethod mref ((tensor <sparse-tensor>) i j)
  (tensor-ref tensor i j))
(defmethod (setf mref) (value (tensor <sparse-tensor>) i j)
  (setf (tensor-ref tensor i j) value))

(defmethod for-each-key ((func function) (tensor <sparse-tensor>))
  (labels ((helper (reversed-index tensor)
	     (case (rank tensor)
	       (1 (loop for index across (indices tensor) do
			(funcall func (reverse (cons index reversed-index)))))
	       (t (loop for index across (indices tensor)
			and entry across (entries tensor) do
			(helper (cons index reversed-index) entry))))))
    (helper () tensor)))

(defmethod for-each-entry ((func function) (tensor <sparse-tensor>))
  (case (rank tensor)
    (1 (loop for entry across (entries tensor) do
	     (funcall func entry)))
    (t (loop for entry across (entries tensor) do
	     (for-each-entry func entry)))))

(defmethod for-each-entry-and-key ((func function) (tensor <sparse-tensor>))
  (let ((full-index (make-list (rank tensor) :initial-element 0)))
    (labels ((helper (tensor current-index)
	       (case (rank tensor)
		 (0 (funcall func tensor))
		 (1 (loop for index across (indices tensor)
			  and entry across (entries tensor) do
			  (setf (car current-index) index)
			  (apply func entry full-index)))
		 (t (loop for index across (indices tensor)
			  and entry across (entries tensor) do
			  (setf (car current-index) index)
			  (helper entry (cdr current-index)))))))
      (helper tensor full-index))))

(defmethod tensor-for-each ((func function) (tensor <sparse-tensor>) &key (job :both) depth)
  (unless depth (setq depth (rank tensor)))
  (cond
    ((not (plusp depth))
     (ecase job
       (:entry (funcall func tensor))
       (:index (funcall func))
       (:both (funcall func tensor))))
    ((= depth 1)
     (loop for index across (indices tensor)
	   and entry across (entries tensor) do
	   (ecase job
	     (:entry (funcall func entry))
	     (:index (funcall func index))
	     (:both (funcall func index entry)))))
    (t
     (loop for index across (indices tensor)
	   and entry across (entries tensor) do
	   (if (eq job :entry)
	       (tensor-for-each func entry :job :entry :depth (1- depth))
	       (tensor-for-each (curry func index) entry :job job :depth (1- depth)))))))

(defmethod tensor-for-each ((func function) (mat standard-matrix)
			    &key (job :both) &allow-other-keys)
  (dotimes (i (nrows mat))
    (dotimes (j (ncols mat))
      (ecase job
	(:entry (funcall func (mref mat i j)))
	(:index (funcall func i j))
	(:both (funcall func i j (mref mat i j)))))))

(defmacro dotensor ((args tensor &key depth) &body body)
  "Usage:
\(dotensor (entry tensor :depth 1) ...)
\(dotensor ((index1 ... . entry) tensor :depth 1) ...)
\(dotensor ((index1 ...) tensor :depth 1) ...)"
  (cond ((atom args)
	 `(tensor-for-each #'(lambda (,args) ,@body) ,tensor :job :entry :depth ,depth))
	((null (cdr (last args)))
	 `(tensor-for-each #'(lambda ,args ,@body)
	   ,tensor :job :index :depth ,(or depth (length args))))
	(t `(tensor-for-each #'(lambda ,(nconc (butlast args 1)
					       (list (car (last args)))
					       (when (cdr (last args))
						 (list (cdr (last args)))))
				 ,@body)
	     ,tensor :job :both :depth ,(or depth (1+ (length (butlast args 1))))))))

(defmethod show ((tensor <sparse-tensor>) &key &allow-other-keys)
  (labels ((show-tensor (tensor level)
	     (dotensor ((index . entry) tensor :depth 1)
	       (if (= level 0)
		   (format t "~&")
		   (format t "~&>~VT" (* level 4)))
	       (format t "~D  --> ~A" index entry)
	       (if (or (typep entry 'full-tensor)  (typep entry '<sparse-tensor>))
		   (show-tensor entry (1+ level))))))
    (show-tensor tensor 0)
    (terpri)
    tensor))

(defmethod ensure-matlisp ((tensor <sparse-tensor>) &optional (type :column))
  (with-slots (indices) tensor
    (ecase (rank tensor)
      (1 (lret* ((n (length (indices tensor)))
		 (result (ecase type
			   (:column (zeros n 1))
			   (:row (zeros 1 n)))))
	   (dotimes (i n)
	     (assert (= i (aref indices i)))
	     (setf (vref result i)
		   (tensor-ref tensor i))))))))

(defmethod m* ((mat <sparse-tensor>) vec)
  (lret ((result (make-instance '<sparse-tensor> :rank 1)))
    (for-each-entry-and-key
     (lambda (entry i j)
       (when (in-pattern-p vec j)
	 (let ((prod (* entry (tensor-ref vec j))))
	   (if (in-pattern-p result i)
	       (incf (tensor-ref result i) prod)
	       (setf (tensor-ref result i) prod)))))
     mat)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; Testing:

(defun test-sparse-tensor ()
  (let ((tensi (sparse-tensor `((,(eye 1) 0)))))
    (for-each-entry-and-key (lambda (value k)
			      (format t "entry=~A k=~A~%" value k))
			    tensi)
    (show (tensor-map #'identity tensi))
    (show (m- tensi tensi)))
  (let ((tensi (make-instance '<sparse-tensor> :rank 2)))
    (setf (tensor-ref tensi 1 2) 1.0)
    (let ((sparse-vec (sparse-tensor `((9.0 2))))
	  (vec #(1.0 2.0 3.0)))
      (show (m* tensi vec))
      (show (m* tensi sparse-vec))
      (show (axpy 1.0 (m* tensi vec) sparse-vec)))
    (show tensi)
    (assert (in-pattern-p tensi 1 2))
    (assert (not (in-pattern-p tensi 2 2)))
    (total-entries tensi)
    (dotensor (entry tensi)
      (format t "~A~%" entry))
    (dotensor ((i) tensi) (format t "~A~%" i))
    (dotensor ((i j . entry) tensi) (format t "~A ~A ~A ~%" i j entry))
    (describe (tensor-ref tensi 1))
    (setf (tensor-ref tensi 1 2) 3.0)
    (describe (tensor-ref tensi 1))
    (assert (= (tensor-ref (tensor-ref tensi 1) 2)
	       (tensor-ref tensi 1 2)
	       (vref tensi '(1 2))))
    (setf (tensor-ref tensi 1 5) 4.0)
    (vref tensi '(1 5))
    (show tensi)
    (show (make-analog tensi))
    (show (copy tensi))
    (show (axpy 1.0 tensi tensi)))
  (let ((tensi (diagonal-sparse-tensor #(3.0 5.0))))
    (show tensi)
    (assert (mzerop (gemm 1.0 tensi #(1.0 2.0) -1.0 #(3.0 10.0)))))
  (let ((x (make-array 2 :initial-element (zeros 1)))
	(y (make-array 2 :initial-element (ones 1)))
	(sigma (diagonal-sparse-tensor (vector #m(0.0) #m(1.0)))))
    (gemm 1.0 sigma x 1.0 y))
  (let ((sigma (diagonal-sparse-tensor (vector 0.0 0.0)))
	(x (double-vec 0.0 0.0))
	(y (double-vec 0.0 0.0)))
    (gemm 1.0 sigma x 1.0 y))
  (let ((sigma (sparse-tensor '((0.0 0) (1.0 1)))))
    (ensure-matlisp sigma :row))
  )

;;; (test-sparse-tensor)
(fl.tests:adjoin-test 'test-sparse-tensor)

