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

(in-package :fl.algebra)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Class definition and basic methods
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  
(defclass <sparse-tensor> ()
  ((rank :reader rank :initarg :rank :type fixnum
	 :documentation "Tensor rank.")
   (dimension :reader dimension :initarg :dimension :type (or null fixnum)
	      :documentation "Dimension of the first slot.")
   (indices :accessor indices :initarg :indices :type vector
	    :initform (make-array 0 :element-type 'fixnum :adjustable t)
	    :documentation "The (nonzero) indices of the first slot.")
   (entries :accessor entries :initarg :entries :type vector
	    :initform (make-array 0 :adjustable t)
	    :documentation "The (nonzero) entries of the first slot."))
  (:documentation "A general sparse tensor class which is implemented as a
sparse vector containing full-or sparse tensor entries."))

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

(defun in-pattern-p (tensor &rest indices)
  "Returns a getter/setter pair for the specified location."
  (if (null indices)
      t
      (whereas ((pos (position (first indices) (indices tensor))))
	(apply #'in-pattern-p (aref (entries tensor) pos) (cdr indices)))))
  
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
	     (adjust-array indices (1+ length) :initial-element -1)
	     (adjust-array entries (1+ length) :initial-element nil)
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

(defmethod for-each-key-and-entry ((func function) (tensor <sparse-tensor>))
  (case (rank tensor)
    (1 (loop for index across (indices tensor)
	     and entry across (entries tensor) do
	     (funcall func index entry)))
    (t (loop for index across (indices tensor)
	     and entry across (entries tensor) do
	     (for-each-key-and-entry (curry func index) entry)))))

(defmethod for-each-key ((func function) (tensor <sparse-tensor>))
  (case (rank tensor)
    (1 (loop for index across (indices tensor) do
	     (funcall func index)))
    (t (loop for index across (indices tensor)
	     and entry across (entries tensor) do
	     (for-each-key (curry func index) entry)))))

(defmethod for-each-entry ((func function) (tensor <sparse-tensor>))
  (case (rank tensor)
    (1 (loop for entry across (entries tensor) do
	     (funcall func entry)))
    (t (loop for entry across (entries tensor) do
	     (for-each-entry func entry)))))

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
	       (if (or (typep entry '<tensor>)  (typep entry '<sparse-tensor>))
		   (show-tensor entry (1+ level))))))
    (show-tensor tensor 0)
    (format t "~%")))

;;; Testing: (test-sparse-tensor)

(defun test-sparse-tensor ()
  (let ((tensi (make-instance '<sparse-tensor> :rank 2)))
    (setf (tensor-ref tensi 1 2) 1.0)
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
    (setf (tensor-ref tensi 1 5) 3.0)
    (describe (tensor-ref tensi 1)))
  )

(fl.tests:adjoin-test 'test-sparse-tensor)

