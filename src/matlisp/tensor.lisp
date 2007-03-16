;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; tensor.lisp
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;
;;; Copyright (C) 2003-2005 Nicolas Neuss, University of Heidelberg.
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

;;;; Provides an extension of vectors and (matlisp) matrices operations to
;;;; tensors of arbitrary dimension.

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Class definition and basic methods
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass tensor ()
  ()
  (:documentation "Tensor superclass."))

;;; Reader and writers
(defgeneric tensor-ref (tensor &rest indices)
  (:documentation "Reader for a tensor entry.")
  (:method (object &rest indices)
     (declare (ignore object indices))
     (assert "Not implemented."))
  (:method :around (object &rest indices)
	   (if (null indices)
	       object
	       (call-next-method))))

(defgeneric (setf tensor-ref) (value x &rest indices)
  (:documentation "Writer for a tensor entry."))

(defgeneric tensor-map (func tensor)
  (:documentation "Maps @arg{tensor} with @arg{func} to a tensor of the
same type.")
  (:method (func tensor)
    "Default method."
    (lret ((result (make-analog tensor)))
      (for-each-entry-and-key
       (lambda (entry &rest indices)
	 (setf (apply #'tensor-ref result indices)
	       (funcall func entry)))
       tensor))))

;;; for vectors

(defmethod tensor-ref ((vec vector) &rest indices)
  (apply #'tensor-ref (aref vec (car indices)) (cdr indices)))

(defmethod (setf tensor-ref) (value (vec vector) &rest indices)
  (destructuring-bind (first . rest) indices
    (if rest
	(setf (apply #'tensor-ref (aref vec first) rest) value)
	(setf (vref vec first) value))))

;;; for matlisp matrices
(defmethod tensor-ref ((mat standard-matrix) &rest indices)
  (if (single? indices)
      (vref mat (first indices))
      (apply #'mref mat indices)))
      
(defmethod (setf tensor-ref) (value (mat standard-matrix) &rest indices)
  (if (single? indices)
      (setf (vref mat (first indices)) value)
      (setf (apply #'mref mat indices) value)))

(defclass full-tensor (tensor)
  ((dimensions :reader dimensions :initarg :dimensions :type fixnum-vec
	       :documentation "The dimensions of the tensor.")
   (offset0 :reader offset0 :initform 0 :initarg :offset0 :type fixnum
	    :documentation "An initial offset into the store-vector which
defaults to 0.")
   (offsets :reader offsets :initarg :offsets :type fixnum-vec
	    :documentation "The offsets for the different dimensions.  This
is internal information computed at tensor construction time."))
  (:documentation "Mixin for full tensors."))

(defun dimensions->offsets (dims)
  (coerce (loop	for k from 0
	     and dim across dims
	     and offset = 1 then (* offset dim) collect offset)
	  'fixnum-vec))

(defun full-tensor (type)
  "Construct a full tensor with entries of @arg{type}."
  (fl.amop:find-programmatic-class
   (list 'full-tensor (store-vector type))
   (intern (format nil "~A" (list 'FULL-TENSOR type)) "FL.MATLISP")))

(defmethod initialize-instance ((tensor full-tensor) &key &allow-other-keys)
  (call-next-method)
  (assert (typep tensor 'static-store-vector))
  (with-slots (store dimensions offsets) tensor
    (unless (slot-boundp tensor 'store)
      (setf store (zero-vector (reduce #'* dimensions) (element-type tensor))))
    (unless (slot-boundp tensor 'offsets)
      (setf offsets (dimensions->offsets dimensions)))))

(defun make-real-tensor (dimensions)
  "Generates an instance of a tensor with DOUBLE-FLOAT entries and the
given @arg{dimensions}."
  (make-instance (full-tensor 'double-float) :dimensions dimensions))

(defmethod make-analog ((tensor full-tensor))
  (make-instance (class-of tensor) :dimensions (dimensions tensor)))

(defun list->real-tensor (tree)
  (assert (tree-uniformp tree))
  (let ((tensor (make-real-tensor (coerce (tree-uniform-number-of-branches tree) 'fixnum-vec)))
	(k -1))
    (on-leaves #'(lambda (leaf) (setf (aref (store tensor) (incf k)) leaf))
	       tree)
    tensor))

(definline compute-offset (index offset0 offsets)
  (declare (type fixnum offset0))
  (declare (type fixnum-vec index offsets))
  (let ((sum offset0))
    (declare (type fixnum sum))
    (dotimes (i (length offsets) sum)
      (setq sum (the fixnum (+ sum (the fixnum (* (aref offsets i) (aref index i)))))))
    sum))

(defmethod tensor-ref ((tensor full-tensor) &rest indices)
  (let ((index (if (and (= 1 (length indices)) (typep (car indices) 'vector))
		   (coerce (car indices) 'fixnum-vec)
		   (coerce indices 'fixnum-vec))))
    (aref (store tensor) (compute-offset index (offset0 tensor) (offsets tensor)))))

(defmethod (setf tensor-ref) (value (tensor full-tensor) &rest indices)
  (let ((index (if (and (= 1 (length indices)) (typep (car indices) 'fixnum-vec))
		   (car indices)
		   (coerce indices 'fixnum-vec))))
    (setf (aref (store tensor) (compute-offset index (offset0 tensor) (offsets tensor)))
	  value)))

(defmethod vref ((tensor tensor) indices)
  (apply #'tensor-ref tensor indices))
(defmethod (setf vref) (value (tensor tensor) indices)
  (setf (apply #'tensor-ref tensor indices) value))

(defmethod rank ((tensor full-tensor)) (length (dimensions tensor)))
(defun make-tensor-index (tensor) (make-fixnum-vec (rank tensor)))

(defparameter *print-tensor*
  nil
  "Maximum number of columns and/or rows to print.  Set this to NIL to
  print no cells (same as *PRINT-ARRAY* set to NIL).  Set this to T
  to print all cells of the tensor.")

(defun print-tensor (tensor stream)
  (let ((dimensions (dimensions tensor))
	(offsets (offsets tensor))
	(store (store tensor)))
    (format stream " ~A~%" dimensions)
    (labels ((print-slice (initial-offset index)
	       (cond ((= index -1) (format stream "~A" (aref store initial-offset)))
		     (t (loop for k below (let ((s (aref dimensions index)))
					    (if (numberp *print-tensor*)
						(min s (1+ *print-tensor*))
						s))
			      for offset from initial-offset by (aref offsets index) do
			      (if (eql *print-tensor* k)
				  (format stream "...")
				  (print-slice offset (1- index)))
			      (if (zerop index)
				  (format stream "  ")
				  (terpri stream)))))))
      (print-slice (offset0 tensor) (1- (rank tensor))))))

(defmethod print-object ((tensor full-tensor) stream)
  (print-unreadable-object (tensor stream :type t :identity (not *print-tensor*))
    (when *print-tensor*
      (print-tensor tensor stream))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Tensor iteration
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; very slow
(defmethod for-each-entry-and-key (func (tensor tensor))
  (dotuple (index (coerce (dimensions tensor) 'list))
    (apply func (apply #'tensor-ref tensor index) index)))

(defmethod for-each-entry-and-vector-index (func (tensor tensor))
  (for-each-entry-and-key
   (lambda (entry &rest indices) (funcall func entry indices))
   tensor))

;;; fast
(definline for-each-index (func dims)
  (declare (type fixnum-vec dims))
  (loop with rank of-type fixnum = (length dims)
	with index of-type fixnum-vec = (make-fixnum-vec rank)
	do (funcall func index)
	while
	(dotimes (i rank nil)
	  (incf (aref index i))
	  (when (< (aref index i) (aref dims i)) (return t))
	  (setf (aref index i) 0))))

(definline for-each-offset-pair (func dims initial-off1 offsets1 initial-off2 offsets2)
  (declare (type fixnum-vec dims offsets1 offsets2))
  (declare (type fixnum initial-off1 initial-off2))
  (loop with off1 of-type fixnum = initial-off1
	and off2 of-type fixnum = initial-off2
	with rank of-type fixnum = (length dims)
	with index of-type fixnum-vec = (make-fixnum-vec rank)
	do
	(funcall func off1 off2)
	while
	(dotimes (i rank nil)
	  (incf (aref index i))
	  (incf off1 (aref offsets1 i))
	  (incf off2 (aref offsets2 i))
	  (when (< (aref index i) (aref dims i)) (return t))
	  (setf (aref index i) 0)
	  (decf off1 (the fixnum (* (aref offsets1 i) (aref dims i))))
	  (decf off2 (the fixnum (* (aref offsets2 i) (aref dims i)))))))

(defun get-remaining-inds (rank cinds)
  (loop with result = (make-fixnum-vec (- rank (length cinds)))
	and k = -1
	for i below rank
	unless (find i cinds) do
	(setf (aref result (incf k)) i)
	finally (return result)))

(defmethod slice ((tensor full-tensor) index-settings)
  (let ((dimensions (dimensions tensor))
	(offsets (offsets tensor)))
    (let ((inds (mapcar #'car index-settings))
	  (vals (mapcar #'cdr index-settings)))
      (assert (set-p inds))
      (let ((rinds (get-remaining-inds (rank tensor) inds)))
	(make-instance
	 (full-tensor (element-type tensor))
	 :dimensions (map 'fixnum-vec (curry #'aref dimensions) rinds)
	 :offset0 (reduce #'+ (map 'fixnum-vec #'(lambda (ind val) (* (aref offsets ind) val)) inds vals)
			  :initial-value (offset0 tensor))
	 :offsets (map 'fixnum-vec (curry #'aref offsets) rinds)
	 :store (store tensor))))))

(defmethod slice-copy ((tensor full-tensor) index-settings)
  (let ((dimensions (dimensions tensor))
	(offsets (offsets tensor))
	(store (store tensor)))
    (let ((inds (mapcar #'car index-settings)))
      (assert (set-p inds))
      (let* ((vals (mapcar #'cdr index-settings))
	     (rinds (get-remaining-inds (rank tensor) inds))
	     (rdims (map 'fixnum-vec (curry #'aref dimensions) rinds))
	     (result (make-instance (full-tensor (element-type tensor))
				    :dimensions rdims))
	     (rstore (store result)))
	(for-each-offset-pair
	 #'(lambda (off1 off2) (setf (aref rstore off1) (aref store off2)))
	 rdims
	 0 (offsets result)
	 (reduce #'+ (map 'fixnum-vec #'(lambda (ind val) (* (aref offsets ind) val)) inds vals)
		 :initial-value (offset0 tensor))
	 (map 'fixnum-vec (curry #'aref offsets) rinds))
	result))))

(defmethod copy ((tensor full-tensor))
  (slice-copy tensor ()))

(defmethod copy! ((a full-tensor) (b full-tensor))
  (assert (equalp (dimensions a) (dimensions b)))
  (let ((store-a (store a))
	(store-b (store b)))
    (for-each-offset-pair
     #'(lambda (off1 off2) (setf (aref store-b off2) (aref store-a off1)))
     (dimensions a) (offset0 a) (offsets a) (offset0 b) (offsets b)))
  b)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Transposition
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun rearrange-tensor (tensor permutation)
  (setq permutation (coerce permutation 'fixnum-vec))
  (assert (permutation-p permutation))
  (let ((dimensions (dimensions tensor))
	(offsets (offsets tensor))
	(offset0 (offset0 tensor))
	(store (store tensor)))
    (declare (type fixnum-vec dimensions offsets))
    (let* ((rdims (map 'fixnum-vec (curry #'aref dimensions) permutation))
	   (result (make-instance (full-tensor (element-type tensor))
				  :dimensions rdims))
	   (new-index (make-tensor-index result)))
      (let ((roff0 (offset0 result))
	    (roffs (offsets result))
	    (rstore (store result)))
	(declare (type fixnum offset0))
	(declare (type fixnum-vec roffs))
	(for-each-index
	 #'(lambda (index)
	     (declare (type fixnum-vec index))
	     (setf (aref rstore (compute-offset index offset0 offsets))
		   (aref store (compute-offset (permute-into permutation index new-index)
					       roff0 roffs))))
	 dimensions))
      result)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Tensor multiplication
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun check-cpairs (dims1 dims2 cpairs)
  "Checks that indices don't appear twice, and that the dimensions of every
contracted index pair fit."
  (labels ((check-rest (cpairs)
	     (or (null cpairs)
		 (let* ((cpair (car cpairs)) (i (car cpair)) (j (cdr cpair)) (rest (cdr cpairs)))
		   (and (= (aref dims1 i) (aref dims2 j))
			(not (assoc i rest)) (not (rassoc j rest))
			(check-rest rest))))))
    (check-rest cpairs)))

(defun consecutive-p (offsets dims)
  (loop for off across offsets
	and dim across dims
	and accumulated-off = 1 then (* accumulated-off dim)
	unless (= off accumulated-off) do (return nil)
	finally (return t)))

(defmethod t* ((tensor1 full-tensor) (tensor2 full-tensor) contraction-pairs)
  (assert (and (typep tensor1 (full-tensor 'double-float))
	       (typep tensor2 (full-tensor 'double-float))))
  (let ((dims1 (dimensions tensor1))
	(off0-1 (offset0 tensor1))
	(offsets1 (offsets tensor1)))
    (let ((dims2 (dimensions tensor2))
	  (off0-2 (offset0 tensor2))
	  (offsets2 (offsets tensor2)))
      (assert (check-cpairs dims1 dims2 contraction-pairs))
      (let* ((rank1 (length dims1))
	     (rank2 (length dims2))
	     (cinds1 (map 'fixnum-vec #'car contraction-pairs))
	     (cinds2 (map 'fixnum-vec #'cdr contraction-pairs))
	     (rinds1 (get-remaining-inds rank1 cinds1))
	     (rinds2 (get-remaining-inds rank2 cinds2))
	     (perm1 (concatenate 'fixnum-vec cinds1 rinds1))
	     (perm2 (concatenate 'fixnum-vec cinds2 rinds2))
	     (shuffled1 (if (and (identity-permutation-p perm1)
				 (zerop off0-1)
				 (consecutive-p offsets1 dims1))
			    tensor1
			    (time (rearrange-tensor tensor1 perm1))))
	     (shuffled2 (if (and (identity-permutation-p perm2)
				 (zerop off0-2)
				 (consecutive-p offsets2 dims2))
			    tensor2 (rearrange-tensor tensor2 perm2)))
	     (store1 (store shuffled1))
	     (store2 (store shuffled2))
	     (rank-part1 (length rinds1))
	     (rank-part2 (length rinds2))
	     (rank (+ rank-part1 rank-part2)))
	(cond
	  ((zerop rank)
	   (let ((x store1)
		 (y store2))
	     (declare (type double-vec x y))
	     (declare (optimize (speed 3) (safety 0)))
	     (let* ((result (make-real-tensor nil))
		    (store (store result)))
	       (declare (type double-vec store))
	       (setf (aref store 0)
		     (let ((sum 0.0))
		       (declare (type double-float sum))
		       (dotimes (i (length x) sum)
			 (incf sum (* (aref x i) (aref y i))))))
	       result)))
	  (t
	   (let* ((crank (length cinds1))
		  (offsets1 (offsets shuffled1))
		  (offsets2 (offsets shuffled2))
		  (rdims1 (vector-map (curry #'aref dims1) rinds1))
		  (rdims2 (vector-map (curry #'aref dims2) rinds2)) 
		  (roffs1 (map 'fixnum-vec (curry #'aref offsets1) (range< crank rank1)))
		  (roffs2 (map 'fixnum-vec (curry #'aref offsets2) (range< crank rank2)))
		  (contraction-length (aref roffs1 0))
		  (dims (concatenate 'fixnum-vec rdims1 rdims2))
		  (result (make-real-tensor dims))
		  (offsets (offsets result))
		  (store (store result)))
	     (declare (type fixnum contraction-length))
	     (loop
	      ;; now the indices from 0,...,crank-1 will be contracted
	      with index1 = (make-fixnum-vec rank-part1)
	      with off1 of-type fixnum = 0
	      with offset of-type fixnum = 0 do
	      (loop with index2 = (make-fixnum-vec rank-part2)
		    with off2 of-type fixnum = 0 do
		    (setf (aref store offset)
			  (let ((x store1)
				(y store2))
			    (declare (type double-vec x y))
			    (declare (optimize (speed 3) (safety 0)))
			    (loop for i from off1 below (the fixnum (+ off1 contraction-length))
				  and j from off2
				  summing (* (aref x i) (aref y j)) double-float)			    
			    ))
		    while
		    (dotimes (i2 rank-part2 nil)
		      (incf (aref index2 i2))
		      (incf off2 (aref roffs2 i2))
		      (incf offset (aref offsets (+ i2 rank-part1)))
		      (when (< (aref index2 i2) (aref rdims2 i2)) (return t))
		      (setf (aref index2 i2) 0)
		      (decf off2 (* (aref roffs2 i2) (aref rdims2 i2)))
		      (decf offset (* (aref offsets (+ i2 rank-part1)) (aref dims (+ i2 rank-part1))))))
	      while
	      (dotimes (i1 rank-part1 nil)
		(incf (aref index1 i1))
		(incf off1 (aref roffs1 i1))
		(incf offset (aref offsets i1))
		(when (< (aref index1 i1) (aref rdims1 i1)) (return t))
		(setf (aref index1 i1) 0)
		(decf off1 (* (aref roffs1 i1) (aref rdims1 i1)))
		(decf offset (* (aref offsets i1) (aref dims i1))))
	      finally (return result)))))
	     ))))

;;; Testing
(defun test-tensor ()
  (let ((x (vector (eye 1))))
    (tensor-ref x 0 0)
    (setf (tensor-ref x 0 0) 2.0)
    (tensor-ref x 0 0))

  (let ((t1 (list->real-tensor '((1.0 2.0) (3.0 4.0))))
	(t2 (list->real-tensor '((1.0 2.0) (3.0 4.0)))))
    (make-analog t1)
    (slice t1 '((0 . 1)))
    (m+ t1 t2)
    (rearrange-tensor t1 #(1 0))
    (t* t1 t2 '((0 . 1))))
  
  (let ((t1 (list->real-tensor '(1.0 2.0)))
	(t2 (list->real-tensor '(1.0 2.0))))
    (t* t1 t2 '((0 . 0)))
    (t* t1 t2 '())
    (tensor-map #'1+ t1)
    (copy! t1 t2))
  
  ;; performance test
  (time
   (let* ((n 50)
	  (t1 (make-real-tensor (make-fixnum-vec 2 n)))
	  (t2 (make-real-tensor (make-fixnum-vec 2 n))))
     (fill! t1 1.0)
     (fill! t2 2.0)
     (loop repeat 10 do
	  (t* t1 t2 '((1 . 0))))))

  #+(or)
  (time
   (let* ((n 50)
	  (t1 (make-real-matrix n))
	  (t2 (make-real-matrix n)))
     (fill! t1 1.0)
     (fill! t2 2.0)
     (loop repeat 10 do (m* t1 t2))))
  )

;;; (test-tensor)
(fl.tests:adjoin-test 'test-tensor)



