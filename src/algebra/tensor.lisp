;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; tensor.lisp
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

;;;; Provides an extension of vectors and (matlisp) matrices operations to
;;;; tensors of arbitrary dimension.

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Class definition and basic methods
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass <tensor> ()
  ((dimensions :reader dimensions :initarg :dimensions :type fixnum-vec)
   (offset0 :reader offset0 :initform 0 :initarg :offset0 :type fixnum)
   (offsets :reader offsets :initarg :offsets :type fixnum-vec)))

(defclass <general-tensor> (<tensor>)
  ((entries :reader entries :initarg :entries :type (array t (*)))))
(defclass <real-tensor> (<tensor>)
  ((entries :reader entries :initarg :entries :type double-vec)))
(defclass <complex-tensor> (<tensor>)
  ((entries :reader entries :initarg :entries :type (array complex (*)))))

(defun dimensions->offsets (dims)
  (loop	with offsets = (make-array (length dims) :element-type 'fixnum)
	for k from 0
	and dim across dims
	and offset = 1 then (* offset dim) do
	(setf (aref offsets k) offset)
	finally
	(return offsets)))

(definline tensor-class-for (type)
  (ecase type
    ((t) '<general-tensor>)
    ((double-float) '<real-tensor>)
    ((complex) '<complex-tensor>)))

(definline tensor-element-type (tensor)
  (ecase (type-of tensor)
    ((<general-tensor>) t)
    ((<real-tensor>) 'double-float)
    ((<complex-tensor>) 'complex)))

(defun make-tensor (dimensions type)
  (let* ((dims (coerce dimensions 'fixnum-vec))
	 (offsets (dimensions->offsets dims))
	 (total (if (zerop (length dims))
		    1
		    (* (vector-last dims) (vector-last offsets)))))
    (make-instance (tensor-class-for type)
	     :dimensions dims :offsets offsets
	     :entries (make-array total :element-type type))))

(defun make-general-tensor (dimensions)
  (make-tensor dimensions t))

(defun make-real-tensor (dimensions)
  (make-tensor dimensions 'double-float))

(defun make-complex-tensor (dimensions)
  (make-tensor dimensions 'complex))

(defun list->real-tensor (tree)
  (assert (tree-uniformp tree))
  (let ((tensor (make-real-tensor (coerce (tree-uniform-number-of-branches tree) 'fixnum-vec)))
	(k -1))
    (on-leaves #'(lambda (leaf) (setf (aref (entries tensor) (incf k)) leaf))
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

(defmethod tensor-ref ((tensor <tensor>) &rest indices)
  (aref (entries tensor) (compute-offset (coerce indices 'fixnum-vec)
				       (offset0 tensor) (offsets tensor))))

(defmethod (setf tensor-ref) (value (tensor <tensor>) &rest indices)
  (setf (aref (entries tensor)
	      (compute-offset (coerce indices 'fixnum-vec)
			      (offset0 tensor) (offsets tensor)))
	value))

(defmethod rank ((tensor <tensor>)) (length (dimensions tensor)))
(defun make-tensor-index (tensor) (make-fixnum-vec (rank tensor)))

(defparameter *print-tensor*
  5
  "Maximum number of columns and/or rows to print.  Set this to NIL to
  print no cells (same as *PRINT-ARRAY* set to NIL).  Set this to T
  to print all cells of the tensor.")

(defun print-tensor (tensor stream)
  (let ((dimensions (dimensions tensor))
	(offsets (offsets tensor))
	(entries (entries tensor)))
    (format stream " ~A~%" dimensions)
    (labels ((print-slice (initial-offset index)
	       (cond ((= index -1) (format stream "~A" (aref entries initial-offset)))
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

(defmethod print-object ((tensor <tensor>) stream)
  (print-unreadable-object (tensor stream :type t :identity (not *print-tensor*))
    (when *print-tensor*
      (print-tensor tensor stream))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Tensor iteration
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

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

(defmethod slice ((tensor <tensor>) index-settings)
  (let ((dimensions (dimensions tensor))
	(offsets (offsets tensor)))
    (let ((inds (mapcar #'car index-settings))
	  (vals (mapcar #'cdr index-settings)))
      (assert (set-p inds))
      (let ((rinds (get-remaining-inds (rank tensor) inds)))
	(make-instance
	 (tensor-class-for (tensor-element-type tensor))
	 :dimensions (map 'fixnum-vec (curry #'aref dimensions) rinds)
	 :offset0 (reduce #'+ (map 'fixnum-vec #'(lambda (ind val) (* (aref offsets ind) val)) inds vals)
			  :initial-value (offset0 tensor))
	 :offsets (map 'fixnum-vec (curry #'aref offsets) rinds)
	 :entries (entries tensor))))))

(defmethod slice-copy ((tensor <tensor>) index-settings)
  (let ((dimensions (dimensions tensor))
	(offsets (offsets tensor))
	(entries (entries tensor)))
    (let ((inds (mapcar #'car index-settings)))
      (assert (set-p inds))
      (let* ((vals (mapcar #'cdr index-settings))
	     (rinds (get-remaining-inds (rank tensor) inds))
	     (rdims (map 'fixnum-vec (curry #'aref dimensions) rinds))
	     (result (make-tensor rdims (tensor-element-type tensor)))
	     (rentries (entries result)))
	(for-each-offset-pair
	 #'(lambda (off1 off2) (setf (aref rentries off1) (aref entries off2)))
	 rdims
	 0 (offsets result)
	 (reduce #'+ (map 'fixnum-vec #'(lambda (ind val) (* (aref offsets ind) val)) inds vals)
		 :initial-value (offset0 tensor))
	 (map 'fixnum-vec (curry #'aref offsets) rinds))
	result))))

(defmethod copy ((tensor <tensor>))
  (slice-copy tensor ()))

(defmethod copy! ((a <tensor>) (b <tensor>))
  (assert (equalp (dimensions a) (dimensions b)))
  (let ((entries-a (entries a))
	(entries-b (entries b)))
    (for-each-offset-pair
     #'(lambda (off1 off2) (setf (aref entries-b off2) (aref entries-a off1)))
     (dimensions a) (offset0 a) (offsets a) (offset0 b) (offsets b)))
  b)

(defmethod tensor-map ((type symbol) (func function) (tensor <tensor>))
  (let* ((result (make-tensor (dimensions tensor) type))
	 (entries (entries tensor))
	 (rentries (entries result)))
    (for-each-offset-pair
     #'(lambda (off roff)
	 (setf (aref rentries roff)
	       (funcall func (aref entries off))))
     (dimensions tensor) (offset0 tensor) (offsets tensor) (offset0 result) (offsets result))
    result))
     
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Tensor addition
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmethod t+ (a b)
  (assert (equalp (dimensions a) (dimensions b)))
  (let* ((result (copy b))
	 (entries1 (entries result))
	 (entries2 (entries a)))
    (for-each-offset-pair
     #'(lambda (off1 off2) (incf (aref entries1 off1) (aref entries2 off2)))
     (dimensions a) (offset0 result) (offsets result) (offset0 a) (offsets a))
    result))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Transposition
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun rearrange-tensor (tensor permutation)
  (declare (type <real-tensor> tensor))
  (setq permutation (coerce permutation 'fixnum-vec))
  (assert (permutation-p permutation))
  (let ((dimensions (dimensions tensor))
	(offsets (offsets tensor))
	(offset0 (offset0 tensor))
	(entries (entries tensor)))
    (declare (type fixnum-vec dimensions offsets))
    (declare (type double-vec entries))
    (let* ((rdims (map 'fixnum-vec (curry #'aref dimensions) permutation))
	   (result (make-tensor rdims (tensor-element-type tensor)))
	   (new-index (make-tensor-index result)))
      (declare (type <real-tensor> result))
      (let ((roff0 (offset0 result))
	    (roffs (offsets result))
	    (rentries (entries result)))
	(declare (type fixnum offset0))
	(declare (type fixnum-vec roffs))
	(declare (type double-vec rentries))
	(declare (optimize (speed 3) (safety 0)))
	(for-each-index
	 #'(lambda (index)
	     (declare (type fixnum-vec index))
	     (setf (aref rentries (compute-offset index offset0 offsets))
		   (aref entries (compute-offset (permute-into permutation index new-index)
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

(defun get-remaining-inds (rank cinds)
  (loop with result = (make-fixnum-vec (- rank (length cinds)))
	and k = -1
	for i below rank
	unless (find i cinds) do
	(setf (aref result (incf k)) i)
	finally (return result)))

(defun consecutive-p (offsets dims)
  (loop for off across offsets
	and dim across dims
	and accumulated-off = 1 then (* accumulated-off dim)
	unless (= off accumulated-off) do (return nil)
	finally (return t)))

(defmethod t* ((tensor1 <real-tensor>) (tensor2 <real-tensor>) contraction-pairs)
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
			    tensor2
			    (time (rearrange-tensor tensor2 perm2))))
	     (entries1 (entries shuffled1))
	     (entries2 (entries shuffled2))
	     (rank-part1 (length rinds1))
	     (rank-part2 (length rinds2))
	     (rank (+ rank-part1 rank-part2)))
	(cond
	  ((zerop rank)
	   (let ((x entries1)
		 (y entries2))
	     (declare (type double-vec x y))
	     (declare (optimize (speed 3) (safety 0)))
	     (let* ((result (make-real-tensor nil))
		    (entries (entries result)))
	       (declare (type double-vec entries))
	       (setf (aref entries 0)
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
		  (entries (entries result)))
	     (declare (type fixnum contraction-length))
	     (loop
	      ;; now the indices from 0,...,crank-1 will be contracted
	      with index1 = (make-fixnum-vec rank-part1)
	      with off1 of-type fixnum = 0
	      with offset of-type fixnum = 0 do
	      (loop with index2 = (make-fixnum-vec rank-part2)
		    with off2 of-type fixnum = 0 do
		    (setf (aref entries offset)
			  (let ((x entries1)
				(y entries2))
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
  (let ((t1 (list->real-tensor '((1.0 2.0) (3.0 4.0))))
	(t2 (list->real-tensor '((1.0 2.0) (3.0 4.0)))))
    (slice t1 '((0 . 1)))
    (t+ t1 t2)
    (rearrange-tensor t1 #(1 0))
    (t* t1 t2 '((0 . 1))))
  
  (let ((t1 (list->real-tensor '(1.0 2.0)))
	(t2 (list->real-tensor '(1.0 2.0))))
    (t* t1 t2 '((0 . 0)))
    (t* t1 t2 '())
    (tensor-map 'double-float #'1+ t1)
    (copy! t1 t2))
  
  )

;;; (test-tensor)
(fl.tests:adjoin-test 'test-tensor)

;;; Finally, a performance test: The new routine t* implements probably the
;;; best way for tensor multiplication by reshuffling the tensors before
;;; multiplication such that the contraction is over successive floats.
;;; Then we are not much slower than Matlisp/LAPACK for large matrices.  An
;;; old version of t* can be found in the CVS tree, which is several times
;;; slower.  Nevertheless, further enhancements might be possible, when
;;; subdividing the matrices in blocks which remain in the first-level
;;; cache.  For comparison, see the values obtained by mflop.lisp.

#+ignore
(time
 (let* ((n 50)
	(t1 (make-real-tensor (make-fixnum-vec 2 n)))
	(t2 (make-real-tensor (make-fixnum-vec 2 n))))
   (fill! (entries t1) 1.0)
   (fill! (entries t2) 2.0)
   (loop repeat 1000 do
	 (t* t1 t2 '((1 . 0))))))
; 100: 0.12, 200: 0.87, 300: 2.52, 400: 5.5, 500: 10.1

#+ignore
(time
 (let* ((n 50)
	(t1 (make-real-matrix n))
	(t2 (make-real-matrix n)))
   (fill! t1 1.0)
   (fill! t2 2.0)
   (loop repeat 1000 do (m* t1 t2))))
; 100: 0.02, 200: 0.31, 300: 1.43, 400: 3.39, 500: 6.7

