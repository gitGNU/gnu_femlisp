;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; sparse-vector.lisp
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;
;;; Copyright (C) 2003-2006 Nicolas Neuss, University of Heidelberg.
;;; Copyright (C) 2007- Nicolas Neuss, University of Karlsruhe.
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

;;;; This module defines a graph sparse matrix format.  It is somewhat
;;;; similar to the format used in UG and outlined in [Neuss1998], but it
;;;; uses hash-table sparse matrices which allow for O(1) random access and
;;;; indexing over arbitrary key sets.

;;; (declaim (optimize (safety 3) (debug 3)))

(in-package :fl.matlisp)

(defvar *parallel-algebra* nil
  "Preliminary switch for allowing parallel linear algebra.  Will hopefully
  become unnecessary in the long run.")

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; <sparse-vector>
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass <sparse-vector> (<vector>)
  ((multiplicity :reader multiplicity :initform 1 :initarg :multiplicity
		 :documentation "Multiplicity of the sparse vector.  A
multiplicity different from 1 is used when handling multiple right-hand sides
and solutions simultaneously."))
  (:documentation "Abstract class for sparse vectors."))

(defmethod element-type ((svec <sparse-vector>))
  (standard-matrix 'double-float))

(defmethod scalar-type ((svec <sparse-vector>))
  'double-float)

(defgeneric key->size (svec)
  (:documentation "Returns NIL, if @arg{svec} cannot extend automatically when
  being accessed.  Otherwise returns a function mapping keys to vector block
  sizes.")
  (:method ((svec <sparse-vector>)) nil))

(defclass <key->size-mixin> ()
  ((key->size :reader key->size :initarg :key->size :type function
	      :documentation "Function determining the dimension of a block."))
  (:documentation "A mixin which is not used up to now, because the same
  functionality is provided using information from the ansatz-space."))

;;; also fl.parallel::locked-region-mixin can be included here

(defgeneric vector-block (svec key)
  (:documentation "Low-level key lookup.  Returns NIL if there is no block at
  this position."))

(defgeneric (setf vector-block) (value svec key)
  (:documentation "Low-level insertion of a block without checks."))

(defmethod entry-allowed-p ((svec <sparse-vector>) &rest indices)
  (destructuring-bind (key) indices
    (or (vector-block svec key)
        (null (key->size svec))
        (plusp (funcall (key->size svec) key)))))

(defmethod vref ((svec <sparse-vector>) key)
  (or (vector-block svec key)
      (setf (vector-block svec key)
            (let ((ks (key->size svec)))
              (cond
                ((null ks) (error "Vector does not extend automatically"))
                ((not (plusp (funcall ks key)))
                 (error "Request for an empty matrix"))
                (t (make-real-matrix (funcall ks key)
                                     (multiplicity svec))))))))

(defmethod (setf vref) (value (svec <sparse-vector>) key)
  (setf (vector-block svec key) value))

(defmethod make-analog ((svec <key->size-mixin>))
  (copy-slots (call-next-method) svec '(key->size)))

(defmethod show ((svec <sparse-vector>) &key keys (zeros t) &allow-other-keys)
  (format t "~&Sparse vector with ~D allocated components:~%"
	  (nr-of-entries svec))
  (dolist (key (or keys (keys svec)))
    (let ((vblock (vref svec key)))
      (when (or zeros (not (mzerop vblock)))
        (format t "~A --> ~A~%" key vblock)))))

(defmethod ncols ((svec <sparse-vector>))
  (multiplicity svec))

(defgeneric remove-key (sobj &rest indices)
  (:documentation "Remove the entry for @arg{key} from the sparse object."))

(defun remove-keys (sobj keys)
  (loop for key in keys do (remove-key sobj key)))

(defun sparse-vector->matlisp (svec &optional keys ranges)
  "Transforms all or a part of @arg{svec} corresponding to the keys in
@arg{keys} and maybe the ranges in 'ranges' to a matlisp matrix."
  (setq keys (coerce (or keys (keys svec)) 'vector))
  (setq ranges (and ranges (coerce ranges 'vector)))
  (lret* ((n (if ranges
                 (reduce #'+ ranges :key #'(lambda (x) (- (cdr x) (car x))))
                 (reduce #'+ keys :key (key->size svec))))
          (multiplicity (multiplicity svec))
          (mm (make-real-matrix n multiplicity)))
    (loop for key across keys and k from 0
	  and offset of-type fixnum = 0 then (+ offset (- end-comp start-comp))
	  for start-comp of-type fixnum = (if ranges (car (aref ranges k)) 0)
	  for end-comp of-type fixnum = (if ranges
					    (cdr (aref ranges k))
					    (funcall (key->size svec) key))
	  for entry = (vref svec key) do
            (when entry
              (extended-minject! entry mm offset 0 start-comp 0 end-comp multiplicity)))))

(defgeneric extract-value-blocks (sobj keys &optional col-keys)
  (:documentation "Extract a vector or array of value blocks from @arg{sobj}."))

(defmethod extract-value-blocks ((svec <sparse-vector>) keys &optional col-keys)
  (assert (null col-keys))
  (map 'vector
       (lambda (key)
         (and (entry-allowed-p svec key)
              (vref svec key)))
       keys))

(defun combine-svec-block (svec local-vec keys ranges operation)
  "Puts a local block in matlisp format into a <sparse-vector>."
  (setq keys (coerce (or keys (keys svec)) 'vector))
  (setq ranges (and ranges (coerce ranges 'vector)))
  (let ((multiplicity (multiplicity svec)))
    (loop for key across keys and k from 0
	  and offset of-type fixnum = 0 then (+ offset (- end-comp start-comp))
	  for start-comp of-type fixnum = (if ranges (car (aref ranges k)) 0)
	  for end-comp of-type fixnum = (if ranges
					    (cdr (aref ranges k))
					    (funcall (key->size svec) key))
	  for entry = (vref svec key) do
	  (loop for i of-type fixnum from start-comp below end-comp do
		(dotimes (j multiplicity)
		  (let ((local-entry (mref local-vec (+ offset i) j)))
		    (ecase operation
		      (:add (incf (mref entry i j) local-entry))
		      (:set (setf (mref entry i j) local-entry)))))))))
  
(defgeneric set-svec-to-local-block (svec local-vec &optional keys ranges)
  (:documentation "Copies a local block in matlisp format into a <sparse-vector>.")
  (:method ((svec <sparse-vector>) local-vec &optional keys ranges)
      (combine-svec-block svec local-vec keys ranges :set)))

(defgeneric add-svec-to-local-block (svec local-vec &optional keys ranges)
  (:documentation "Copies a local block in matlisp format into a <ht-sparse-vector>.")
  (:method ((svec <sparse-vector>) local-vec &optional keys ranges)
      (combine-svec-block svec local-vec keys ranges :add)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; matlisp operations for the <sparse-vector> class
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defgeneric extract-if (test svec &key &allow-other-keys)
  (:documentation "Extract a subvector or submatrix from a sparse
  vector/matrix.")
  (:method ((test function) (svec <vector>) &key &allow-other-keys)
    "Default method on vectors."
    (lret ((sub-vec (make-analog svec)))
      (dovec ((key entry) svec)
        (when (funcall test key entry)
          (setf (vref sub-vec key) entry))))))

(defmethod minject! ((x <sparse-vector>) (y <sparse-vector>) row-off col-off)
  "This routine can only inject in the multiplicity dimension."
  (assert (or (null row-off) (zerop row-off)))
  (dovec ((xc i) x)
    (minject! xc (vref y i) 0 col-off)))

(defmethod mextract! ((x <sparse-vector>) (y <sparse-vector>) row-off col-off)
  "This routine can only extract from the multiplicity dimension."
  (assert (or (null row-off) (zerop row-off)))
  (dovec ((yc i) y)
    (mextract! (vref x i) yc 0 col-off)))

(defun make-analog-with-multiplicity (x multiplicity)
  (lret ((result (make-analog x)))
    (setf (slot-value result 'multiplicity) multiplicity)))

(defmethod matrix-slice ((x <sparse-vector>) &key
			 from-row (from-col 0) nrows
			 (ncols (- (ncols x) from-col)))
  (assert (and (null from-row) (null nrows)))
  (lret ((result (make-analog-with-multiplicity x ncols)))
    (mextract! result x from-row from-col)))

(defmethod join-instance (orientation (x <sparse-vector>) &rest vecs)
  (unless (eq orientation :horizontal)
    (error "Only the horizontal join of sparse vectors is allowed."))
  (make-analog-with-multiplicity
   x (reduce #'+ (cons x vecs) :key #'multiplicity)))

(defmethod join-horizontal! ((result <sparse-vector>) &rest vectors)
  (loop with k = 0
     for vec in vectors do
     (minject! vec result 0 k)
     (incf k (multiplicity vec))
     finally (assert (= k (multiplicity result))))
  result)

(defmethod m*-tn-product-instance ((x <sparse-vector>) (y <sparse-vector>))
  (make-instance (standard-matrix (scalar-type x))
		 :nrows (ncols x) :ncols (ncols y)))

(defmethod gemm-tn! (alpha (x <sparse-vector>) (y <sparse-vector>) beta (z standard-matrix))
  (scal! beta z)
  (dovec ((xc i) x)
    (gemm-tn! alpha xc (vref y i) 1.0 z))
  z)

(defmethod copy! ((x <sparse-vector>) (y <sparse-vector>))
  ;; because the entries of x might not yet exist in y,
  ;; they have to be generated in parallel.  This is problematic
  ;; for all containers we use at the moment.
  (let ((*parallel-algebra* nil))
    (call-next-method)))

(defmethod m* ((x <sparse-vector>) (y standard-matrix))
  (lret ((result (make-analog-with-multiplicity x (ncols y))))
    (dovec ((xc i) x)
      (setf (vref result i) (m* xc y)))))

(defmethod dot ((x <sparse-vector>) (y <sparse-vector>))
  (with-accumulators (sum 0 #'+)
    (for-each-entry-and-key
     (lambda (xc i) (incf sum (dot xc (vref y i))))
     x)))

(defmethod norm ((x <sparse-vector>) &optional (p 2))
  (case p
    (:inf
     (with-accumulators (sum 0 #'max)
       (for-each-entry
        (lambda (xc)
          (setf sum (max sum (norm xc :inf))))
        x)))
    (2 (sqrt (dot x x)))
    (t (expt (with-accumulators (sum 0 #'+)
               (for-each-entry
                (lambda (xc)
                  (incf sum (expt (norm xc p) p)))
                x))
             (/ p)))))


(defmethod print-object :after ((svec <sparse-vector>) stream)
  (format stream "{nrows=~A, mult=~A}"
          (nr-of-entries svec) (multiplicity svec)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; <sparse-dictionary-vector>
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass <sparse-dictionary-vector> (<sparse-vector>)
  ()
  (:documentation "Sparse block vector class indexed with general keys."))

(defgeneric blocks (vec)
  (:documentation "Returns a dictionary mapping keys to entries for @arg{vec}."))

(defmethod vector-block ((svec <sparse-dictionary-vector>) key)
  (dic-ref (blocks svec) key))

(defmethod (setf vector-block) (value (svec <sparse-dictionary-vector>) key)
  (setf (dic-ref (slot-value svec 'blocks) key) value))

(defmethod for-each-key (fn (svec <sparse-dictionary-vector>))
  (dic-for-each-key fn (blocks svec) :parallel *parallel-algebra*))

(defmethod for-each-entry-and-key (fn (svec <sparse-dictionary-vector>))
  (dic-for-each (lambda (key value) (funcall fn value key))
                (blocks svec) :parallel *parallel-algebra*))

(defmethod for-each-entry (fn (svec <sparse-dictionary-vector>))
  (dic-for-each-value fn (blocks svec) :parallel *parallel-algebra*))

(defmethod remove-key ((svec <sparse-dictionary-vector>) &rest indices)
  (destructuring-bind (key) indices
    (dic-remove key (blocks svec))))

#+(or)  ; Wrong!  We want a deep copy!
(defmethod copy ((svec <sparse-dictionary-vector>))
  (lret ((result (make-analog svec)))
    (for-each-entry-and-key (lambda (value key)
                              (setf (vref result key) value))
                            svec)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; <ht-sparse-vector>
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass <ht-sparse-vector> (locked-region-mixin <sparse-dictionary-vector>)
  ((test :initform 'eql :initarg :test)
   (blocks :reader blocks :type hash-table
	   :documentation "Table of blocks."))
  (:documentation "Sparse block vector class implemented using a hash-table."))

(defmethod initialize-instance :after ((svec <ht-sparse-vector>) &key &allow-other-keys)
  (with-slots (test blocks fl.parallel::locked-region) svec
    (setf blocks (make-hash-table :test test)
          fl.parallel::locked-region (make-hash-table :test test))))

(defun make-sparse-vector
    (&rest args
     &key (type '(<key->size-mixin> <ht-sparse-vector>))
     &allow-other-keys)
  (apply #'fl.amop::make-programmatic-instance type
         (sans args '(:type))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Tests
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun test-sparse-vector ()
  (flet ((constantly-1 (x) (declare (ignore x)) 1))
    (let* ((x (make-sparse-vector :key->size #'constantly-1))
	   (y (make-sparse-vector :key->size #'constantly-1))
	   (z (make-sparse-vector :key->size #'constantly-1 :multiplicity 2)))
      (setf (vref x 1) #m((1.0)))
      (setf (vref y 1) #m((1.0)))
      (vref x 2)
      ;; a parallelization test
      (time
       (with-workers ((lambda (i)
                        (with-region (x (list i))
                          (loop repeat 10 do
                               (scal! 2.0 (vref x i))
                               (sleep 0.1)))))
         (work-on 1)
         (work-on 1)
         (work-on 2)
         ))
      (copy! y x)
      (assert (= (mref (vref x 1) 0 0) (mref (vref y 1) 0 0)))
      (x<-0 x)
      (assert (mzerop (vref x 1)))
      (axpy! 2.0 y x)
      (assert (= (norm x 1) 2.0))
      (minject! x z nil 1)
      (show z)
      (m*-tn z z)
      (mextract! y z nil 1)
      (show x)
      (show z)
      (show (join :horizontal x z))
      ))
  )

;;; (test-sparse-vector)
(fl.tests:adjoin-test 'test-sparse-vector)

