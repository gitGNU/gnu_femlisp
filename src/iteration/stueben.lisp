;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; stueben.lisp - a variant of Stueben's AMG
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

(in-package :fl.multigrid)

;;;; This file contains coarsening and interpolation routines for an
;;;; algebraic multigrid of selection type.  Such AMG algorithms were first
;;;; proposed in 1983 by Ruge and Stueben, for an intermediate reference on
;;;; an often-used variant see their 1987 paper on AMG which appeared in
;;;; 'Multigrid methods' edited by S. F. McCormick.

(defclass <stueben> (<selection-amg>)
  ((theta :accessor theta :initform 0.25 :initarg :theta))
  (:documentation "This provides something like Stueben's variant for
selection amg.  The original Ruge-Stueben algorithm was developed further
since 1987 by Klaus Stueben and collaborators.  These developments are
published in ....  The algorithm implemented here uses their ideas, but
does not claim to be equivalent to their code which can be bought at SCAI,
St. Augustin, Germany.  At this point, I want to thank Tanja Clees
for several discussions on AMG."))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Ruge-Stueben prolongation
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; Filtering of strong connections - this may be of interest also for
;;; other AMG algorithms.

(defun min-in-row (mat row &optional min)
  (for-each-key-and-entry-in-row
   #'(lambda (col entry)
       (unless (eq col row)
	 (let ((current (mref entry 0 0)))
	   (when (or (eq min nil) (< current min))
	     (setq min current)))))
   mat row)
  (unless min (setq min 0.0))
  (when (plusp min) (setq min 0.0))
  ;;(assert (not (plusp min)))
  min)

(defmethod filtered-matrix ((amg <stueben>) blackboard)
  "Special filtering used by Ruge-Stueben."
  (let ((matrix (getbb blackboard :matrix))
	(theta (theta amg))
	(filtered-keys (make-hash-table))
	(pivot-table (make-hash-table)))
    
    (flet ((strong-p (row col entry)
	     (and (not (eq row col))
		  (< (mref entry 0 0) (gethash row pivot-table)))))
      
      ;; fill pivot-table
      (for-each-row-key
       #'(lambda (row-key)
	   (setf (gethash row-key pivot-table)
		 (if theta
		     (* theta (min-in-row matrix row-key))
		     (min-in-row matrix row-key 0.0))))
       matrix)
      
      ;; set filtered-keys
      (for-each-row-key
       #'(lambda (row-key)
	   (unless (zerop (gethash row-key pivot-table))
	     (setf (gethash row-key filtered-keys) t)))
       matrix)
      
      #+(or)
      (for-each-key-and-entry
       #'(lambda (row-key col-key entry)
	   (and (gethash row-key filtered-keys)
		(gethash col-key filtered-keys)
		(or (not theta) (strong-p row-key col-key entry))
		(setf (mref filtered-matrix row-key col-key) entry)))
       matrix)
      
      ;; return result
      (setf (getbb blackboard :filtered-keys) filtered-keys)
      (setf (getbb blackboard :filtered-matrix)
	    (extract-if
	      #'(lambda (row-key col-key entry)
		  (and (gethash row-key filtered-keys)
		       (gethash col-key filtered-keys)
		       (or (not theta) (strong-p row-key col-key entry))))
	      matrix))))
  blackboard)

;;; Priority-table - this is used for the Ruge-Stueben coarsening

(defclass priority-table ()
  ((table :reader pt-table :initform
	  (make-array 32 :initial-element () :adjustable t :fill-pointer 0))
   (dictionary :reader pt-dictionary :initform (make-hash-table)))
  (:documentation "Objects of this class are a dictionary for object/priority
pairs.  It allows to pop objects of maximum or minimum priority.  Objects can
also be deleted.

Internally, table is an adjustable vector with fill-pointer being either 0 or
larger than 0, if there is an entry in the doubly-linked list at the
corresponding index."))

(defmethod pt-insert ((pt priority-table) obj index)
  "Puts an object in the priority-table."
  (let ((table (pt-table pt))
	(dictionary (pt-dictionary pt)))
    (assert (not (gethash obj dictionary)))
    (loop until (< index (fill-pointer table)) do
	  (if (< index (array-dimension table 0))
	      (setf (fill-pointer table) (1+ index))
	      (adjust-array table (* 2 (array-dimension table 0)) :initial-element ())))
    (unless (aref table index)
      (setf (aref table index) (make-dll)))
    ;; we insert obj in the suitable slot of the table, the dll-item returned is
    ;; put in the translation table together with the index
    (setf (gethash obj dictionary)
	  (cons index
		(dll-front-insert obj (aref table index))))))

(defun adapt-fill-pointer (table)
  (loop
   while (and (plusp (fill-pointer table))
	      (let ((dll (aref table (1- (fill-pointer table)))))
		(or (null dll) (dll-empty-p dll))))
   do (decf (fill-pointer table))))

(defmethod pt-remove ((pt priority-table) obj)
  "Removes a node from the priority table."
  (let ((table (pt-table pt))
	(dictionary (pt-dictionary pt)))
    (destructuring-bind (index . item)
      (gethash obj dictionary)
    (dll-remove-item item (aref table index))
    (remhash obj dictionary)
    (when (= (fill-pointer table) (1+ index))
      (adapt-fill-pointer table)))))

(defmethod pt-shift-priority ((pt priority-table) obj shift)
  "Changes the priority of some node."
  (let ((index (car (gethash obj (pt-dictionary pt)))))
    (unless (zerop shift)
      (pt-remove pt obj)
      (pt-insert pt obj (+ index shift)))))

(defmethod pt-pop ((pt priority-table))
  "Gets the node with the highest priority from the priority table."
  (let ((table (pt-table pt))
	(dictionary (pt-dictionary pt)))
    (and (plusp (fill-pointer table))
	 (let ((obj (dll-pop-first (aref table (1- (fill-pointer table))))))
	   (remhash obj dictionary)
	   (adapt-fill-pointer table)
	   obj))))

(defmethod pt-in-table-p ((pt priority-table) obj)
  (gethash obj (pt-dictionary pt)))

;;; Coarse-grid construction

(defmethod choose-coarse-grid ((amg <stueben>) blackboard)
  "Ruge-Stueben coarse grid selection."
  (let* ((filtered-mat (getbb blackboard :filtered-matrix))
	 (undecided (make-instance 'priority-table))
	 (coarse ())
	 (fine ()))
    
    ;; distribute points to undecided
    (for-each-row-key
     #'(lambda (row-key)
	 (pt-insert undecided row-key
		    (aif (matrix-column filtered-mat row-key)
			 (hash-table-count it)
			 0)))
     filtered-mat)
    
    ;; loop through undecided building a tentative coarse grid by trying to
    ;; choose coarse-grid points such that as many fine-grid points as
    ;; possible depend on them.
    (loop for i = (pt-pop undecided) until (null i) do
	  (push i coarse)
	  (for-each-key-in-col
	   #'(lambda (j)
	       (when (pt-in-table-p undecided j)
		 (pt-remove undecided j)
		 (push j fine)
		 (for-each-key-in-row
		  #'(lambda (l)
		      (when (pt-in-table-p undecided l)
			(pt-shift-priority undecided l 1)))
		  filtered-mat j)))
	   filtered-mat i)
	  (for-each-key-in-row
	   #'(lambda (j)
	       (when (pt-in-table-p undecided j)
		 (pt-shift-priority undecided j -1)))
	   filtered-mat i))
    
    ;; We drop the second pass which is recommended in the 1987
    ;; Ruge-Stueben paper, because Stueben's new prolongation does not need
    ;; it according to Tanja.

    ;; augment result with coarsening data
    (setf (getbb blackboard :coarse-nodes) coarse
	  (getbb blackboard :fine-nodes) fine)
    ))


;;;; Testing

(defun test-stueben ()
  (let ((pt (make-instance 'priority-table)))
    (pt-insert pt 'A 1)
    (pt-insert pt 'B 1)
    (pt-remove pt 'B)
    (pt-pop pt)
    (describe pt))
    
  (let* ((A (laplace-sparse-matrix 3))
	 (b (make-image-vector-for A)))
    (fill-random! b 1.0)
    (let ((amg (make-instance '<stueben>
			      :max-depth 1 :cg-max-size 1 :output t)))
      (coarsen amg A)
      (solve (make-instance '<linear-solver> :iteration amg :success-if '(> :step 0) :output t)
	     (blackboard :problem (lse :matrix A :rhs b)))))
  )

;;; (test-stueben)
(fl.tests:adjoin-test 'test-stueben)