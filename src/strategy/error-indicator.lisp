;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; error-indicator.lisp - Refinement indication
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

(in-package :strategy)

(defclass <refinement-indicator> ()
  ((ensure-mesh-quality :accessor ensure-mesh-quality
			:initform t :initarg :ensure-mesh-quality))
  (:documentation "An indicator is used as first argument in the generic
functions indicate which works on a blackboard.  Based on the
quantities computed by an error estimator, i.e. eta, indicate puts a list
of elements to be refined on the blackboard.  When ensure-mesh-quality
is t, the indicator ensures that the difference of mesh widths of
neighboring cells does not become larger than a factor of 4."))

(defgeneric indicate (indicator blackboard)
  (:documentation "Puts a list of elements to be refined on the blackboard."))

(defmethod indicate :around (indicator blackboard)
  "This around-method refines a cell if it detects a side that is refined
more than two times.  Thus, a maximum of two refinement levels between
side-adjacent cells is ensured."
  (call-next-method)
  (let ((mesh (getbb blackboard :mesh))
	(refinement-table (getbb blackboard :refinement-table)))
    (when (ensure-mesh-quality indicator)
      (doskel (cell mesh :dimension :highest :where :surface)
	(unless (gethash cell refinement-table)
	  (when (some #'(lambda (side)
			  (whereas ((children (children side mesh)))
			    (some (rcurry #'refined-p mesh) children)))
		      (boundary cell))
	    (setf (gethash cell refinement-table) t)))))
    ;; pass on the blackboard
    blackboard))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; uniform-refinement-indicator
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass <uniform-refinement-indicator> (<refinement-indicator>)
  ()
  (:documentation "Marks all cells for refinement."))

(defmethod indicate ((indicator <uniform-refinement-indicator>) blackboard)
  "Marks all cells for refinement which have no parent or for which the
error estimator yields a large eta."
  (with-items (&key mesh refinement-table)
      blackboard
    (let ((hash-table (make-hash-table)))
      (doskel (cell mesh :dimension :highest :where :surface)
	(setf (gethash cell hash-table) cell))
    (setf refinement-table hash-table))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; region-indicator
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass <region-indicator> (<refinement-indicator>)
  ((in-region :reader in-region :initarg :in-region))
  (:documentation "Marks all cells in a region for refinement."))

(defmethod indicate ((indicator <region-indicator>) blackboard)
  "Marks all cells for refinement which have no parent or for which the
error estimator yields a large eta."
  (with-items (&key mesh refinement-table)
      blackboard
    (let ((hash-table (make-hash-table))
	  (in-region (in-region indicator)))
      (doskel (cell mesh :dimension :highest :where :surface)
	(when (funcall in-region cell)
	  (setf (gethash cell hash-table) cell)))
    (setf refinement-table hash-table))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; largest-eta-indicator
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass <largest-eta-indicator> (<refinement-indicator>)
  ((fraction :accessor fraction :initform nil :initarg :fraction)
   (pivot-factor :accessor pivot-factor :initform 0.2 :initarg :pivot-factor)
   (from-level :accessor from-level :initform 0 :initarg :from-level)
   (block-p :reader block-p :initform nil :initarg :block-p))
  (:documentation "Puts the fraction of the cells with the largest error
contributions in the refinement table.  Note that a fraction of 1.0 yields
uniform refinement.  Below from-level, global refinement is used.  block-p
indicates that all children of a parent cell have to be refined at once."))

(defmethod indicate ((indicator <largest-eta-indicator>) blackboard)
  "Marks all cells for refinement which have no parent or for which the
error estimator yields a large eta."
  (with-items (&key refinement-table) blackboard
    (let ((mesh (getbb blackboard :mesh))
	  (eta (getbb blackboard :eta)))
      (with-slots (fraction pivot-factor from-level block-p)
	indicator
	(let ((hash-table (make-hash-table)))
	  (if (and (>= (top-level mesh) from-level) eta)
	      ;; error estimator could work out a value
	      (let ((sorted (sort (coerce (hash-table-keys eta) 'simple-vector) #'>
				  :key #'(lambda (key) (abs (gethash key eta))))))
		(loop with upper-part = (and fraction (ceiling (* (length sorted) fraction)))
		      with pivot = (and pivot-factor
					(* pivot-factor
					   (abs (gethash (aref sorted 0) eta))))
		      for i from 1
		      and cell across sorted
		      for parent = (parent cell mesh)
		      unless (gethash cell hash-table) do
		      (when (or (and upper-part (<= i upper-part))
				(and pivot (>= (abs (gethash cell eta)) pivot)))
			(setf (gethash cell hash-table) t)
			(when (and (block-p indicator) parent)
			  (loop for child across (children parent mesh)
				unless (refined-p child mesh) do
				(setf (gethash child hash-table) t))))))
	      ;; otherwise we mark everything globally
	      (for-each-cell-of-highest-dimension-on-surface
	       #'(lambda (cell) (setf (gethash cell hash-table) cell))
	       mesh))
	  (setf refinement-table hash-table))))))


;;;; Testing

(defun test-error-indicator ()
  (make-instance '<largest-eta-indicator>)
  )

(fl.tests:adjoin-test 'test-error-indicator)
