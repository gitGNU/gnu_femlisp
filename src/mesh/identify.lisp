;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; identify.lisp - Identification of cells in skeletons
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

(in-package :fl.mesh)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Identification
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass identification ()
  ((skeleton :initarg :skeleton)
   (cells :initarg :cells :documentation "A list of identified cells."))
  (:documentation "This object describes a cluster of identified cells.  At
the moment, it is used to describe identified boundaries and might be used
later on for parallelization purposes.  Those identifications are refered
to in the properties of each cell."))

(defmethod print-object :after ((id identification) stream)
  "Print the cells in the identification."
  (when *print-cell*
    (format stream "{~{~A~^,~}}" (slot-value id 'cells))))

(definline cell-identification (cell skel)
  "Returns @arg{cell}'s identification in @arg{skel} or NIL."
  (get-cell-property cell skel 'IDENTIFIED))

(definline (setf cell-identification) (identification cell skel)
  (setf (get-cell-property cell skel 'IDENTIFIED)
	identification))

(defun identified-p (cell skel)
  "Returns @arg{cell}'s identification in @arg{skel} or NIL."
  (cell-identification cell skel))

(defgeneric representative (obj)
  (:documentation "Returns a representative for this object.")
  (:method (obj)
    "This default method returns the object itself."
    obj)
  (:method ((id identification))
    "Gets a representative cell for this identification."
    (first (slot-value id 'cells))))

(defmethod identified-cells (cell skel)
  (aif (cell-identification cell skel)
       (slot-value it 'cells)
       (list cell)))

(defun identify (identified-cells skel)
  "Identifies all cells in @arg{identified-cells} within @arg{skel}."
  (dbg :mesh "Identifying: ~A" identified-cells)
  (when (> (length identified-cells) 1)
    (let ((identification
	   (make-instance 'identification :skeleton skel :cells identified-cells)))
      (dolist (cell identified-cells)
	(setf (cell-identification cell skel) identification)))))

(defun key-is-subcell-p (key1 key2)
  "Checks if @arg{key1} occurs as subcell of @arg{key2}.  The keys can be
either cells or identifications."
  (let ((cell-subcells (subcells (representative key2))))
    (some #'(lambda (subcell) (find subcell cell-subcells))
	  (typecase key1
	    (identification (slot-value key1 'cells))
	    (t (list key1))))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Helper routines for constructing a cell domain with identified
;;; boundaries
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun boundary-identifications (identifications)
  "Returns an identification list for the boundaries of the cells in
identifications."
  (mappend
   #'(lambda (identified)
       (apply #'map 'list #'list (mapcar #'boundary identified)))
   identifications))

(defun combine-identifications (sets)
  "Reduces identifications to maximally connected sets."
  (and sets
       (multiple-value-bind (set disconnected)
	   (maximally-connected (car sets) (cdr sets) :test #'intersection :combine #'union)
	 (cons set (combine-identifications disconnected)))))

(defun iterate-identifications (initial-identifications)
  "Generates all identifications of the skeleton from the identifications
of some higher-dimensional cells."
  (loop for identifications = initial-identifications
	then (combine-identifications (boundary-identifications identifications))
	until (null identifications)
	appending identifications))

(defun check-identification (skel)
  "Should checks correct identification.  At the moment, it tests only if
all identification-entries are eq."
  (doskel (cell skel)
    (whereas ((id (cell-identification cell skel)))
      (let ((cells (slot-value id 'cells)))
	(assert (not (single? cells)))
	(dolist (cell2 cells)
	  (assert (eq id (cell-identification cell2 skel))))))))

(defun identify-unit-cell-faces (skel &key (indices :all))
  "This routines identifies boundary cells in skel which correspond to
boundary cells in the unit cube.  Warning: exact arithmetic is used to
recognize identified cells.  This should work for skeletons derived from
the unit cell, but may create problems in other situations."
  (let ((table (make-hash-table :test 'equalp)))
    ;; sort cells wrt a pivot
    (doskel (side (skeleton-boundary skel))
      (let ((key
	     (if (eq indices :all)
		 (vector-map #'(lambda (coord) (if (= coord 1.0) 0.0 coord))
			     (midpoint side))
		 (loop with result = (copy-seq (midpoint side))
		       for index in indices
		       when (= (aref result index) 1.0)
		       do (setf (aref result index) 0.0)
		       finally (return result)))))
	;; identify cells on the unit cube boundary
	(when (some #'zerop key)
	  (push side (gethash key table)))))
    ;; identify all lists of identified cells
    (mapc (rcurry #'identify skel) (hash-table-values table))
    (check-identification skel)
    skel))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Interfacing
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass <interface> (<skeleton>)
  ((neighbors :accessor neighbors :initarg :neighbors
	      :type (simple-array <skeleton> (*))))
  (:documentation "Not in use up to now.  Interface structure for domain
decomposition approaches.  Might be used also for parallelization in the
future.  Cells are mapped to identified cells in the neighboring
skeletons."))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Tests
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(in-suite mesh-suite)

(test identify
  ;;; nothing yet
  )
