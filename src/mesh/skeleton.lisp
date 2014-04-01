;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;  skeleton.lisp
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

(defclass <skeleton> ()
  ((dimension :accessor dimension :initarg :dimension :type (integer -1))
   (etables :accessor etables))
  (:documentation "A skeleton is a vector of hash-tables containing the
cells of a certain dimension as keys.  The information stored in the
values is different depending on the subclass derived from skeleton."))

(defmethod initialize-instance :after ((skel <skeleton>) &key cells &allow-other-keys)
  "The etables-list has to be initialized.  Furthermore, we allow a
constructor by giving a cell-list."
  ;; ensure dimension slot
  (unless (slot-boundp skel 'dimension)
    (setf (dimension skel) (apply #'max -1 (mapcar #'dimension cells))))
  ;; setup cell-tables
  (loop with tables = (make-array (1+ (dimension skel)))
	for i from 0 upto (dimension skel) do
	(setf (aref tables i) (make-hash-table :test #'eq))
	finally (setf (etables skel) tables))
  ;; and insert cells
  (insert-cells! skel cells))

;;; analog constructor
(defmethod make-analog ((skel <skeleton>))
  (make-instance '<skeleton> :dimension (dimension skel)))

;;; an extended access
(defmethod etable ((skel <skeleton>) dim)
  (aref (etables skel) dim))

(inlining
 (defun etable-of-highest-dim (skel)
   (etable skel (dimension skel))))

(defun member-of-skeleton? (cell skel)
  "Returns T if @arg{cell} is in @arg{skel}, NIL otherwise."
  (and (<= (dimension cell) (dimension skel))
       (nth-value 1 (gethash cell (etable skel (dimension cell))))))

(defun skel-ref (skel cell)
  "Returns the properties of @arg{cell} in @arg{skel}."
  (nth-value 0 (gethash cell (etable skel (dimension cell)))))

(defun (setf skel-ref) (value skel cell)
  "Setter for the properties of @arg{cell} in @arg{skel}."
  (setf (gethash cell (etable skel (dimension cell))) value))

(definline get-cell-property (cell skel property)
  "Returns the value of the property."
  (getf (skel-ref skel cell) property))
(definline (setf get-cell-property) (value cell skel property)
  "Sets the value of the property."
  (setf (getf (skel-ref skel cell) property) value))

;;; Identification

(defgeneric identified-cells (cell skel)
  (:documentation "Returns a list of cells in @arg{skel} which are
identified with @arg{cell}."))

;;; Building tools

(defun insert-cell! (skel cell &optional properties)
  "Inserts @arg{cell} and if necessary also its boundary into @arg{skel}.
If properties are given those are used for @arg{cell}."
  (setf (skel-ref skel cell) properties)
  (loop for side across (boundary cell)
     unless (member-of-skeleton? side skel) do
       (insert-cell! skel side)))

(defun insert-cells! (skel cells)
  "Inserts a list of cells into a skeleton."
  (dolist (cell cells skel)
    (insert-cell! skel cell)))

;;; further constructors
(defun cells->skeleton (cells)
  (let ((skel (make-instance '<skeleton> :dimension
			     (apply #'max -1 (mapcar #'dimension cells)))))
    (insert-cells! skel cells)
    skel))

(defmethod skeleton ((cells sequence))
  (make-instance '<skeleton> :cells (coerce cells 'list)))
(defmethod skeleton ((cell <cell>))
  (skeleton (list cell)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; basic methods
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;(defmethod dimension ((skel <skeleton>)) (1- (length (etables skel))))

(defmethod skel-for-each (func (skel <skeleton>) &key direction dimension where with-properties)
  "Loops through a skeleton applying func.  When direction is :down then loops
with dimension of the cells decreasing, otherwise increasing."
  (flet ((etable-for-each (etable)
	   (maphash #'(lambda (cell props)
			(unless (and (eq where :surface)
				     (refined-p cell skel))
			  (if with-properties
			      (funcall func cell props)
			      (funcall func cell))))
		    etable)))
    (unless (minusp (dimension skel))  ; nothing there
      (when (eq dimension :highest)
        (setf dimension (dimension skel)))
      (multiple-value-bind (first-dim last-dim)
          (cond ((consp dimension)
                 (values (car dimension) (cdr dimension)))
                ((numberp dimension) (values dimension dimension))
                (t (values 0 (dimension skel))))
        (let ((step (if (eq direction :down) -1 1)))
          (loop repeat (1+ (- last-dim first-dim))
             for dim = (if (eq direction :down) last-dim first-dim)
             then (+ dim step) do
               (etable-for-each (etable skel dim))))))))

(defmacro doskel ((looping-var skel &key (direction :up) where dimension) &body body)
  "Loop through a skeleton.  If looping-var is an atom, it loops through
all cells, otherwise it loops through cells and properties."
  (let* ((with-properties (listp looping-var))
	 (arguments (if with-properties looping-var (list looping-var))))
    `(block nil
      (skel-for-each #'(lambda ,arguments ,@body)
       ,skel :direction ,direction :dimension ,dimension :where ,where
       :with-properties ,with-properties))))

(defun nr-of-cells (skel &optional dimension)
  "Returns number of cells in a skeleton."
  (mapper-count #'skel-for-each skel :dimension dimension))

(defun skel-empty-p (skel)
  (zerop (nr-of-cells skel)))

(defun cells-of-dim (skel dim)
  "Returns the cells of @arg{skel} of dimension @arg{dim} in form of a
list."
  (mapper-collect #'skel-for-each skel :dimension dim))

(defun cells-of-highest-dim (skel)
  "Returns the cells of @arg{skel} of highest dimension in form of a list."
  (cells-of-dim skel (dimension skel)))

(defun mark-skeleton (skel prop value)
  "Marks all cells of @arg{skel} with the given @arg{prop}/@arg{value}
pair."
  (doskel (cell skel)
    (setf (get-cell-property cell skel prop) value))
  skel)

(defun skel-map (func skel)
  (lret ((new-skel (make-analog skel)))
    (doskel ((cell value) skel)
      (setf (skel-ref new-skel cell) (funcall func cell value)))))

(defun find-cells (test skel &key dimension with-properties where)
  "Returns a list of cells contained in skel and satisfying test."
  (lret ((result ()))
    (doskel ((cell props) skel :dimension dimension :where where)
      (when (if with-properties
		(funcall test cell props)
		(funcall test cell))
	(push cell result)))))

(defun find-cell (test skel &rest rest)
  (let ((cells (apply #'find-cells test skel rest)))
    (assert (<= (length cells) 1))
    (car cells)))

(defgeneric find-cell-from-position (skel pos)
  (:documentation
   "Finds a cell from @arg{skel} which contains the global position @arg{pos}."))

(defmethod find-cell-from-position ((skel <skeleton>) (pos array))
  (doskel (cell skel :dimension :highest)
    (when (inside-cell? cell pos)
      (return-from find-cell-from-position cell))))

(defun find-cell-from-corners (skel corners)
  (let ((cells (find-cells
		#'(lambda (cell)
		    (not (set-exclusive-or corners (corners cell) :test #'equalp)))
		skel)))
    (assert (<= (length cells) 1))
    (car cells)))

(defmethod embedded-dimension ((skel <skeleton>))
  (lret (dim)
    (doskel (cell skel)
      (let ((dim1 (embedded-dimension cell)))
        (if dim
            (unless (= dim dim1)
              (setq dim nil)
              (return))
            (setq dim dim1))
        (unless *check-well-defined-embedded-dimension*
          (return))))))

;;;; Parts

(defgeneric dimension-of-part (skel part)
  (:documentation
   "Parts of a skeleton can be named with the property @symbol{:part}.")
  (:method ((skel <skeleton>) part)
    (lret ((dim -1))
      (doskel (cell skel)
        (when (eql (get-cell-property cell skel :part) part)
          (setq dim (max (dimension cell) dim)))))))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Printing
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defparameter *print-skeleton* :cells
  "May be nil, :cells, :cells-of-highest-dimension.")
(defparameter *print-skeleton-values* t
  "When t, prints also the values of a skeleton.")

(defmethod describe-object :after ((skel <skeleton>) stream)
  (unless (minusp (dimension skel))
    (when *print-skeleton*
      (format stream "~&Cells:~%")
      (skel-for-each
       #'(lambda (cell value)
	   (princ cell stream)
	   (if *print-skeleton-values*
	       (format stream " ->~%~{  ~S ~S~%~}~%" value)
	       (terpri stream)))
       skel :with-properties t :direction :up :dimension
       (when (eq *print-skeleton* :cells-of-highest-dimension)
	 (dimension skel))))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Checking skeletons
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defparameter *optional-skeleton-checks*
  '(:embedded-dimension)
  "Some optional checks for skeletons.")

(defmethod check progn ((skel <skeleton>))
  "Checks the skeleton.  An error is signaled if the skeleton looks bad."
  ;; check dimensions
  (let ((emb-dim (let ((*check-well-defined-embedded-dimension* t))
                   (embedded-dimension skel))))
    (when (member :embedded-dimension *optional-skeleton-checks*)
      (assert emb-dim () "Embedded dimensions differ for cells."))
    (loop for dim from 0 upto (dimension skel) do
         (doskel (cell skel :dimension dim)
           (check cell)
           (cond
             ((not (= (dimension cell) dim))
              (error "cell dimension does not fit"))
             ((and emb-dim (not (= (embedded-dimension cell) emb-dim)))
              (error "embedded-dimension does not fit"))))))
  ;; check completeness
  (doskel (cell skel)
    (loop for subcell across (subcells cell) do
         (assert (member-of-skeleton? subcell skel)))))

(defun check-for-multiples (skel &key (output t))
  "Checks for different vertices at the same location and
for cells of higher dimensions connecting the same vertices
in a possibly different order.
Warning: Up to now a quadratic algorithm is used."
  ;; warning: quadratic algorithm
  (let (multiples)
    ;; find multiple vertices
    (let* ((vertices (coerce (cells-of-dim skel 0) 'vector))
           (n (length vertices)))
      (loop for i below n
         for x = (midpoint (aref vertices i))
         do
         (loop for j from (1+ i) below n
            for y = (midpoint (aref vertices j))
            do
            (when (mzerop (m- x y) 1.0e-8)
              (push (cons x y) multiples)))))
    ;; find multiple higher-dimensional cells
    (loop for k from 1 upto (dimension skel)
       for cells = (cells-of-dim skel k)
       for boundaries = (mapcar (lambda (cell)
                                  (cons cell (coerce (boundary cell) 'list)))
                                cells)
       do
       (loop for ((cell1 . bdry1) . others) on boundaries do
            (loop for (cell2 . bdry2) in others do
                 (unless (set-difference bdry1 bdry2)
                   (push (cons cell1 cell2) multiples)))))
    (setf multiples (nreverse multiples))
    ;; output
    (when output
      (loop for (cell1 . cell2) in multiples do
           (format t "~A~30T~A~%" cell1 cell2)))
    multiples))

(defmethod skeleton-substance ((skel <skeleton>))
  "The substance of a skeleton are those cells which are not boundary of a
cell in @arg{cell}.  This function returns those cells in the form of a
hash-table."
  (lret ((non-bdry-cells (make-hash-table)))
    (doskel (cell skel)
      (setf (gethash cell non-bdry-cells) t))
    (doskel (cell skel)
      (loop for side across (boundary cell) do
           (remhash side non-bdry-cells)))))

(defmethod skeleton-boundary ((skel <skeleton>))
  "Returns a skeleton consisting of cells of skel of dimension n-1 which
have only one neighbor."
  (when (zerop (dimension skel))
    (return-from skeleton-boundary
      (make-instance '<skeleton> :dimension -1))) ; an empty skeleton
  ;; find boundary of the "flesh" and build up a table of those cells to
  ;; the flesh cells
  (let ((bdry-ht (make-hash-table)))
    (loop for cell being each hash-key of (skeleton-substance skel) do
         (loop for side across (boundary cell) do
              (dolist (side1 (identified-cells side skel))
                (push cell (gethash side1 bdry-ht ())))))
    ;; and build a skeleton from the boundary hyperfaces
    (make-instance
     '<skeleton> :cells
     (loop for face being each hash-key of bdry-ht using (hash-value neighbors)
        for nr-neighbors = (length neighbors)
        when (= nr-neighbors 1) collecting face
        when (> nr-neighbors 2) do
        (error "Hyperface ~A shared by more than two cells:~%~A" face neighbors)))))

(defmethod closed? ((skel <skeleton>))
  "Checks if the boundary of skel is empty.  This should be the case for
boundaries of skeletons."
  (skel-empty-p (skeleton-boundary skel)))

;;;; Testing:

(defun test-skeleton ()
  "Tests are done when initializing the classes."
  (make-instance '<skeleton> :dimension -1)
  (macroexpand-1 '(doskel (cell (skeleton *unit-interval*) :direction :down)
		   (format t "~A~%" (corners cell))))
  (let ((v1 (make-vertex #d(0.0)))
        (v2 (make-vertex #d(0.0 0.0)))
        (*check-well-defined-embedded-dimension* t))
    (let ((line (make-line v1 v2)))
      (assert (= (embedded-dimension v1) 1))
      (assert (null (embedded-dimension line))))
    (let ((skel (make-instance '<skeleton> :dimension 0)))
      (insert-cells! skel (list v1 v2))
      (assert (null (embedded-dimension skel)))))
  )

;;; (test-skeleton)
(fl.tests:adjoin-test 'test-skeleton)

