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

(in-package :mesh)

(defclass <skeleton> ()
  ((dim :accessor dimension :initarg :dimension :type (integer -1))
   (etables :accessor etables :type (array t (*))))
  (:documentation "A skeleton is a vector of hash-tables containing the
cells of a certain dimension as keys.  The information stored in the
values is different depending on the subclass derived from skeleton."))

(defmethod initialize-instance :after ((skel <skeleton>) &key cells &allow-other-keys)
  "The etables-list has to be initialized.  Furthermore, we allow a
constructor by giving a cell-list."
  ;; ensure dimension slot
  (unless (slot-boundp skel 'dim)
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
(defmethod etable ((skel <skeleton>) (dim fixnum))
  (aref (etables skel) dim))

(definline etable-of-highest-dim (skel)
  (etable skel (dimension skel)))

(defun nr-of-cells (skel &key dimension)
  "Returns number of cells in a skeleton."
  (cond
    ((numberp dimension) (hash-table-count (etable skel dimension)))
    ((eq dimension :top)
     (hash-table-count (etable skel (dimension skel))))
    (t (loop for etable across (etables skel)
	     summing (hash-table-count etable)))))

(defun cells-of-highest-dim (skel)
  (hash-table-keys (etable-of-highest-dim skel)))

(defun getskel (cell skel)
  (and (<= (dimension cell) (dimension skel))
       (gethash cell (etable skel (dimension cell)))))

(defun (setf getskel) (value cell skel)
  (and (<= (dimension cell) (dimension skel))
       (setf (gethash cell (etable skel (dimension cell)))
	     value)))

(defun member-of-skeleton? (cell skel)
  (nth-value 1 (getskel cell skel)))

(defun skel-ref (skel cell)
  (nth-value 0 (gethash cell (etable skel (dimension cell)))))

(defun (setf skel-ref) (value skel cell)
  (setf (gethash cell (etable skel (dimension cell))) value))

(defun get-cell-property (cell skel property)
  "Returns the value of the property."
  (getf (skel-ref skel cell) property))
(defun (setf get-cell-property) (value cell skel property)
  "Sets the value of the property."
  (setf (getf (skel-ref skel cell) property) value))

(defun insert-cell! (skel cell)
  "Inserts a cell and its boundary into a skeleton."
  (setf (skel-ref skel cell) nil)
  (loop for side across (boundary cell) do
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

(defun skel-empty-p (skel)
  (every (compose #'zerop #'hash-table-count) (etables skel)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; basic methods
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;(defmethod dimension ((skel <skeleton>)) (1- (length (etables skel))))

(defmethod manifold-dimension ((skel <skeleton>))
  (whereas ((cell (get-arbitrary-key-from-hash-table (etable skel 0))))
    (manifold-dimension cell)))

(defun skel-for-each (func skel &key direction dimension)
  "Loops through a skeleton applying func.  When direction is :down then loops
with dimension of the cells decreasing, otherwise increasing."
  (if (and dimension (>= dimension 0))
      (maphash func (etable skel dimension))
      (loop for etable across (if (eq direction :down)
				  (reverse (etables skel))
				  (etables skel))
	    do (maphash func etable))))

(defun skel-for-each-cell (func skel &key (direction :up) dimension)
  "Loops through a skeleton applying func on each cell.  When direction
is :down then loops with dimension of the cells decreasing, otherwise
increasing."
  (flet ((etable-for-each-cell (etable)
	   (loop for cell being each hash-key of etable do
		 (funcall func cell))))
    (if (and dimension (>= dimension 0))
	(etable-for-each-cell (etable skel dimension))
	(loop for etable across (if (eq direction :down)
				    (reverse (etables skel))
				    (etables skel))
	do (etable-for-each-cell etable)))))

(defmacro doskel ((looping-var skel &key (direction :up) where dimension) &body body)
  "Loops through a skeleton.  If looping-var is an atom, it loops through
all cells, otherwise it loops through cells and properties."
  (let ((looping-var (if (consp looping-var) looping-var (list looping-var)))
	(skel-for-each (if (consp looping-var) `skel-for-each `skel-for-each-cell))
	(skel-var (gensym "SKEL")))
    `(block nil
      (let ((,skel-var ,skel))
	(,skel-for-each
	 #'(lambda ,looping-var
	     ,@(if (eq where :surface)
		   `((unless (refined-p ,(car looping-var) ,skel-var)
		       ,@body))
		   body))
	 ,skel-var
	 :direction ,direction :dimension
	 ,(if (eq dimension :highest) `(dimension ,skel) dimension))))))


(definline skel-for-each-cell-of-highest-dimension (func skel)
  (skel-for-each-cell func skel :dimension (dimension skel)))

(defun skel-map (func skel)
  (let ((new-skel (make-analog skel)))
    (doskel ((cell value) skel)
      (setf (skel-ref new-skel cell) (funcall func cell value)))
    new-skel))

(defun find-cells (test skel &key dimension with-properties)
  "Returns a list of cells contained in skel and satisfying test."
  (let ((result ()))
    (doskel ((cell props) skel :dimension dimension)
      (when (if with-properties
		(funcall test cell props)
		(funcall test cell))
	(push cell result)))
    result))

(defun find-cell (test skel &rest rest)
  (let ((cells (apply #'find-cells test skel rest)))
    (assert (single? cells))
    (car cells)))

(defmethod find-cell-from-position ((skel <skeleton>) (pos array))
  (doskel (cell skel :dimension :highest)
    (when (inside-cell? cell pos)
      (return-from find-cell-from-position cell))))

(defun find-cell-from-corners (skel corners)
  (let ((cells (find-cells
		#'(lambda (cell)
		    (not (set-exclusive-or corners (corners cell) :test #'equalp)))
		skel)))
    (assert (single? cells))
    (car cells)))

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
	       (format stream " -> ~A~%" value)
	       (terpri stream)))
       skel :direction :up
       :dimension (when (eq *print-skeleton* :cells-of-highest-dimension)
		    (dimension skel))))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Checking skeletons
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmethod check ((skel <skeleton>))
  "Checks the skeleton.  An error is signaled if the skeleton appears bad."
  ;; check dimensions
  (loop with mfd-dim = (manifold-dimension skel)
	for dim from 0 to (dimension skel) do
	(loop for cell being the hash-keys of (etable skel dim) do
	      (check cell)
	      (cond
	       ((not (= (dimension cell) dim))
		(error "cell dimension does not fit"))
	       ((not (= (manifold-dimension cell) mfd-dim))
		(error "manifold-dimension does not fit")))))
  ;; check completeness
  (doskel (cell skel)
    (loop for subcell across (subcells cell) do
	  (assert (member-of-skeleton? subcell skel)))))

(defmethod skeleton-boundary ((skel <skeleton>))
  "Returns a skeleton consisting of cells of skel of dimension n-1 which
have only one neighbor."
  (if (zerop (dimension skel))
      (make-instance '<skeleton> :dimension -1)		; an empty skeleton
      (let ((bdry-ht
	     (loop with bdry-ht = (make-hash-table)
		   for cell being each hash-key of (etable skel (1- (dimension skel))) do
		   (setf (gethash cell bdry-ht) 0)
		   finally (return bdry-ht))))
	;; count neighbors for all hyperfaces
	(dohash (cell (etable skel (dimension skel)))
	  (loop for side across (boundary cell) do
		(dolist (side1 (or (cell-identification side skel)
				   (list side)))
		  (incf (gethash side1 bdry-ht)))))

	;; and build a skeleton from the boundary hyperfaces
	(make-instance
	 '<skeleton> :cells
	 (loop for cell being each hash-key of bdry-ht using (hash-value count)
	       when (= count 1) collecting cell
	       when (> count 2) do (error "Hyperface shared by more than two cells."))))))

(defmethod closed? ((skel <skeleton>))
  "Checks if the boundary of skel is empty.  This should be the case for
boundaries of skeletons."
  (minusp (dimension (skeleton-boundary skel))))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Cell-class activation part
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmemo refcell-skeleton (refcell)
  (assert (reference-cell? refcell))
  (skeleton refcell))

;; add this memoized functions to the cell-class activation procedure
(adjoin-cell-class-initialization-function 'refcell-skeleton)

;;;; Testing: (test-skeleton)

(defun test-skeleton ()
  "Tests are done when initializing the classes."
  (make-instance '<skeleton> :dimension -1)
  (macroexpand-1 '(doskel (cell (skeleton *unit-interval*) :direction :down)
		   (format t "~A~%" (corners cell))))
  (dimension (skeleton *reference-vertex*))
  (skeleton-boundary (skeleton *reference-vertex*))
  (make-instance '<skeleton> :dimension -1))
    
(tests:adjoin-femlisp-test 'test-skeleton)

