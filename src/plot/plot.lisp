;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; plot.lisp - General plotting routines
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

(in-package :plot)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Public interface
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defgeneric plot (object &rest rest &key &allow-other-keys)
  (:documentation "Plot is a generic function which dispatches depending on
the type of object it receives.  Its behaviour can additionally be modified
by keyword parameters."))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; graphic-write-data
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(definline make-key (cell local-vtx)
  "Returns a suitable hash-table key for the vertex of cell."
  (cons cell local-vtx))

(defmemo sorted-1d-refinement-vertices (depth)
  "Because DX does not use connections in 1D, we sort our input."
  (sort (copy-seq (refcell-refinement-vertices *unit-interval* depth))
	#'< :key (lambda (vtx) (aref (vertex-position vtx) 0))))

(defun compute-position-indices (cells depth)
  "Collects all vertices of surface cells into a hash-table and indexes
them.  When type is :continuous, then the key is the vertex, otherwise the
key is the pair cell/vertex, which allows for the plotting of functions
that are discontinuous across cell boundaries."
  (let ((position-indices (make-hash-table :test 'equal))
	(global-index -1)
	(flag-1d (every (lambda (cell) (= (dimension cell) 1)) cells)))
    (when flag-1d
      (setq cells (sort (copy-seq cells) #'<
			:key (lambda (cell) (aref (midpoint cell) 0)))))
    (dolist (cell cells)
     (loop for local-vtx across
	   (if (= (dimension cell) 1)
	       (sorted-1d-refinement-vertices depth)
	       (refcell-refinement-vertices (reference-cell cell) depth))
	   do
	   (setf (gethash (make-key cell local-vtx) position-indices)
		 (incf global-index))))
    position-indices))

(defun connections (cells position-indices depth)
  "Returns the connections in a list.  Each list item is again a list
consisting of cell and vertex indices."
  (let ((connection-list ()))
    (dolist (cell cells connection-list)
      (skel-for-each-cell-of-highest-dimension
       #'(lambda (mini-cell)
	   (push
	    (mapcar
	     #'(lambda (local-vtx) (gethash (make-key cell local-vtx) position-indices))
	     (let ((vertices (vertices mini-cell)))
	       ;; check orientation, if necessary swap vertices
	       (if (and (simplex? mini-cell)
			(minusp (det (local->Dglobal
				      mini-cell (local-coordinates-of-midpoint mini-cell)))))
		   (cons (cadr vertices) (cons (car vertices) (cddr vertices)))
		   vertices)))
	    connection-list))
       (mesh::refcell-refinement-skeleton (reference-cell cell) depth)))))

(defun compute-all-position-values (cells position-indices depth cell->values)
  "Computes the values at the positions."
  (let* ((nr-positions (hash-table-count position-indices))
	 (data (make-double-vec nr-positions)))
    ;; generate data array
    (loop for cell in cells
	  for values = (funcall cell->values cell)
	  for vertices = (refcell-refinement-vertices
			  (reference-cell cell) depth)
	  do
	  (dotimes (i (length vertices))
	    (setf (aref data (gethash (make-key cell (aref vertices i))
				      position-indices))
		  (matrix-ref values i))))
    ;; return data field
    data))


(defun position-array (cells position-indices depth)
  "Collects all positions into an array.  For the moment, we ignore hanging
nodes."
  (let ((position-array (make-array (hash-table-count position-indices))))
    (dolist (cell cells position-array)
      (loop
       for local-vtx across (refcell-refinement-vertices
			     (reference-cell cell) depth)
       do
       (setf (aref position-array (gethash (make-key cell local-vtx) position-indices))
	     (local->global cell (vertex-position local-vtx)))))))


