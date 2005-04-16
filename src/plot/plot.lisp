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

(in-package :fl.plot)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Public interface
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defgeneric plot (object &key &allow-other-keys)
  (:documentation "Plot is a generic function which dispatches depending on
the type of object it receives.  Its behaviour can additionally be modified
by keyword parameters."))

(defparameter *plot* fl.graphic::*dx-pathname*
  "If set to NIL, plotting is disabled.  If set to :message, a message is
printed to *trace-output* instead of plotting.")

(defmethod plot :around (object &key &allow-other-keys)
  "Handles the *plot* parameter."
  (case *plot*
    ((t) (call-next-method))
    (:message (format *trace-output* "~&<Plotting suppressed>~%")))
  ;; and generally we return the object itself
  object)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; graphic-write-data
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun vertex-index-table (skel)
  "Returns a hash-table vertices->indices for refinement skeletons of
reference cells.  This is needed for plotting."
  (let ((index-table (make-hash-table))
	(index -1))
    (skel-for-each #'(lambda (vtx) (setf (gethash vtx index-table) (incf index)))
		   skel :dimension 0)
    index-table))

(let ((table (make-hash-table :test 'equalp)))
  (defun refcell-refinement-vertices (refcell level)
    "Transforms refcell-refinement-index-table into a vector."
    (assert (reference-cell-p refcell))
    (let ((skeleton (refcell-refinement-skeleton refcell level 0)))
      (whereas ((entry (gethash (cons refcell level) table)))
	;; safety check that nothing has changed
	(when (eq (car entry) skeleton)
	  (return-from refcell-refinement-vertices (cdr entry))))
      ;; generate a new entry
      (let* ((index-table (vertex-index-table skeleton))
	     (vertex-array (make-array (hash-table-count index-table))))
	(maphash #'(lambda (vtx index) (setf (aref vertex-array index) vtx))
		 index-table)
	;; because DX does not use connections in 1D
	(when (= (dimension refcell) 1)
	  (setf vertex-array
		(sort vertex-array #'<
		      :key (lambda (vtx) (aref (vertex-position vtx) 0)))))
	(setf (gethash (cons refcell level) table)
	      (cons skeleton vertex-array))
	vertex-array))))

(defun nr-of-refinement-vertices (cell depth)
  (length (refcell-refinement-vertices
	   (reference-cell cell)
	   (if (numberp depth) depth (funcall depth cell)))))

(definline make-key (cell local-vtx)
  "Returns a suitable hash-table key for the vertex of cell."
  (cons cell local-vtx))

(defun compute-position-indices (cells depth)
  "Collects all vertices of surface cells into a hash-table and indexes
them.  When type is :continuous, then the key is the vertex, otherwise the
key is the pair cell/vertex, which allows for the plotting of functions
that are discontinuous across cell boundaries.  @arg{depth} can be a number
representing the number of cell subdivisions, or a function which yields
this number for the actual cell."
  (let ((position-indices (make-hash-table :test 'equal))
	(global-index -1)
	(flag-1d (every (lambda (cell) (= (dimension cell) 1)) cells)))
    (when flag-1d
      (setq cells (sort (copy-seq cells) #'<
			:key (lambda (cell) (aref (midpoint cell) 0))))
      (dbg :plot "Cells=~A" cells))
    (dolist (cell cells)
     (loop with depth = (if (numberp depth) depth (funcall depth cell))
	for local-vtx across
	  (let ((vertices (refcell-refinement-vertices (reference-cell cell) depth)))
	    (when (= (dimension cell) 1)
	      ;; the following is a cludge because DX strangely breaks for
	      ;; unordered vertices
	      (when (minusp (mref (l2Dg cell #d(0.5)) 0 0))
		(setq vertices (reverse vertices))))
	    vertices)
	do
	  (setf (gethash (make-key cell local-vtx) position-indices)
		(incf global-index))))
    position-indices))

(defun connections (cells position-indices depth)
  "Returns the connections in a list.  Each list item is again a list
consisting of cell and vertex indices.  @arg{depth} can be a number
representing the number of cell subdivisions, or a function which yields
this number for the actual cell."
  (let ((connection-list ()))
    (dolist (cell cells connection-list)
      (let ((depth (if (numberp depth) depth (funcall depth cell))))
	(skel-for-each
	 #'(lambda (mini-cell)
	     (unless (simplex-p mini-cell)
	       (setq mini-cell (cell->cube mini-cell)))
	     (push
	      (mapcar
	       #'(lambda (local-vtx) (gethash (make-key cell local-vtx) position-indices))
	       (let ((vertices (vertices mini-cell)))
		 ;; check orientation, if necessary swap vertices
		 (if (and (simplex-p mini-cell)
			  (minusp (det (local->Dglobal
					mini-cell (local-coordinates-of-midpoint mini-cell)))))
		     (cons (cadr vertices) (cons (car vertices) (cddr vertices)))
		     vertices)))
	      connection-list))
	 (fl.mesh::refcell-refinement-skeleton (reference-cell cell) depth 0)
	 :dimension :highest)))))

(defun compute-all-position-values (cells position-indices depth cell->values)
  "Computes the values at the positions.  @arg{depth} can be a number
representing the number of cell subdivisions, or a function which yields
this number for the actual cell."
  (let* ((nr-positions (hash-table-count position-indices))
	 (data (make-double-vec nr-positions)))
    ;; generate data array
    (dolist (cell cells)
      (let* ((values (funcall cell->values cell))
	     (depth (if (numberp depth) depth (funcall depth cell)))
	     (vertices (refcell-refinement-vertices
			(reference-cell cell) depth)))
	(dotimes (i (length vertices))
	  (setf (aref data (gethash (make-key cell (aref vertices i))
				    position-indices))
		(vref values i)))))
    ;; return data field
    data))


(defun position-array (cells position-indices depth &optional transformation)
  "Collects all positions in the hash-table POSITION-INDICES into an
array."
  (let ((position-array (make-array (hash-table-count position-indices))))
    (dolist (cell cells)
      (let ((depth (if (numberp depth) depth (funcall depth cell))))
	(loop
	   for local-vtx across (refcell-refinement-vertices
				 (reference-cell cell) depth)
	   do
	   (setf (aref position-array (gethash (make-key cell local-vtx) position-indices))
		 (let ((pos (local->global cell (vertex-position local-vtx))))
		   (if transformation
		       (evaluate transformation pos)
		       pos))))))
      (dbg :plot "Positions=~A" position-array)
    position-array))


