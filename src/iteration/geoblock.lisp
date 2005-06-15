;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; geoblock.lisp
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

(in-package :fl.geomg)

(defclass <geometric-blocking-mixin> ()
  ((type :initform :vertex-centered :initarg :type))
  (:documentation "Determines if the block choice is centered on cells or
vertices.  The latter choice can be shown to be robust in the order p for
finite elements of Lagrange type [Pavarino 1994]."))
  
(defclass <geometric-psc> (<geometric-blocking-mixin> <psc>)
  ()
  (:documentation "PSC with geometry-based block choice."))

(defclass <geometric-ssc> (<geometric-blocking-mixin> <ssc>)
  ()
  (:documentation "SSC with geometry-based block choice."))

;;; Constructors

(defun geometric-psc (&rest rest)
  "Constructor of a geometric parallel subspace correction."
  (apply #'make-instance '<geometric-psc> rest))
		 
(defun geometric-ssc (&rest rest)
  "Constructor of a geometric successive subspace correction."
  (apply #'make-instance '<geometric-ssc> rest))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Block setup for geometric subspace correction schemes
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun find-vertex-centered-block (vertex-key asa)
  "Collects a block of keys for cells containing the vertex.  Handles
identification and hanging nodes."
  (let ((block-keys (list vertex-key))
	(candidates (remove-if (rcurry #'slave-or-dirichlet-dof-p asa)
			       (keys-of-row asa vertex-key)))
	(constraints-Q (getf (properties (ansatz-space asa)) :constraints-Q)))
    ;; first direct neighbors are chosen
    (dolist (candidate candidates)
      (when (key-is-subcell-p vertex-key candidate)
	(setq block-keys (adjoin candidate block-keys))))
    ;; remove the positive choices
    (setq candidates (set-difference candidates block-keys))
    ;; next we check if dependent cells are subcells
    (let ((dependent-keys
	   (reduce #'union block-keys :initial-value ()
		   :key (curry #'keys-of-column constraints-Q))))
      ;; second loop checks if some candidate has a subcell in dependent-keys
    (dolist (candidate candidates)
      (when (some #'(lambda (dependent-key)
		      (key-is-subcell-p dependent-key candidate))
		  dependent-keys)
	(push candidate block-keys))))
    ;; return result
    block-keys))

(defun find-connected-blocks-in-table (table asa)
  (flet ((extend (keys)
	   (let ((new ()))
	     (dolist (rk keys (append new keys))
	       (for-each-key-in-row
		#'(lambda (ck)
		    (when (gethash ck table)
		      (unless (member ck keys)
			(push ck new)
			(remhash ck table))))
		asa rk)))))
    (loop until (zerop (hash-table-count table)) collecting
	  (loop with key = (get-arbitrary-key-from-hash-table table)
		with keys = (list key)
		initially (remhash key table)
		for extend = (extend keys)
		until (eq extend keys) do
		(setq keys extend)
		finally (return (coerce keys 'vector))))))

(defmethod setup-blocks ((blockit <geometric-blocking-mixin>) (asa <sparse-matrix>))
  "Collects blocks consisting either of all subcells of cells in the
cell-centered case or all cells to which a vertex belongs in the
vertex-centered case."
  (let ((mesh-dim (dimension (mesh asa)))
	(type (slot-value blockit 'type))
	(blocks ())
	(remaining (make-hash-table)))
    (for-each-row-key
     #'(lambda (key)
	 (let ((block-keys ()))
	   ;; enlarge support under some circumstances
	   (cond
	     ((slave-dof-p key asa)
	      (setq block-keys (list key)))
	     #+(or) ((dirichlet-dof-p key asa)  ; is treated as everything else
		     (setq block-keys (list key)))
	     ((= (dimension (representative key))
		 (ecase type (:vertex-centered 0) (:cell-centered mesh-dim)))
	      (ecase type 
		 (:cell-centered
		  (for-each-key-in-row
		   #'(lambda (col-key) (push col-key block-keys))
		   asa key))
		 (:vertex-centered
		  (setq block-keys (find-vertex-centered-block key asa)))))
	     (t ;; may remain at the boundaries (e.g. on triangles)
	      (setf (gethash key remaining) t)))
	   (when block-keys
	     (push (sort (coerce block-keys 'vector) #'>
			 :key (compose #'dimension #'representative))
		   blocks))))
     asa)
    ;; remove all from remaining which are covered
    (loop for block in blocks do
	  (loop for key across block do
		(remhash key remaining)))
    (dohash (key remaining)
      (assert (dirichlet-dof-p key asa)))
    #+(or)(setq blocks (nconc blocks (find-connected-blocks-in-table remaining asa)))
    blocks))

