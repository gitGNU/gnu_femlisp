;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; vanka.lisp
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

(in-package :geomg)

(defclass <local-bgs> (<block-gauss-seidel>)
  ((type :initarg :type))
  (:documentation "Iteration where the blocks in a block Gauss-Seidel
iteration are determined by combining the unknowns of several geometric
entities.  We provide two versions: cell-centered where all subcells of a
highest-dimensional cell constitute the block and vertex-centered where all
non-vertex cells surrounding a vertex constitute the subcell.  This latter
version can be shown to be robust in p, which seems to be a new result.  I
reported it first at the Oberwolfach meeting on fast solvers in June,
2003."))

(defun in-subcells-support (subcells cells)
  "Checks if one of the subcells is a subcell of one of the cells.
Usually, this will be one element lists.  The general case occurs for
identifications and hanging nodes."
  (some #'(lambda (cell)
	    (let ((cell-subcells (subcells cell)))
	      (some #'(lambda (subcell) (find subcell cell-subcells))
		    subcells)))
	cells))

(defun find-vertex-centered-block (vertex-key asa)
  "Collects a block of keys for cells containing the vertex.  Handles
identification and hanging nodes."
  (let ((block-keys (list vertex-key))
	(candidates (remove-if (rcurry #'slave-or-dirichlet-dof-p asa)
			       (keys-of-row asa vertex-key)))
	(constraints-Q (getf (structure-information (ansatz-space asa)) :constraints-Q)))
    ;; first direct neighbors are chosen
    (dolist (candidate candidates)
      (when (in-subcells-support (mklist vertex-key) (mklist candidate))
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
		      (in-subcells-support (mklist dependent-key) (mklist candidate)))
		  dependent-keys)
	(push candidate block-keys))))
    ;; return result
    block-keys))

(defmethod setup-blocks ((bgs <local-bgs>) (asa <sparse-matrix>))
  "Collects blocks consisting either of all subcells of cells in the
cell-centered case or all cells to which a vertex belongs in the
vertex-centered case."
  (dbg :iter "<local-bgs>-setup-blocks: starting")
  (let ((mesh-dim (dimension (mesh asa)))
	(type (slot-value bgs 'type))
	(blocks ()))
    (for-each-row-key
     #'(lambda (key)
	 (let ((block-keys ()))
	   ;; enlarge support under some circumstances
	   (if (slave-or-dirichlet-dof-p key asa)
	       (setq block-keys (list key))
	       (ecase type
		 (:cell-centered
		  (when (= (dimension (representative key)) mesh-dim)
		    (for-each-key-in-row
		     #'(lambda (col-key) (push col-key block-keys))
		     asa key)))
		 (:vertex-centered
		  (when (= (dimension (representative key)) 0)
		    (setq block-keys (find-vertex-centered-block key asa))))))
	   (when block-keys
	     (push (sort (coerce block-keys 'vector) #'>
			 :key (compose #'dimension #'representative))
		   blocks))))
     asa)
    (dbg :iter "<local-bgs>-blocks: ~A" blocks)
    blocks))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; <vanka>
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass <vanka> (<local-bgs>)
  ((type :initform :vertex-centered))
  (:documentation "New Vanka smoother for Q^{k+1}/Q^k discretizations of
Navier-Stokes."))

(defun extended-block (asa keys)
  "Collect the next layer around keys in the matrix graph of asa."
  (let ((inner (make-hash-table))
	(outer (make-hash-table)))
    ;; set inner hash table
    (loop for key across keys do (setf (gethash key inner) t))
    ;; if it is a standard non-slave block extend it by the outer layer
    (unless (slave-or-dirichlet-dof-p (aref keys 0) asa)
      (loop for key across keys do
	    (for-each-key-in-row
	     #'(lambda (col-key)
		 (unless (gethash col-key inner)
		   (setf (gethash col-key outer) t)))
	     asa key)))
    ;; finally return the block together with the components
    (let ((extended-keys (concatenate 'vector keys (hash-table-keys outer)))
	  (fe-class (fe-class asa)))
      (values
       extended-keys
       (vector-map
	#'(lambda (key)
	    (let ((fe (get-fe fe-class (representative key))))
	      (if (gethash key inner)
		  (cons 0 (nr-of-inner-dofs fe))
		  (let ((subcell-offsets (getf (properties fe) 'SUBCELL-OFFSETS)))
		    (cons 0 (aref (vector-last subcell-offsets) 0))))))
	extended-keys)))))

(defmethod setup-blocks ((vanka <vanka>) (asa <ansatz-space-automorphism>))
  "Adds to the usual vertex centered blocks all surrounding velocity
degrees of freedom."
  (let ((inner-blocks (call-next-method)) ; blocks from vertex-centered BGS
	(result-blocks ())
	(result-components ()))
    (dolist (inner inner-blocks)
      (multiple-value-bind (keys components)
	  (extended-block asa inner)
	(push keys result-blocks)
	(push components result-components)))
    (dbg :iter "Vanka blocks: ~A" result-blocks)
    (values result-blocks result-components)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; <heuveline-vanka>
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass <ns-vanka> (<local-bgs>)
  ()
  (:documentation "Special Vanka smoother used by Volker John and Vincent
Heuveline to solve Q^{k+1}/Q^k discretizations of Navier-Stokes."))


(defmethod setup-blocks ((vanka <ns-vanka>) (asa <ansatz-space-automorphism>))
  "Sets the block slot by collecting blocks consisting either of all
subcells in the cell-centered case or all cells to which a vertex belongs
in the vertex-centered case."
  (let ((fedisc (fe-class asa))
	(blocks ())
	(components ()))
    (for-each-row-key
     #'(lambda (key)
	 (unless (slave-or-dirichlet-dof-p key asa)
	   (let* ((fe (get-fe fedisc (representative key)))
		  (pressure-fe (vector-last (components fe))))
	     (when (plusp (nr-of-inner-dofs pressure-fe))
	       ;; in inner-block we collect small blocks of the form (key . (from . to))
	       (let ((inner-block (list (cons key (cons 0 (nr-of-inner-dofs fe))))))
		 (for-each-key-in-row
		  #'(lambda (col-key)
		      (unless (eq key col-key)
			(let* ((fe2 (get-fe fedisc (representative col-key)))
			       (subcell-offsets (getf (properties fe2) 'SUBCELL-OFFSETS))
			       (pressure-off (aref (vector-last subcell-offsets) 0)))
			  (push (cons col-key (cons 0 pressure-off)) inner-block))))
		  asa key)
		 (push (map 'vector #'car inner-block) blocks)
		 (push (map 'vector #'cdr inner-block) components))))))
     asa)
    (dbg :iter "Heuveline-Vanka blocks: ~D" (length blocks))
    (values blocks components)))

