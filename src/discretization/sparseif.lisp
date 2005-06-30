;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; sparseif.lisp - Interface between discretization and algebra
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

(in-package :fl.discretization)

;;; This file provides the interface between finite-cell
;;; discretization and sparse-matrix representation.

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; generic function interface
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; First series of transfer routines: local vector arrays can be
;;; transfered to and from a global vector.  Changing the global
;;; vector can be done either by incrementing or setting.

(defgeneric cell->value-blocks (cell fe svec)
  (:documentation "Gets all value-blocks associated with the
subcells dofs.  If necessary, those value-blocks are generated."))

(defgeneric fill-local-from-global-vec (cell fe global-vec local-vec)
  (:documentation "Copies the region in global-vec determined by cell and
fe to local-vec."))

(defgeneric get-local-from-global-vec (cell fe global-vec)
  (:documentation "Maps the region in global-vec determined by cell and fe
to a local vector."))

(defgeneric set-global-to-local-vec (cell fe global-vec local-vec)
  (:documentation "Sets the region in global-vec determined by cell and fe
to the values of the local vector array."))

(defgeneric increment-global-by-local-vec (cell fe global-vec local-vec)
  (:documentation "Increments the region in global-vec determined by cell
and fe to the values of the local vector array."))

;;; Second group of transfer routines: local matrix arrays are
;;; transfered to and from a global sparse matrix.  Changing the
;;; global matrix can be done either by incrementing or setting.

(defgeneric cell->matrix-value-blocks (cell fe svec)
  (:documentation "Gets all value-blocks associated with the local
stiffness matrix.  If necessary, those value-blocks are generated."))

(defgeneric get-local-from-global-mat (cell fe global-mat)
  (:documentation "Maps the region in the global stiffness matrix
determined by cell and fe to a local matrix array."))

(defgeneric set-global-to-local-mat (cell fe global-mat local-mat)
  (:documentation "Sets the region in global-mat determined by cell and fe
to the values of the local matrix array."))

(defgeneric increment-global-by-local-mat (cell fe global-mat local-mat)
  (:documentation "Increments the region in global-mat determined by cell
and fe to the values of local-mat."))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; The actual implementation for <sparse-vector> and
;;;; <sparse-matrix>.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; Remark: Some routines might also be usable for interfaces to other data
;;; representations.

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; <sparse-vector> interface
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmethod cell->value-blocks ((cell <cell>) fe (svec <sparse-vector>))
  (let* ((mesh (mesh svec))
	 (subcells (vector-map (rcurry #'cell-key mesh) (subcells cell)))
	 (result (make-array (nr-of-subcells cell) :initial-element nil)))
    (loop for subcell-index in (subcell-indices fe) do
	  (setf (aref result subcell-index)
		(vref svec (aref subcells subcell-index))))
    result))

(defun do-fe-dofs-vblocks (vblocks fe local-vec operation &optional subcell-offset)
  (declare (type (simple-array t (*)) vblocks)
	   (type <scalar-fe> fe)
	   (type symbol operation)
	   (type (or null fixnum-vec) subcell-offset))
  (dotimes (i (ncols local-vec))
    (do-dof (dof fe)
      ;;(declare (optimize (speed 3) (safety 1)))
      (let* ((vblock-index (dof-subcell-index dof))
	     (vblock (aref vblocks vblock-index))
	     (in-vblock-index
	      (if subcell-offset
		  (the fixnum (+ (aref subcell-offset vblock-index)
				 (dof-in-vblock-index dof)))
		  (dof-in-vblock-index dof))))
	(symbol-macrolet ((local (mref local-vec (dof-index dof) i))
			  (global (mref vblock in-vblock-index i)))
	  (ecase operation
	    (:local<-global (setq local global))
	    (:global<-local (setq global local))
	    (:global+=local (incf global local))
	    (:global-=local (decf global local))))))))

(defmethod do-fe-dofs ((cell <cell>) (fe <scalar-fe>) (svec <sparse-vector>)
		       local-vec (operation symbol))
  (let ((vblocks (cell->value-blocks cell fe svec)))
    (do-fe-dofs-vblocks vblocks fe local-vec operation)))

(defmethod do-fe-dofs ((cell <cell>) (vecfe <vector-fe>) (svec <sparse-vector>)
		       (local-vec array) (operation symbol))
  (let ((vblocks (cell->value-blocks cell vecfe svec))
	(components (components vecfe))
	(subcell-offsets (subcell-offsets vecfe)))
    (dotimes (i (length components))
      (do-fe-dofs-vblocks vblocks (aref components i)
			  (aref local-vec i) operation
			  (aref  subcell-offsets i)))))

(defmethod fill-local-from-global-vec ((cell <cell>) fe (svec <sparse-vector>) local-vec)
  (do-fe-dofs cell fe svec local-vec :local<-global))

(defmethod get-local-from-global-vec ((cell <cell>) fe (svec <sparse-vector>))
  (let ((local-vec (make-local-vec fe (multiplicity svec))))
    (fill-local-from-global-vec cell fe svec local-vec)
    local-vec))

(defmethod set-global-to-local-vec ((cell <cell>) fe (svec <sparse-vector>) local-vec)
  (do-fe-dofs cell fe svec local-vec :global<-local))

(defmethod increment-global-by-local-vec ((cell <cell>) fe (svec <sparse-vector>) local-vec)
  (do-fe-dofs cell fe svec local-vec :global+=local))

(defun set-lagrange-ansatz-space-vector (asv func)
  "Sets an ansatz-space-vector to interpolate a given function.  This is
still a suboptimal implementation for vector functions."
  (doskel (cell (mesh asv) :where :surface :dimension :highest)
    (let* ((fe (get-fe (fe-class asv) cell))
	   (local-vec (make-local-vec fe (multiplicity asv))))
      (do-dof (dof fe)
	(let ((value (funcall func (local->global cell (dof-gcoord dof)))))
	  (when (numberp value) (setq value (make-real-matrix `((,value)))))
	  (dotimes (i (multiplicity asv))
	    (setf (mref local-vec (dof-index dof) i)
		  (mref (if (typep dof 'vector-dof)
			    (aref value (dof-component dof))
			    value)
			0 i)))))
      (set-global-to-local-vec cell fe asv local-vec))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; <sparse-matrix> interface
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmethod cell->matrix-value-blocks ((cell <cell>) fe (smat <sparse-matrix>))
  (let* ((mesh (mesh smat))
	 (subcells (vector-map (rcurry #'cell-key mesh) (subcells cell)))
	 (result (make-array (twice (nr-of-subcells cell)) :initial-element nil)))
    (loop for subcell-index-1 in (subcell-indices fe) do
	  (loop for subcell-index-2 in (subcell-indices fe) do
		(setf (aref result subcell-index-1 subcell-index-2)
		      (mref smat (aref subcells subcell-index-1)
			       (aref subcells subcell-index-2)))))
    result))

(defun do-fe-dofs-mblocks (mblocks fe1 fe2 local-mat operation
			   &optional subcell-offset1 subcell-offset2)
  ;; (declare (optimize (speed 3) (safety 1)))
  (declare (type (simple-array standard-matrix (* *)) mblocks)
	   (type <scalar-fe> fe1 fe2)
	   (type symbol operation)
	   (type (or null fixnum-vec) subcell-offset1 subcell-offset2))
  (do-dof (dof1 fe1)
    (do-dof (dof2 fe2)
      (let* ((mblock-index1 (dof-subcell-index dof1))
	     (mblock-index2 (dof-subcell-index dof2))
	     (mblock (aref mblocks mblock-index1 mblock-index2))
	     (in-mblock-index1
	      (if subcell-offset1
		  (the fixnum (+ (aref subcell-offset1 mblock-index1)
				 (dof-in-vblock-index dof1)))
		  (dof-in-vblock-index dof1)))
	     (in-mblock-index2
	      (if subcell-offset2
		  (the fixnum (+ (aref subcell-offset2 mblock-index2)
				 (dof-in-vblock-index dof2)))
		  (dof-in-vblock-index dof2))))
	(symbol-macrolet
	    ((local (mref local-mat (dof-index dof1) (dof-index dof2)))
	     (global (mref mblock in-mblock-index1 in-mblock-index2)))
	  (ecase operation
	    (:local<-global (setq local global))
	    (:global<-local (setq global local))
	    (:global+=local (incf global local))
	    (:global-=local (decf global local))))))))

(defmethod do-fe-dofs-mat ((cell <cell>) (fe <scalar-fe>) (smat <sparse-matrix>)
			   local-mat (operation symbol))
  (let ((mblocks (cell->matrix-value-blocks cell fe smat)))
    (do-fe-dofs-mblocks mblocks fe fe local-mat operation)))

(defmethod do-fe-dofs-mat ((cell <cell>) (vecfe <vector-fe>) (smat <sparse-matrix>)
			   (local-mat array) (operation symbol))
  (let ((mblocks (cell->matrix-value-blocks cell vecfe smat))
	(components (components vecfe))
	(subcell-offsets (subcell-offsets vecfe)))
    (dotimes (i (length components))
      (dotimes (j (length components))
	(do-fe-dofs-mblocks mblocks (aref components i) (aref components j)
			    (aref local-mat i j) operation
			    (aref subcell-offsets i) (aref subcell-offsets j))))))

(defmethod fill-local-from-global-mat ((cell <cell>) fe (smat <sparse-matrix>) local-mat)
  (do-fe-dofs-mat cell fe smat local-mat :local<-global))

(defmethod get-local-from-global-mat ((cell <cell>) fe (smat <sparse-matrix>))
  (let ((local-mat (make-local-mat fe)))
    (fill-local-from-global-mat cell fe smat local-mat)
    local-mat))

(defmethod set-global-to-local-mat ((cell <cell>) fe (smat <sparse-matrix>) local-mat)
  (do-fe-dofs-mat cell fe smat local-mat :global<-local))

(defmethod increment-global-by-local-mat ((cell <cell>) fe (smat <sparse-matrix>) local-mat)
  (do-fe-dofs-mat cell fe smat local-mat :global+=local))

(defmethod extract-ip-data ((cell <cell>) qrule property-list)
  "Converts all ansatz-space objects in the parameters list into local
value arrays corresponding to the finite element."
  (loop for (key object) on property-list by #'cddr
	when (typep object '<ansatz-space-vector>)
	collect key and collect
	(let* ((fe (get-fe (fe-class object) cell))
	       (local-vec (get-local-from-global-vec cell fe object)))
	  (map 'vector #'(lambda (shape-vals) (m*-tn shape-vals local-vec))
	       (ip-values fe qrule)))))

