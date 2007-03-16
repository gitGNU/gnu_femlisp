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
;;;; interface
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; generation of local vectors and matrices

(defun make-local-vec (svec cell)
  "Generates a local vector for local discretization."
  (let ((fe (get-fe (ansatz-space svec) cell))
	(multiplicity (multiplicity svec)))
    (map 'simple-vector
	 #'(lambda (fe) (make-real-matrix (nr-of-dofs fe) multiplicity))
	 (components fe))))

(defun make-local-mat (smat cell)
  "Generates a local matrix for local discretization."
  (let* ((fe-1 (get-fe (image-ansatz-space smat) cell))
	 (fe-2 (get-fe (domain-ansatz-space smat) cell))
	 (comps-1 (components fe-1))
	 (comps-2 (components fe-2)))
    (lret ((result (make-array (list (length comps-1) (length comps-2)))))
      (for-each-key
       (lambda (i j)
	 (setf (aref result i j)
	       (make-real-matrix (nr-of-dofs (aref comps-1 i))
				 (nr-of-dofs (aref comps-2 j)))))
       result))))

;;; transfer between local and global vector

(defgeneric local-value-blocks (cell svec)
  (:documentation "Gets all value-blocks associated with the subcells.
If necessary, those value-blocks are generated."))

(defgeneric fill-local-from-global-vec (cell global-vec local-vec)
  (:documentation "Copies the region in global-vec determined by cell to
local-vec."))

(defgeneric get-local-from-global-vec (cell global-vec)
  (:documentation "Maps the region in global-vec determined by cell to a
local vector."))

(defgeneric set-global-to-local-vec (cell global-vec local-vec)
  (:documentation "Sets the region in global-vec determined by cell to the
values of the local vector array."))

(defgeneric increment-global-by-local-vec (cell global-vec local-vec)
  (:documentation "Increments the region in global-vec determined by cell
to the values of the local vector array."))

;;; transfer between local and global matrix

(defgeneric get-local-from-global-mat (cell global-mat)
  (:documentation "Maps the region in the global stiffness matrix
determined by cell to a local matrix array."))

(defgeneric set-global-to-local-mat (cell global-mat local-mat)
  (:documentation "Sets the region in global-mat determined by cell to the
values of the local matrix array."))

(defgeneric fill-local-from-global-mat (cell global-vec local-vec)
  (:documentation "Copies the region in global-mat determined by cell to
local-mat."))

(defgeneric increment-global-by-local-mat (cell global-mat local-mat)
  (:documentation "Increments the region in global-mat determined by cell
to the values of local-mat."))

;;; The actual implementation for <sparse-vector> and <sparse-matrix>.

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; <sparse-vector> interface
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; the following routine should be extended for extracting subvectors, if
;;; partial equations are to be assembled on subdomains

(with-memoization (:type :global)
  (defun fe-secondary-information (fe)
    "Computes for a (vector) finite element the secondary information which
is important when extracting data from ansatz-space vectors and matrices.
For each dof in the fe, the access information to the global matrix in the
form vblock/in-vblock-index and the access information to the local matrix
in the form component-index/in-component-index is computed."
    (memoizing-let ((fe fe))
      (let* ((nr-dofs (reduce #'+ (components fe) :key #'nr-of-dofs))
	     (nr-subcells (length (subcells (reference-cell fe)))) 
	     (component-index (make-array nr-dofs))
	     (in-component-index (make-array nr-dofs))
	     (vblock-index (make-array nr-dofs))
	     (in-vblock-index (make-array nr-dofs))
	     (in-vblock-pos (make-array nr-subcells :initial-element 0))
	     (k 0))
	(loop+ ((scalar-fe (components fe)) comp) do
	       (loop+ (in-comp (dof (fe-dofs scalar-fe))) do
		      (let ((sci (dof-subcell-index dof)))
			(setf (aref component-index k) comp
			      (aref in-component-index k) in-comp
			      (aref vblock-index k) sci
			      (aref in-vblock-index k) (aref in-vblock-pos sci)
			      k (1+ k))
			(incf (aref in-vblock-pos sci))
			)))
	(assert (= k nr-dofs))  ; consistency check
	(list :nr-of-components (nr-of-components fe) :nr-dofs nr-dofs
	      :component-index component-index :in-component-index in-component-index
	      :vblock-index vblock-index :in-vblock-index in-vblock-index)))))

(with-memoization (:type :global :test 'equalp)
  (defun fe-extraction-information (fe indices)
    "Computes information for extracting components out of a vector finite
 element."
    (memoizing-let ((components (components fe)) (indices indices))
      (assert (and (every (lambda (i) (< i (length components))) indices)
		   (<= (length indices) (length components))
		   (equalp indices (remove-duplicates indices))))
      (destructuring-bind (&key component-index vblock-index in-vblock-index
				&allow-other-keys)
	  (fe-secondary-information fe)
	(coerce (loop for index across indices appending
		      (loop for c across component-index
			    and v across vblock-index
			    and iv across in-vblock-index
			    when (and (= c index) (zerop v))
			    collecting iv))
		'vector)))))

(defmethod local-value-blocks ((svec <sparse-vector>) (cell <cell>))
  (let ((mesh (mesh svec)))
    (map 'vector
	 (lambda (cell)
	   (let ((key (cell-key cell mesh)))
	     (and (in-pattern-p svec key)
		  (vref svec (cell-key cell mesh)))))
	 (subcells cell))))

(defmethod global-local-operation ((cell <cell>) (svec <sparse-vector>) local-vec operation)
  (let ((vblocks (local-value-blocks svec cell))
	(fe (get-fe (ansatz-space svec) cell)))
    (destructuring-bind (&key nr-dofs component-index in-component-index
			      vblock-index in-vblock-index &allow-other-keys)
	(fe-secondary-information fe)
      (let ((multiplicity (multiplicity svec)))
	(with-mutual-exclusion (svec)
	  (dotimes (i nr-dofs)
	    (let ((component-index (aref component-index i))
		  (in-component-index (aref in-component-index i))
		  (vblock-index (aref vblock-index i))
		  (in-vblock-index (aref in-vblock-index i)))
	      (dotimes (j multiplicity)
		(symbol-macrolet
		      ((local (mref (aref local-vec component-index)
				    in-component-index j))
		       (global (mref (aref vblocks vblock-index)
				     in-vblock-index j)))
		    (ecase operation
		      (:local<-global (setq local global))
		      (:global<-local (setq global local))
		      (:global+=local (incf global local))
		      (:global-=local (decf global local))))))))))))

(defmethod fill-local-from-global-vec ((cell <cell>) (svec <sparse-vector>) local-vec)
  (global-local-operation cell svec local-vec :local<-global))

(defmethod get-local-from-global-vec ((cell <cell>) (svec <sparse-vector>))
  (lret ((local-vec (make-local-vec svec cell)))
    (fill-local-from-global-vec cell svec local-vec)))

(defmethod set-global-to-local-vec ((cell <cell>) (svec <sparse-vector>) local-vec)
  (global-local-operation cell svec local-vec :global<-local))

(defmethod increment-global-by-local-vec ((cell <cell>) (svec <sparse-vector>) local-vec)
  (global-local-operation cell svec local-vec :global+=local))

(defun set-lagrange-ansatz-space-vector (asv func)
  "Sets an ansatz-space-vector to interpolate a given function.  This is
still a suboptimal implementation for vector functions."
  (doskel (cell (mesh asv) :where :surface :dimension :highest)
    (let ((local-vec (make-local-vec asv cell)))
      (let ((fe (get-fe (ansatz-space asv) cell)))
	(loop for comp from 0
	      and comp-fe across (components fe) do
	      (do-dof (dof comp-fe)
		(let ((value (funcall func (local->global cell (dof-gcoord dof)))))
		  (setq value (fl.cdr::ensure-1-component-vector value))
		  (dotimes (i (multiplicity asv))
		    (setf (mref (aref local-vec comp) (dof-index dof) i)
			  (mref (aref value comp) 0 i)))))))
      (set-global-to-local-vec cell asv local-vec))))

(defun multiple-evaluate-local-fe (local-vec shape-values)
  "Evaluates the vector given in @arg{local-vec} at multiple points.  Here
@arg{local-vec} should be a data vector obtained with
@function{get-local-from-global-vec} and @arg{ip-values} should be a vector
obtained from @function{ip-values}."
  (vector-map (lambda (point-values)
		(vector-map #'m*-tn point-values local-vec))
	      shape-values))

(defmethod extract-ip-data ((cell <cell>) qrule property-list)
  "Converts all ansatz-space objects in the parameters list into local
value arrays corresponding to the finite element."
  (loop for (key object) on property-list by #'cddr
	when (typep object '<ansatz-space-vector>)
	collect key and
	collect (multiple-evaluate-local-fe
		 (get-local-from-global-vec cell object)
		 (ip-values (get-fe (ansatz-space object) cell) qrule))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; <sparse-matrix> interface
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmethod local-value-blocks ((smat <sparse-matrix>) (cell <cell>))
  (let* ((subcell-keys (vector-map (rcurry #'cell-key (mesh smat)) (subcells cell)))
	 (n (length subcell-keys)))
    (lret ((result (make-array (twice n) :initial-element nil)))
      (dotimes (i n)
	(let ((row-key (aref subcell-keys i)))
	  (dotimes (j n)
	    (let ((col-key (aref subcell-keys j)))
	      (when (in-pattern-p smat row-key col-key)
		(setf (aref result i j)
		      (mref smat row-key col-key))))))))))

(defmethod global-local-operation ((cell <cell>) (smat <sparse-matrix>) local-mat operation)
  (declare (optimize debug))
  (let ((mblocks (local-value-blocks smat cell))
	(domain-fe (get-fe (domain-ansatz-space smat) cell))
	(image-fe (get-fe (image-ansatz-space smat) cell)))
    (destructuring-bind
	  (&key ((:nr-dofs nr-dofs-1))
		((:component-index component-index-1)) ((:in-component-index in-component-index-1))
		((:vblock-index vblock-index-1)) ((:in-vblock-index in-vblock-index-1))
		&allow-other-keys)
	(fe-secondary-information image-fe)
      (destructuring-bind
	  (&key ((:nr-dofs nr-dofs-2))
		((:component-index component-index-2)) ((:in-component-index in-component-index-2))
		((:vblock-index vblock-index-2)) ((:in-vblock-index in-vblock-index-2))
		&allow-other-keys)
	  (fe-secondary-information domain-fe)
	(dotimes (i nr-dofs-1)
	  (let ((comp-1 (aref component-index-1 i))
		(in-comp-1 (aref in-component-index-1 i))
		(vblock-index-1 (aref vblock-index-1 i))
		(in-vblock-index-1 (aref in-vblock-index-1 i)))
	    (dotimes (j nr-dofs-2)
	      (let ((comp-2 (aref component-index-2 j))
		    (in-comp-2 (aref in-component-index-2 j))
		    (vblock-index-2 (aref vblock-index-2 j))
		    (in-vblock-index-2 (aref in-vblock-index-2 j)))
		(symbol-macrolet
		      ((local (mref (mref local-mat comp-1 comp-2)
				    in-comp-1 in-comp-2))
		       (global (mref (aref mblocks vblock-index-1 vblock-index-2)
				     in-vblock-index-1 in-vblock-index-2)))
		    (ecase operation
		      (:local<-global (setq local global))
		      (:global<-local (setq global local))
		      (:global+=local (incf global local))
		      (:global-=local (decf global local))))))))))))

(defmethod fill-local-from-global-mat ((cell <cell>) (smat <sparse-matrix>) local-mat)
  (global-local-operation cell smat local-mat :local<-global))

(defmethod get-local-from-global-mat ((cell <cell>) (smat <sparse-matrix>))
  (lret ((local-mat (make-local-mat smat cell)))
    (fill-local-from-global-mat cell smat local-mat)))

(defmethod set-global-to-local-mat ((cell <cell>) (smat <sparse-matrix>) local-mat)
  (global-local-operation cell smat local-mat :global<-local))

(defmethod increment-global-by-local-mat ((cell <cell>) (smat <sparse-matrix>) local-mat)
  (global-local-operation cell smat local-mat :global+=local))


;;;; Testing
(defun test-sparseif ()
  
  )

;;; (test-sparseif)
(fl.tests:adjoin-test 'test-sparseif)
