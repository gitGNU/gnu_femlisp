;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;  mesh.lisp
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

(in-package mesh)
;;;(declaim (optimize (speed 0) (safety 3)))
;;;(declaim (optimize (debug 3)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Mesh class
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass <mesh> (<skeleton>)
  ((domain :accessor domain :initarg :domain :type <domain>)
   (parametric :accessor parametric :initform nil :initarg :parametric))
  ;;
  (:documentation "A <mesh> is a special <skeleton> mapping cells to
property lists with properties of the cell.  The most important property of
a cell is its patch in the domain.  Another one could be a list of possibly
identified cells.  The slot parametric determines which kind of cell
mappings are used for approximating the domain.  These can be the nonlinear
mappings used in the domain definition, but also arbitrary approximations,
to those mappings, e.g. isoparametric mappings.  The special value NIL
means that multilinear mappings are used for all cells outside the
boundaries."))

(defmethod initialize-instance ((mesh <mesh>) &key domain &allow-other-keys)
  "When a mesh is constructed from a domain its dimension is taken as the
domain dimension by default."
  (call-next-method)
  (assert domain)
  (unless (slot-boundp mesh 'dim)
    (setf (dimension mesh) (dimension domain))))

(defmethod make-analog ((mesh <mesh>))
  (make-instance '<mesh> :domain (domain mesh) :parametric (parametric mesh)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; cell property list handling
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(definline patch-of-cell (cell mesh)
  (get-cell-property cell mesh 'PATCH))
(definline (setf patch-of-cell) (value cell mesh)
  (setf (get-cell-property cell mesh 'PATCH) value))

(definline get-patch-property (cell mesh property)
  "Returns the value of the property."
  (getf (skel-ref (domain mesh) (patch-of-cell cell mesh)) property))
(definline (setf get-patch-property) (value cell mesh property)
  "Sets the value of the property."
  (setf (getf (skel-ref (domain mesh) (patch-of-cell cell mesh))
	      property)
	value))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Mesh refinement
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defgeneric refine (mesh &key test)
  (:documentation "Refines a mesh or hierarchical-mesh object either locally or
globally depending on the refinement criterion function 'test'.  Later: if
a region of type skeleton is provided the work is restricted to that
region."))

(defgeneric update-refinement! (mesh refined-mesh &key region test refined)
  (:documentation "Lower level interface to refinement."))


;;; Modifications for mesh cell refinement

(defmethod refine-cell! :after ((cell <cell>) (mesh <mesh>) (refined-mesh <mesh>) refined-region)
  "First, the domain patch is copied to the children.  Then, cell mappings
are modified for boundary-approximating meshes."
  ;; set patch
  (loop with patch = (patch-of-cell cell mesh)
	for child across (children cell mesh) do
	(setf (patch-of-cell child refined-mesh) patch))
  
  ;; For meshes approximating smooth boundaries, we have additionally to
  ;; modify cell mapping.  More precisely, we abandon the mappings from
  ;; interior children and create finer parametric ones for the boundary
  ;; neighbors.
  (let ((parametric (parametric mesh))
	(patch (patch-of-cell cell mesh)))
    ;; when not polygonal or directly from domain patches replace mappings
    (when (and parametric
	       (not (eq parametric :from-domain))
	       (= (dimension patch) (dimension mesh)))
      (loop for child across (children cell mesh) do
	    (unless (vertex? child)
	      (if (every #'(lambda (side) ; tests if inside some patch
			    (eq (patch-of-cell side refined-mesh) patch))
			(boundary child))
		  (change-class child (unmapped-cell-class (class-of child)))
		  (setf (mapping child) (funcall parametric child))))))))

(defmethod update-refinement! ((mesh <mesh>) (refined-mesh <mesh>)
			       &key region test refined)
  (assert (or region (not (eq mesh refined-mesh))))
  (loop	for cell being the hash-keys of
	(etable-of-highest-dim (or region mesh))
	when (and (not (refined-p cell mesh)) (funcall test cell)) do
	(refine-cell! cell mesh refined-mesh refined)))

  
(defmethod refine ((mesh <mesh>) &key (test (constantly t)))
  "Refines a mesh.  Returns two values: the first is the refined mesh,
i.e. a skeleton of the new mesh referencing the patches of the domain,
the second is the refinement, i.e. a skeleton for the old mesh
referencing the refinement vectors.

This refinement algorithm is more or less like skeleton refinement.
It usually makes sense only for global refinements.  Local refinements
should be done within hierarchical-mesh structures."
  (let* ((refined-mesh (make-analog mesh))
	 (refinement (make-instance '<skeleton> :dimension (dimension mesh))))
    (update-refinement! mesh refined-mesh :region mesh :test test :refined refinement)
    ;; decouple the refined mesh
    (doskel (cell refined-mesh)
      (remf (skel-ref refined-mesh cell) :PARENT))
    ;; return it
    (values refined-mesh refinement)))

(defmethod check :after ((mesh <mesh>))
  "Performs some additional checks for mesh."
  nil)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; h-mesh class
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass <hierarchical-mesh> (<mesh>)
  ((levels :accessor levels :initarg :levels :type (array <skeleton> (*))))
  (:documentation "Hierarchical-meshes are those meshes which will be used
most often, because they remember the refinement history and therefore
allow for refinement and coarsening.  The slot levels is an array of
skeletons containing the cells for different levels."))

(defmethod describe-object :after ((h-mesh <hierarchical-mesh>) stream)
  (loop initially (terpri stream)
	for level upto (top-level h-mesh) do
	(format stream "Level ~D~%" level)
	(describe (cells-on-level h-mesh level))))

(definline cells-on-level (mm level) (aref (levels mm) level))
(definline nr-of-levels (mm) (length (levels mm)))
(definline bottom-level-cells (mm) (cells-on-level mm 0))
(definline top-level (mm) (1- (nr-of-levels mm)))
(definline top-level-cells (mm) (vector-last (levels mm)))
(definline hierarchical-mesh-p (obj) (typep obj '<hierarchical-mesh>))
(definline flat-mesh-p (obj) (and (typep obj '<mesh>) (not (typep obj '<hierarchical-mesh>))))

(defmethod update-instance-for-different-class :after
    ((mesh <mesh>) (h-mesh <hierarchical-mesh>) &rest initargs)
  "One constructor for an h-mesh is by changing the class of a one-level
mesh.  This method definition fills the level slot appropriately."
  (declare (ignore initargs))
  (setf (slot-value h-mesh 'levels)
	(make-array 1 :element-type '<mesh>
		    :initial-element (skel-map #'(lambda (cell value)
						   (declare (ignore cell))
						   value)
					       mesh)
		    :adjustable t)))

(defmethod refine ((h-mesh <hierarchical-mesh>) &key (test (constantly t)))
  "Refine a hierarchical-mesh.  When the argument 'test' is supplied, all
cells satisfying the test are refined."
  (let ((levels (levels h-mesh))
	(dim (dimension h-mesh)))
    ;; possible extensions of refined regions
    (loop with refined-region = (make-instance '<skeleton> :dimension dim)
	  for level from 0 upto (top-level h-mesh)
	  for level-cells = (aref levels level)
	  for level+1-cells = (if (= level (top-level h-mesh))
				  (make-instance '<skeleton> :dimension dim)
				  (aref levels (1+ level)))
	  for refined = (make-instance '<skeleton> :dimension dim) do
	  (update-refinement! h-mesh h-mesh :region level-cells :test test :refined refined)
	  do
	  (unless (skel-empty-p refined)
	    (doskel ((cell children) refined)
	      ;; copy information from h-mesh to levels (should be
	      ;; non-nil, what else could go there?)
	      (setf (skel-ref level-cells cell)
		    (skel-ref h-mesh cell))
	      (loop for child across children do
		    (setf (skel-ref level+1-cells child)
			  (skel-ref h-mesh child)))
	      (setf (skel-ref refined-region cell) t)))
	  finally
	  ;; possible creation of a new refinement level
	  (unless (skel-empty-p refined)
	    (setf (levels h-mesh)
		  (adjust-array levels (1+ (length levels)) :initial-element level+1-cells)))
	  (return (values h-mesh refined-region)))))

(defun refinement-interface (h-mesh &key level)
  "Returns the refined boundary subcells of unrefined cells in a skeleton.
At the moment, this is a global operation.  Later on, it should probably be
localized."
  (make-instance
   '<skeleton> :cells
   (loop for level from (or level 0) upto (or level (1- (top-level h-mesh)))
	 nconcing
	 (when (< level (top-level h-mesh))
	   (let ((domain-boundary (domain-boundary (domain h-mesh)))
		 (level-boundary (skeleton-boundary (cells-on-level h-mesh (1+ level))))
		 (interface (make-hash-table)))
	     (doskel (cell level-boundary :dimension :highest)
	       (let ((parent (parent cell h-mesh)))
		 (when (or (identified-p cell h-mesh)
			   (not (member-of-skeleton? (patch-of-cell parent h-mesh)
						     domain-boundary)))
		   (setf (gethash parent interface) t))))
	     (hash-table-keys interface))))))

(defun hierarchically-ordered-cells (h-mesh &key level)
  "Sorts the cells up to the given level (defaulting to the last
level) hierarchically for use in something similar to the nested
disection method.  Returns a list of the sorted cells."
  (let ((elist ()))
    ;; collect the cells of the level 0 mesh in the correct order
    (doskel (cell (bottom-level-cells h-mesh)) (push cell elist))
    ;; sort the cells of the refinements
    (loop for k from 0 upto (or level (- (nr-of-levels h-mesh) 2))
	  for level-cells = (aref (levels h-mesh) k) do
	  (setq elist
		(loop for cell in elist appending
		      (loop with children = (children cell h-mesh)
			    for i from (1- (length children)) downto 0
			    collect (aref children i))
		      )))
    ;; return the resulting cell elist
    elist))

(defun for-each-cell-of-highest-dimension-on-surface (func h-mesh)
  "Calls func for each cell on the hierarchical-mesh surface."
  (loop for cell being each hash-key of (etable-of-highest-dim h-mesh)
	unless (refined-p cell h-mesh)
	do (funcall func cell)))

(defun surface-cells-of-highest-dim (h-mesh)
  "This function is needed mainly for plotting.  It returns the surface
cells of a locally refined hierarchical-mesh structure."
  (loop for cell being each hash-key of (etable-of-highest-dim h-mesh)
	unless (refined-p cell h-mesh)
	collect cell))

(defun nr-of-surface-cells (h-mesh)
  (let ((sum 0))
    (doskel (cell h-mesh :dimension :highest :where :surface)
      (incf sum))
    sum))

(defmethod find-cell-from-position ((h-mesh <hierarchical-mesh>) (pos array))
  "Hierarchical search for a leaf cell containing the given position."
  (loop for level from 0 below (nr-of-levels h-mesh)
	for cell = (find-cell-from-position (cells-on-level h-mesh level) pos)
	unless (null cell) do
	(return
	  (loop for children = (children cell h-mesh) while children do
		(let ((child (loop for child across children
				   when (and (= (dimension child) (dimension h-mesh))
					     (inside-cell? child pos))
				   do (return child))))
		  (if child
		      (setq cell child)
		      (return nil))) ; position not covered on higher level
		finally (return cell)))))

(defmethod check :after ((h-mesh <hierarchical-mesh>))
  "Performs some additional checks for hierarchical meshes."
  (loop for level-mesh across (levels h-mesh) do
	(doskel (cell level-mesh)
	  ;; check properties
	  (let ((props-1 (skel-ref level-mesh cell))
		(props-2 (skel-ref h-mesh cell)))
	    (unless (eq props-1 props-2)
	      (error "~&~A: Level-Properties do not agree:~%~A~%~A~%"
		     cell props-1 props-2)))
	  ;; check refinements
	  (whereas ((children (children cell h-mesh)))
	    (loop for child across children do
		  (assert (eq (parent child h-mesh) cell))
		  (let ((patch-1 (patch-of-cell child h-mesh))
			(patch-2 (patch-of-cell cell h-mesh)))
		    (unless (eq patch-1 patch-2)
		      (error "~&~A: Patch of parent ~A different:~%~A~%~A~%"
			     child cell patch-1 patch-2)))))
	  ;; check patches
	  )))

;;;; Testing
(defun test-mesh ()
  "More tests can be found in meshgen.lisp."
  (let ((*print-skeleton-values* t))
    (describe
     (make-instance '<mesh> :domain *unit-interval-domain*))))

(tests::adjoin-femlisp-test 'test-mesh)
