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

(in-package :fl.mesh)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Mesh class
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass <mesh> (<skeleton>)
  ((domain :reader domain :initarg :domain :type <domain>
	   :documentation "The domain of the mesh.")
   (parametric :reader parametric :initform nil :initarg :parametric))
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
  (unless (slot-boundp mesh 'dimension)
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

(defmethod dimension-of-part ((mesh <mesh>) part)
  "Parts of a skeleton can be named with the property @symbol{:part}."
  (dimension-of-part (domain mesh) part))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; cell refinement
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmethod refine-cell! :after (rule (cell <cell>) (mesh <mesh>)
				(refined-mesh <mesh>) refined-region)
  "This :AFTER-method for cell refinement copies the domain patch to the
children.  For meshes approximating smooth boundaries, the mappings from
interior children are abandoned and finer parametric ones for the boundary
neighbors are generated."
  (declare (ignore rule refined-region))
  ;; set patch
  (loop with patch = (patch-of-cell cell mesh)
	for child across (children cell mesh) do
	(setf (patch-of-cell child refined-mesh) patch))
  (let ((parametric (parametric mesh))
	(patch (patch-of-cell cell mesh)))
    ;; when not polygonal or directly from domain patches replace mappings
    (when (and parametric
	       (not (eq parametric :from-domain))
	       #-(or)(= (dimension patch) (dimension mesh))
	       #+(or)(plusp (dimension cell)))
      (loop for child across (children cell mesh) do
	    (unless (vertex? child)
	      (if (every #'(lambda (side) ; tests if inside some patch
			    (eq (patch-of-cell side refined-mesh) patch))
			(boundary child))
		  (change-class child (unmapped-cell-class (class-of child)))
		  (change-class child (mapped-cell-class (class-of child))
				:mapping (funcall parametric child))))))))

;;; (untrace :methods 'refine-cell!)
(defmethod check progn ((mesh <mesh>) &key &allow-other-keys)
  "Performs some additional checks for mesh."
  nil)

(defgeneric meshsize (mesh)
  (:documentation "Computes a meshsize.  Please refer to the method
documentations for the exact definition.")
  (:method ((mesh <mesh>))
    "Computes a meshsize as the size of the longest edge in the mesh."
    (let ((size 0.0))
      (doskel (cell mesh :where :surface :dimension 1)
	(let ((diam (diameter cell)))
	  (when (> diam size)
	    (setf size diam))))
      size)))

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

(defun cells-on-level (mm level) (aref (levels mm) level))
(defun nr-of-levels (mm) (length (levels mm)))
(defun bottom-level-cells (mm) (cells-on-level mm 0))

(defgeneric top-level (mh)
  (:documentation "Top-level of a mesh hierarchy.")
  (:method ((mm <hierarchical-mesh>))
      (1- (nr-of-levels mm))))

(defun top-level-cells (mm) (vector-last (levels mm)))
(defun hierarchical-mesh-p (obj) (typep obj '<hierarchical-mesh>))
(defun flat-mesh-p (obj) (and (typep obj '<mesh>) (not (typep obj '<hierarchical-mesh>))))

(defun level-of-cell (cell h-mesh)
  "Returns the level of @arg{cell} in the hirearchical mesh @arg{h-mesh}."
  (loop for level from 0 upto (top-level h-mesh)
     when (member-of-skeleton? cell (cells-on-level h-mesh level))
     do (return level)))

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

(defmethod refine ((h-mesh <hierarchical-mesh>) &key (indicator (constantly t)))
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
	  (update-refinement! h-mesh h-mesh :region level-cells
			      :indicator indicator :refined refined)
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
Those cells are found as all refined cells which are not part of the domain
boundary.  At the moment, this is a global operation.  Later on, it should
probably be localized."
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
			   ;; in the following (patch-of-cell parent
			   ;; h-mesh) does not work for porous domains (see
			   ;; hom-mg.lisp)
			   (not (member-of-skeleton? (patch-of-cell cell h-mesh)
						     domain-boundary)))
		   (setf (gethash parent interface) t))))
	     (hash-table-keys interface))))))

(defun hierarchically-ordered-cells (h-mesh &key (where :surface) level)
  "Sorts the cells up to the given level (defaulting to the last
level) hierarchically for use in something similar to the nested
disection method.  Returns a list of the sorted cells."
  (let ((table (make-hash-table))
        (result ()))
    (labels ((collect (cell current-level dim)
               (if (and (or (null level) (= current-level level))
                        (or (null where)
                            (ecase where
                              (:surface (not (refined-p cell h-mesh)))
                              (:refined (refined-p cell h-mesh))
                              (:all t))))
                   (loop+ ((subcell (subcells cell))) do
                        (unless (gethash subcell table)
                          (push (setf (gethash subcell table) subcell) result)))
                   (when (or (null level) (< current-level level))
                     (for-each (rcurry #'collect (1+ current-level) dim)
                               (children cell h-mesh))))))
    (doskel (cell (bottom-level-cells h-mesh) :dimension :highest)
      (let* ((patch (patch-of-cell cell h-mesh))
             (dim (dimension patch)))
        (when (and (member :substance (patch-classification patch (domain h-mesh)))
                   (= (dimension cell) dim))
          (collect cell 0 dim))))
    (nreverse result))))

(defun for-each-cell-of-highest-dimension-on-surface (func h-mesh)
  "Calls func for each cell on the hierarchical-mesh surface."
  (skel-for-each func h-mesh :dimension :highest :where :surface))

(defun surface-cells-of-dim (h-mesh dim)
  "This function returns the surface cells of a locally refined
hierarchical-mesh structure."
  (mapper-collect #'skel-for-each h-mesh :dimension dim :where :surface))

(defun surface-cells-of-highest-dim (h-mesh)
  "This function returns the surface cells of highest dimension of a
locally refined hierarchical-mesh structure."
  (surface-cells-of-dim h-mesh (dimension h-mesh)))

(defun nr-of-surface-cells (h-mesh)
  (mapper-count #'skel-for-each h-mesh :dimension :highest :where :surface))

(defun hierarchical-search (h-mesh test &key level)
  "Hierarchical search for a subtree of cells in @arg{h-mesh} satisfying
@arg{test}.  A leaf cell is returned, if successful, otherwise NIL."
  (labels ((h-mesh-search (cell current-level)
             (and (funcall test cell current-level)
                  (let ((children (children cell h-mesh)))
                    (if (and level (= current-level level))
                        cell
                        (if children
                            (loop for child across children
                               when (h-mesh-search child (1+ current-level))
                               return it)
                            (unless level cell)))))))
    (loop for level from 0 upto (or level (top-level h-mesh))
         thereis
         (doskel (cell (cells-on-level h-mesh level))
           (awhen (h-mesh-search cell level)
             (return it))))))

(defun find-cell-on-patch-with-midpoint (h-mesh patch midpoint level)
  (hierarchical-search
   h-mesh (lambda (cell cell-level)
            (and (eq (patch-of-cell cell h-mesh) patch)
                 (if (< cell-level level)
                     (inside-cell? cell midpoint 1.0e-10)
                     (equalp (midpoint cell) midpoint))))))

(defvar *find-cell-base-level* 0
  "First level on which a cell is searched for containing a given position.")

(defmethod find-cell-from-position ((h-mesh <hierarchical-mesh>) (pos array))
  "Hierarchical search for a leaf cell containing the given position.  A result
of NIL is given if no cell covering @arg{pos} is found."
  (loop for level from *find-cell-base-level* below (nr-of-levels h-mesh)
     for cell = (find-cell-from-position
                 (cells-on-level h-mesh level) pos)
     when cell do
       (return
         (loop for children = (children cell h-mesh) while children do
              (let ((child (loop for child across children
                              when (and (= (dimension child) (dimension h-mesh))
                                        (inside-cell? child pos))
                              do (return child))))
                (if child
                    (setq cell child)
                    ;; In the following case pos is not covered by cell's
                    ;; children.  This may happen for non-polygonal
                    ;; domains, but unfortunately also due to rounding
                    ;; errors.
                    (return nil)))
              finally (return cell)))))

(defvar *allow-child-patch-change* nil
  "When T it is allowed for a child to have a different patch from its
father.  This is used at the moment for homogenization of porous domains,
but it is not a standard situation.")

(defmethod check progn ((h-mesh <hierarchical-mesh>) &key &allow-other-keys)
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
		      (unless *allow-child-patch-change*
			(error "~&~A: Patch of parent ~A different:~%~A~%~A~%"
			       child cell patch-1 patch-2))))))
	  ;; check patches
	  )))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Tests
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(in-suite mesh-suite)

(test mesh
  (is (= 16641 (length
                (hierarchically-ordered-cells
                 (uniformly-refined-hierarchical-mesh (n-cube-domain 2) 6)))))
  (finishes
    (let ((*print-skeleton-values* t))
      (describe
       (make-instance '<mesh> :domain (n-simplex-domain 1)))))
  )

