;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;  meshgen.lisp
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

;;;; meshgen.lisp contains routines for mesh and hierarchical-mesh
;;;; generation.

(defgeneric make-mesh-from-domain (domain &key &allow-other-keys)
  (:documentation "Transforms a domain which is specified sufficiently well
into a mesh."))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Meshes
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmethod make-mesh-from-domain ((domain <domain>) &key parametric &allow-other-keys)
  "This function creates a mesh on the domain skeleton.  The mesh either
uses the nonlinear mappings of the domain patches (for parametric =
:from-domain) or approximates those with polygonal (parametric = nil) or
isoparametric mappings.  Note that the boundary always uses the domain
mappings."
  (let ((domain-dim (dimension domain)))
    (multiple-value-bind (mesh patch->cell)
	(copy-skeleton domain)
      ;; post-processing
      (dohash ((patch cell) patch->cell)
	;; set patch
	(setf (get-cell-property cell mesh 'PATCH) patch)
	(unless (eq parametric :from-domain)
	  ;; drop possible mappings for all cells of highest dimension
	  (unless (or (vertex? cell) (< (dimension cell) domain-dim))
	    (setf (mapping cell)
		  (and parametric
		       (funcall parametric cell))))))
      (change-class mesh '<mesh> :domain domain :parametric parametric))))

(defun copy-mesh (mesh)
  "Copies a mesh.  Properties copied are only patch and identification.  If
necessary, one might add further properties to be copied as a keyword
argument."
  (let ((new-mesh (make-analog mesh)))
    (doskel (cell mesh)
      (setf (patch-of-cell cell new-mesh)
	    (patch-of-cell cell mesh))
      (setf (cell-identification cell new-mesh)
	    (cell-identification cell mesh)))
    ;; return result
    new-mesh))



(defun uniformly-refined-mesh (domain n &key parametric)
  "Generates a mesh by refining the domain partition uniformly."
  (do ((mesh (make-mesh-from-domain domain :parametric parametric)
	     (refine mesh))
       (k 0 (1+ k)))
      ((= k n) mesh)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Mesh construction in the UG way
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun insert-cell-from-corners (mesh corners->cell cell-class corners properties)
  "Creates a cell of type cell-class with corners given by corners.
corners->cell has to be an equalp hash-table mapping corners to the
corresponding cell.  It is updated by this function."
  (assert (flat-mesh-p mesh))
  (or (gethash corners corners->cell)
      (let* ((refcell (cell-class-reference-cell cell-class))
	     (refcell-corners (corners refcell)))
	(assert (= (length refcell-corners) (length corners)))
	(flet ((refcell-corner->corner (refcell-corner)
		 (nth (position refcell-corner refcell-corners :test #'equal) corners)))
	  (let ((cell (if (vertex? refcell)
			  (make-vertex (car corners))
			  (copy-cell refcell))))
	    (setf (boundary cell)
		  (vector-map
		   #'(lambda (refcell-side)
		       (let ((side-corners (mapcar #'refcell-corner->corner
						   (corners refcell-side))))
			 (insert-cell-from-corners
			  mesh corners->cell
			  (cell-class refcell-side) side-corners properties)))
		   (boundary refcell)))
	    (setf (skel-ref mesh cell) properties)
	    (setf (gethash corners corners->cell) cell))))))

(defun uniform-mesh-on-box-domain (domain N)
  "Creates a uniform mesh on a box domain consisting of N_1 x ... x N_dim)
cells."
  (let* ((dim (dimension domain))
	 (corners->cell (make-hash-table :test 'equalp))
	 (mesh (make-instance '<mesh> :domain domain))
	 (N (if (numberp N) (make-fixnum-vec dim N) N)))
    (doskel (patch domain)
      (let* ((patch-dim (dimension patch))
	     (cell-class (cell-class (n-cube patch-dim)))
	     (patch-corners (corners patch))
	     (N-patch
	      (coerce
	       (loop for i below dim
		     unless (apply #'= (mapcar (rcurry #'elt i) patch-corners))
		     collect (elt N i))
	       'fixnum-vec)))
	(multi-for (ivec (make-fixnum-vec patch-dim 1) N-patch)
	  (let ((corners ()))
	    (multi-for (jvec (make-fixnum-vec patch-dim -1) (make-fixnum-vec patch-dim 0))
	      (push (l2g patch
			 (map 'double-vec
			      #'(lambda (c N_i) (float (/ c N_i) 1.0d0))
			      (ivec+ ivec jvec) N-patch))
		    corners))
	    (insert-cell-from-corners mesh corners->cell cell-class
				      (nreverse corners) (list 'PATCH patch))))))
    mesh))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Hierarchical meshes
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defgeneric make-hierarchical-mesh-from-domain (domain &key &allow-other-keys)
  (:documentation "Construct a hierarchical-mesh from a domain."))

(defmethod make-hierarchical-mesh-from-domain (domain &key parametric &allow-other-keys)
  "Construct a hierarchical-mesh from a domain."
  (let ((mesh (make-mesh-from-domain domain :parametric parametric)))
    (change-class mesh '<hierarchical-mesh>)))

(defun uniformly-refined-hierarchical-mesh (domain level &key parametric)
  (let ((mm (make-hierarchical-mesh-from-domain domain :parametric parametric)))
    (dotimes (k level mm)
      (refine mm))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Miscellaneous
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun sort-lexicographically (elist &key (fuzzy 1.0d-12))
  "Sorts a cell list lexicographically by the coordinates of their
midpoint."
  (flet ((my-< (pos1 pos2)
	   (loop for x across pos1 and y across pos2
		 when (< x (- y fuzzy)) do (return t)
		 when (> x (+ y fuzzy)) do (return nil)
		 finally (return nil))))
    (sort elist #'my-< :key #'midpoint)))

(defun triangulize (mesh)
  "Transforms a tensorial mesh into a simplex mesh."
  (when (typep mesh '<hierarchical-mesh>)
    (error "NYI for hierarchical meshes."))
  (when (> (dimension mesh) 2)
    (error "NYI.  Should be a Kuhn decomposition."))
  (let ((new-mesh (make-instance '<mesh> :domain (domain mesh)))
	(copy-table (make-hash-table)))
    (doskel ((cell properties) mesh)
      (if (simplex? cell)
	  (let ((new-cell (copy-structure cell)))
	    (setf (gethash cell copy-table) new-cell)
	    (setf (boundary new-cell)
		  (vector-map (rcurry #'gethash copy-table)
			      (boundary cell)))
	    (setf (skel-ref new-mesh new-cell) properties))
	  ;; must be a quadrangle, triangulate it
	  (let ((subcells (subcells cell)))
	    (let ((n00 (aref subcells 8))
		  (n11 (aref subcells 5))
		  (e00-10 (aref subcells 4))
		  (e01-11 (aref subcells 3))
		  (e00-01 (aref subcells 2))
		  (e10-11 (aref subcells 1)))
	      (let ((diagonal (make-line (gethash n00 copy-table)
					 (gethash n11 copy-table))))
		(setf (skel-ref new-mesh diagonal) properties)
		(let ((triangle1 (make-simplex
				  (vector (gethash e10-11 copy-table)
					  diagonal (gethash e00-10 copy-table))))
		      (triangle2 (make-simplex
				  (vector (gethash e01-11 copy-table)
					  diagonal (gethash e00-01 copy-table)))))
		  (setf (skel-ref new-mesh triangle1) properties)
		  (setf (skel-ref new-mesh triangle2) properties)))))))
    ;; return new triangle mesh
    new-mesh))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Testing: (test-meshgen)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun test-meshgen ()
  (describe (uniformly-refined-mesh *unit-interval-domain* 1))
  (let ((h-mesh (uniformly-refined-hierarchical-mesh *unit-interval-domain* 2)))
    (refine h-mesh)
    (loop repeat 1 do (refine h-mesh :test (rcurry #'inside-cell? #(0.25))))
    (describe (refinement-interface h-mesh)))
  (check-identification (make-mesh-from-domain (n-cell-domain 2)))
  (describe (copy-skeleton (n-cell-domain 1)))
  (describe (refine (make-mesh-from-domain (n-cell-domain 1))))
  (let ((h-mesh (uniformly-refined-hierarchical-mesh (n-cell-domain 1) 1)))
    (loop repeat 1 do (refine h-mesh :test (rcurry #'inside-cell? #(0.25))))
    (describe (refinement-interface h-mesh)))
  )
(tests::adjoin-femlisp-test 'test-meshgen)

