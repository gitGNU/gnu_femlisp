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

(in-package :fl.mesh)

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
	    (if parametric
		(change-class cell (mapped-cell-class (class-of cell))
			      :mapping (funcall parametric cell))
		(change-class cell (unmapped-cell-class (class-of cell)))))))
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

(defun special-mesh-on-box-domain (domain patch->mesh-sizes)
  "Creates a uniform mesh on a box domain consisting of N_1 x ... x N_dim)
cells."
  (let ((corners->cell (make-hash-table :test 'equalp))
	(mesh (make-instance '<mesh> :domain domain)))
    (doskel (patch domain)
      (let* ((patch-dim (dimension patch))
	     (cell-class (class-of (n-cube patch-dim)))
	     (N-patch (funcall patch->mesh-sizes patch)))
	(multi-for (ivec (make-fixnum-vec patch-dim 1) N-patch)
	  (let ((corners ()))
	    (multi-for (jvec (make-fixnum-vec patch-dim -1) (make-fixnum-vec patch-dim 0))
	      (push (l2g patch
			 (map 'double-vec
			      #'(lambda (c N_i) (float (/ c N_i) 1.0))
			      (m+ ivec jvec) N-patch))
		    corners))
	    (insert-cell-from-corners mesh corners->cell cell-class
				      (nreverse corners) (list 'PATCH patch))))))
    mesh))

(defun uniform-mesh-on-box-domain (domain N)
  "Creates a uniform mesh on a box domain consisting of N_1 x ... x N_dim)
cells."
  (let ((dim (dimension domain)))
    (special-mesh-on-box-domain
     domain
     #'(lambda (patch)
	 (let ((patch-corners (corners patch)))
	   (coerce
	    (loop for i below dim
		  unless (apply #'= (mapcar (rcurry #'elt i) patch-corners))
		  collect (if (numberp N) N (aref N i)))
	    'fixnum-vec))))))

(defun isotropic-mesh-on-rectangle-domain (domain)
  "Generates a rather isotropic mesh on a domain consisting of rectangular
patches."
  (special-mesh-on-box-domain
   domain
   #'(lambda (patch)
       (if (vertex? patch)
	   #()
	   (let ((Dphi (l2Dg patch (local-coordinates-of-midpoint patch))))
	     (if (= (dimension patch) 1)
		 (loop for i from 0 below 2 for D = (mref Dphi i 0)
		       unless (zerop D) do (return (vector (truncate D))))
		 (coerce (loop for i from 0 below 2 for D = (mref Dphi i i)
			       unless (zerop D) collecting (truncate D))
			 'vector)))))))


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

(defun compare-lexicographically (&key (fuzzy 1.0d-12) direction)
  "Returns a function which compares two vectors lexicographically."
  #'(lambda (pos1 pos2)
      (loop for x across pos1 and y across pos2
	 for dirs = direction then (cdr dirs)
	 for dir = (or (car dirs) :up)
	 when (< x (- y fuzzy)) do (return (eq dir :up))
	 when (> x (+ y fuzzy)) do (return (eq dir :down))
	 finally (return nil))))

(defun sort-lexicographically (elist &key (fuzzy 1.0d-12))
  "Sorts a cell list lexicographically by the coordinates of their
midpoint."
  (sort elist (compare-lexicographically :fuzzy fuzzy) :key #'midpoint))

(defun triangulize (mesh)
  "Transforms a tensorial mesh into a simplex mesh."
  (when (typep mesh '<hierarchical-mesh>)
    (error "NYI for hierarchical meshes."))
  (when (> (dimension mesh) 2)
    (error "NYI.  Should be a Kuhn decomposition."))
  (let ((new-mesh (make-instance '<mesh> :domain (domain mesh)))
	(copy-table (make-hash-table)))
    (doskel ((cell properties) mesh :dimension :highest)
      (if (simplex-p cell)
	  (let ((new-cell (copy-cell cell)))
	    (setf (gethash cell copy-table) new-cell)
	    (setf (slot-value new-cell 'boundary)
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
;;;; Testing:
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun test-meshgen ()
  (describe (uniformly-refined-mesh *unit-interval-domain* 1))
  (describe (uniform-mesh-on-box-domain (n-cube-domain 2) #(2 2)))
  (let ((h-mesh (uniformly-refined-hierarchical-mesh *unit-interval-domain* 2)))
    (refine h-mesh)
    (loop repeat 1 do (refine h-mesh :test (rcurry #'inside-cell? #d(0.25))))
    (describe (refinement-interface h-mesh)))
  (check-identification (make-mesh-from-domain (n-cell-domain 2)))
  (describe (copy-skeleton (n-cell-domain 1)))
  (describe (refine (make-mesh-from-domain (n-cell-domain 1))))
  (let ((h-mesh (uniformly-refined-hierarchical-mesh (n-cell-domain 1) 1)))
    (loop repeat 1 do (refine h-mesh :test (rcurry #'inside-cell? #d(0.25))))
    (describe (refinement-interface h-mesh)))
  (let* ((domain (n-ball-domain 2))
	 (mesh (make-mesh-from-domain domain)))
    (describe mesh))
  )

;;; (test-meshgen)
(fl.tests:adjoin-test 'test-meshgen)
