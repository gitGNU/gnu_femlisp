;;; -*- mode: lisp; fill-column: 70 -*-

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

(defgeneric make-mesh-from (object &key &allow-other-keys)
  (:documentation "Makes a mesh from or for the given object."))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Meshes
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmethod make-mesh-from ((domain <domain>)
                           &key parametric (initial-mesh-refinements 0)
                           &allow-other-keys)
  "This function creates a mesh on the domain skeleton.

The mesh either uses the nonlinear mappings of the domain patches (for
parametric = :from-domain) or approximates those with
polygonal (parametric = nil) or isoparametric mappings.  Note that the
boundary always uses the domain mappings.

The parameter @arg{base-level} can be used when starting from a
refined mesh as base-level is desired.  This can be used when the mesh
derived from the domin definition is too coarse which may be the case
especially when it is used in the context of a domain decomposition."
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
      (lret ((mesh (change-class mesh '<mesh> :domain domain :parametric parametric)))
        ;; start from refined mesh if desired
        (loop repeat initial-mesh-refinements
              do (setf mesh (refine mesh)))
        ))))

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
  (do ((mesh (make-mesh-from domain :parametric parametric)
	     (refine mesh))
       (k 0 (1+ k)))
      ((= k n) mesh)))

(defun special-mesh-on-box-domain (domain patch->mesh-sizes)
  "Creates a uniform mesh consisting of N_1 x ... x N_dim cells on a box
domain."
  (let ((corners->cell (make-hash-table :test 'equalp))
	(mesh (make-instance '<mesh> :domain domain)))
    (doskel (patch domain)
      (let* ((patch-dim (dimension patch))
	     (cell-class (class-of (n-cube patch-dim)))
	     (N-patch (funcall patch->mesh-sizes patch)))
	(multi-for (ivec (make-fixnum-vec patch-dim 1) N-patch)
	  (let ((corners ()))
	    (multi-for (jvec (make-fixnum-vec patch-dim -1)
                             (make-fixnum-vec patch-dim 0)
                             :from-end t)
	      (push (l2g patch
			 (map 'double-vec
			      #'(lambda (c N_i) (float (/ c N_i) 1.0))
                              (m+ ivec jvec) N-patch))
		    corners))
            (dbg-show :mesh (reverse corners))
	    (insert-cell-from-corners mesh corners->cell cell-class
				      (nreverse corners)
                                      (list 'PATCH patch))))))
    mesh))

(defun uniform-mesh-on-box-domain (domain N)
  "Creates a uniform mesh consisting of N_1 x ... x N_dim cells on a box
domain."
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

(defmethod make-hierarchical-mesh-from-domain (domain &rest key-args &key &allow-other-keys)
  "Construct a hierarchical-mesh from a domain."
  (let ((mesh (apply #'make-mesh-from domain key-args)))
    (change-class mesh '<hierarchical-mesh>)))

(defun make-hierarchical-mesh-from (&rest args)
  "Creates a hierarchical-mesh from the given arguments.  See @function{MAKE-MESH-FROM}."
  (change-class (apply #'make-mesh-from args) '<hierarchical-mesh>))

(defmethod make-hierarchical-mesh-from-domain
    (domain &rest key-args &key &allow-other-keys)
  "Construct a hierarchical-mesh from a domain."
  (let ((mesh (apply #'make-mesh-from domain key-args)))
    ;; and generate the hierarchical mesh
    (change-class mesh '<hierarchical-mesh>)))

(defun uniformly-refined-hierarchical-mesh (domain level &rest key-args &key &allow-other-keys)
  (let ((mm (apply #'make-hierarchical-mesh-from-domain domain key-args)))
    (dotimes (k level mm)
      (refine mm))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Miscellaneous
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

#+(or)
(defun compare-lexicographically (&key (fuzzy 1.0d-12) direction)
  "Returns a function which compares two vectors lexicographically."
  #'(lambda (pos1 pos2)
      (loop for x across pos1 and y across pos2
	 for dirs = direction then (cdr dirs)
	 for dir = (or (car dirs) :up)
	 when (< x (- y fuzzy)) do (return (eq dir :up))
	 when (> x (+ y fuzzy)) do (return (eq dir :down))
	 finally (return nil))))

(defun compare-lexicographically (&key (fuzzy 1.0d-12) direction)
  "Returns a function which compares two vectors lexicographically."
  #'(lambda (pos1 pos2)
      (block nil
        (let ((dirs direction))
          (map nil (lambda (x y)
                     (let ((dir (or (car dirs) :up)))
                       (when (< x (- y fuzzy))
                         (return (eq dir :up)))
                       (when (> x (+ y fuzzy))
                         (return (eq dir :down)))
                       (setq dirs (cdr dirs))))
               pos1 pos2)))))

(defun sort-lexicographically (elist &key (fuzzy 1.0d-12))
  "Sorts a cell list lexicographically by the coordinates of their
midpoint."
  (sort elist (compare-lexicographically :fuzzy fuzzy) :key #'midpoint))

(defun triangulize (mesh)
  "Transforms a product-cell mesh into a simplex mesh."
  (declare (optimize debug))
  (when (typep mesh '<hierarchical-mesh>)
    (error "NYI for hierarchical meshes."))
  (when (> (dimension mesh) 2)
    (error "NYI.  Should be a Kuhn decomposition."))
  (lret ((new-mesh (make-instance '<mesh> :domain (domain mesh))))
    (doskel ((cell properties) mesh)
      (if (or (< (dimension cell) 2)
              (simplex-p cell))
          (setf (skel-ref new-mesh cell) properties)
	  ;; must be a quadrangle, triangulate it
	  (let ((subcells (subcells cell)))
	    (let ((n00 (aref subcells 8))
		  (n11 (aref subcells 5))
		  (e00-10 (aref subcells 4))
		  (e01-11 (aref subcells 3))
		  (e00-01 (aref subcells 2))
		  (e10-11 (aref subcells 1)))
	      (let ((diagonal (make-line n00 n11)))
                (loop for new-cell in
                     (list diagonal
                           (make-simplex (vector e10-11 diagonal e00-10))
                           (make-simplex (vector e01-11 diagonal e00-01)))
                     do
                     (setf (skel-ref new-mesh new-cell) properties)))))))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Tests
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(in-suite mesh-suite)

(test meshgen
  (is (= 8 (length (corners (n-cube 3)))))
  (finishes
    (uniform-mesh-on-box-domain (n-cube-domain 3) #(1 1 1))
    (triangulize (make-mesh-from (n-cube-domain 2)))
    (describe (uniformly-refined-mesh (n-simplex-domain 1) 1))
    (describe (uniform-mesh-on-box-domain (n-cube-domain 2) #(2 2)))
    (let ((h-mesh (uniformly-refined-hierarchical-mesh (n-simplex-domain 1) 2)))
      (refine h-mesh)
      (loop repeat 1 do (refine h-mesh :indicator (rcurry #'inside-cell? #d(0.25))))
      (describe (refinement-interface h-mesh)))
    (check-identification (make-mesh-from (n-cell-domain 2)))
    (describe (copy-skeleton (n-cell-domain 1)))
    (describe (refine (make-mesh-from (n-cell-domain 1))))
    (let ((h-mesh (uniformly-refined-hierarchical-mesh (n-cell-domain 1) 1)))
      (loop repeat 1 do (refine h-mesh :indicator (rcurry #'inside-cell? #d(0.25))))
      (describe (refinement-interface h-mesh)))
    (let* ((domain (n-ball-domain 2))
           (mesh (make-mesh-from domain)))
      (describe mesh))
    )
  )
