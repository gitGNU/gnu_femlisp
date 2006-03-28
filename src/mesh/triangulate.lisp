;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;  triangulate.lisp - basic triangulation stuff 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;
;;; Copyright (C) 2004 Nicolas Neuss, University of Heidelberg.
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
;;; Boundary cells
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass <boundary-cell> (<standard-cell>)
  ((dimension :reader dimension :initarg :dimension)
   (midpoint :reader midpoint :initform nil :initarg :midpoint)
   (holes :reader holes :initform () :initarg :holes))
  (:documentation "This cell is only fuzzily defined.  Its use is mostly
for defining domains by their boundary.  The slot @slot{midpoint} can be
useful for the graphical output of the cell, the slot @slot{holes} contains
a list of points lying inside holes.  This is intended as help for
triangulation programs."))

(defmethod subcells ((cell <boundary-cell>))
  (let ((skel (skeleton (boundary cell))))
    (coerce (cons cell (find-cells (constantly t) skel)) 'vector)))

(defmethod copy-cell ((cell <boundary-cell>))
  (let ((copy (call-next-method)))
    (setf (slot-value copy 'dimension) (dimension cell))
    (setf (slot-value copy 'midpoint) (midpoint cell))
    (setf (slot-value copy 'holes) (holes cell))
    copy))

(defmethod transform-cell! :after ((cell <boundary-cell>) transformation)
  (with-slots (midpoint holes)
    cell
    (setf midpoint (evaluate transformation midpoint))
    (setf holes (mapcar (curry #'evaluate transformation) holes))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Triangulation interface
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun triangulate (domain &rest args &key parametric &allow-other-keys)
  "Triangulate @arg{domain} by successively building a mesh on the domain
skeleton starting from the 0-dimensional patches."
  (assert (not (eq parametric :from-domain)))
  (let ((mesh (make-instance '<mesh> :domain domain :parametric parametric)))
    (loop for dim upto (dimension domain) do
	  (setq args (apply #'extend-triangulation mesh dim args)))
    (check mesh)
    mesh))

(defgeneric extend-triangulation (mesh dim &key &allow-other-keys)
  (:documentation "Extend an existing triangulation on the
@arg{dim}-1-skeleton to the @arg{dim}-skeleton."))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Triangulation of 0-dimensional patches
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmethod extend-triangulation (mesh (dim (eql 0)) &rest args)
  (let* ((domain (domain mesh))
	 (patches (cells-of-dim domain 0)))
    (loop for patch = (first patches) while patch do
     (let* ((identified-patches (identified-cells patch domain))
	    (identified-cells (mapcar #'(lambda (patch)
					  (make-vertex (vertex-position patch)))
				      identified-patches)))
       (loop for cell in identified-cells
	     and patch in identified-patches do
	     (setf (patch-of-cell cell mesh) patch))
       (identify identified-cells mesh)
       (setq patches (nset-difference patches identified-patches)))))
  ;; pass on args
  args)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Triangulation of 1-dimensional patches
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun compute-raster (patch &key (threshold 0.08) meshsize indicator
		       patch->raster (minimal-depth 0) (maximal-depth nil)
		       &allow-other-keys)
  "Computes a curvature-adapted raster on a curved line.  @arg{threshold}
is controlling the allowed deviation from linearity, @arg{meshsize} is a
maximally allowed mesh width.  The user may also pass an @arg{indicator}
which should be a function of patch, local position and meshwidth, or a
function @arg{patch->raster} mapping patches to complete rasters."
  (awhen patch->raster  ; user-supplied raster
    (return-from compute-raster
      (let ((result (funcall it patch)))
	(if (listp result)
	    result
	    (loop for i from 1 below result collect (float (/ i result) 1.0))))))
  ;; compute a raster for a differentiable curve
  (labels ((curve (x) (local->global patch (double-vec x)))
	   (raster (alpha x beta y depth)
	     (let* ((gamma (/ (+ alpha beta) 2))
		    (z (curve gamma))
		    (delta (norm (m- z (scal 0.5 (m+ x y)))))
		    (h (norm (m- x y)))
		    (eps (* threshold h)))
	       (when (block nil
		       (when (< depth minimal-depth)
			 (return t))
		       (when (and maximal-depth (> depth maximal-depth))
			 (return nil))
		       (whereas ((result (and indicator (funcall indicator patch z h))))
			 (ecase result
			   (:yes (return t))
			   (:no (return nil))))
		       (awhen meshsize
			 (when (> h it) (return t))
			 (when (and (< h (* 0.1 it))
				    (plusp depth))   ; beware of circular segment
			   (return nil)))
		       (> delta eps))
		 (append (raster alpha x gamma z (1+ depth))
			 (list gamma)
			 (raster gamma z beta y (1+ depth)))))))
    (raster 0.0 (curve 0.0) 1.0 (curve 1.0) 0)))

(defun fix-identification-for-positive-dimension (mesh dim)
  "Copies domain identification to cells of dimension larger 0."
  (assert (plusp dim))
  (let ((domain (domain mesh))
	(table (make-hash-table :test 'equalp)))
    (flet ((key (cell)
	     (list (cell-identification (patch-of-cell cell mesh) domain)
		   (map 'list #'(lambda (side)
				  (cell-identification side mesh))
			(boundary cell)))))
      (doskel (cell mesh :dimension dim)
	(when (identified-p (patch-of-cell cell mesh) domain)
	  (push cell (gethash (key cell) table))))
      (loop for cluster being each hash-value of table
	   for id = (make-instance 'identification :skeleton mesh :cells cluster) do
	   (dolist (cell cluster)
	     (setf (cell-identification cell mesh) id))))))

(defmethod extend-triangulation (mesh (dim (eql 1)) &rest keys)
  "Generates a mesh on the 1-skeleton for domain.  Identification is taken
care of."
  (declare (optimize debug))
  (let ((patch->nodes (make-hash-table)) ; table of (s . node(s)) entries
	(vertex-table (make-hash-table))
	(domain (domain mesh)))
    ;; initialize tables of patch vertices
    (doskel (patch domain :dimension 1)
      (push (cons 0.0 (find-cell
		       #'(lambda (cell)
			   (eq (patch-of-cell cell mesh) (aref (boundary patch) 1)))
		       mesh))
	    (gethash patch patch->nodes)))
    ;; generate internal vertices on patches of dimension 1
    (doskel (patch domain :dimension 1)
      (let ((identified-patches (identified-cells patch domain)))
	(when (eq patch (car identified-patches))
	  (let ((raster (apply #'compute-raster patch keys)))
	    (dolist (s raster)
	      (let ((identified-cells
		     (mapcar #'(lambda (patch)
				 (make-vertex (local->global patch (double-vec s))))
			     identified-patches)))
		(loop for cell in identified-cells
		      and patch in identified-patches do
		      (setf (patch-of-cell cell mesh) patch)
		      (unless (single? identified-cells)
			(identify identified-cells mesh))
		      (push (cons s cell) (gethash patch patch->nodes)))))))))
    ;; finalize tables of patch vertices
    (doskel (patch domain :dimension 1)
      (push (cons 1.0 (find-cell
		       #'(lambda (cell)
			   (eq (patch-of-cell cell mesh) (aref (boundary patch) 0)))
		       mesh))
	    (gethash patch patch->nodes)))
    ;; number all those vertices
    (let ((count 0))
      (doskel (vertex mesh :dimension 0)
	(setf (gethash vertex vertex-table) (incf count))
	(dbg :triangulate "Node ~D: ~A"  count (vertex-position vertex))))
    ;; generate the edges
    (dolist (patch (cells-of-dim (domain mesh) 1))
      (loop for (a b) on (gethash patch patch->nodes) while b do
	    (when (> (gethash (cdr a) vertex-table) (gethash (cdr b) vertex-table))
	      (rotatef a b))
	    (let ((mapping
		   (aand (mapping patch)
			 (transform-function
			  it :domain-transform
			  (list (make-real-matrix `(,(- (car b) (car a))))
				(double-vec (car a)))))))
	    (setf (skel-ref mesh (make-line (cdr a) (cdr b) :mapping mapping))
		  (list 'PATCH patch))
	    (dbg :triangulate "Edge : ~A" (list a b)))))
    (fix-identification-for-positive-dimension mesh 1)
    ;; we return the vertex-table for further processing.  This might be
    ;; removed when indices are incorporated into the mesh data structure.
    (list* :vertex-table vertex-table keys)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Triangulation of 2-dimensional patches
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defgeneric triangulate-2d (mesh program &key &allow-other-keys)
  (:documentation "Extends the mesh on the 1-skeleton to the 2-skeleton
using @arg{program}."))

(defmethod extend-triangulation (mesh (dim (eql 2)) &rest keys
				 &key (program :triangle) vertex-table
				 &allow-other-keys)
  "Calls the generic @function{triangulate-2d} which dispatches on the
triangulation program."
  (triangulate-2d mesh program :vertex-table vertex-table)
  keys)

(defmethod extend-triangulation :after (mesh (dim (eql 2)) &key &allow-other-keys)
  "Does the setup of isoparametric mapping of boundary cells."
  (when (parametric mesh)
    (doskel (cell mesh :dimension 2)
      (when (some #'(lambda (side)
		      (and (not (eq (patch-of-cell cell mesh)
				    (patch-of-cell side mesh)))
			   (mapped-p (patch-of-cell side mesh))))
		  (boundary cell))
	(change-class cell (mapped-cell-class (class-of cell))
		      :mapping (funcall (parametric mesh) cell))))))

;;; Testing

(defun test-triangulate ()
  (check (triangulate (n-cube-domain 0)))
  (check (triangulate (n-cube-domain 1)))
  (check (triangulate (n-cell-domain 1)))
  
  (let* ((vtx (make-vertex #d(1.0 0.0)))
	 (circle (make-line vtx vtx :mapping (circle-function 1.0 #d(0.0 0.0) (* 2 pi))))
	 (ball (make-instance '<boundary-cell>
			      :dimension 2
			      :boundary (vector circle)
			      :midpoint #d(0.0 0.0)))
	 (domain (change-class (skeleton ball) '<domain>))
	 (mesh (triangulate domain :meshsize 0.3)))
    #+(or)(fl.plot:plot mesh)
    (check domain)
    (check mesh)
    (assert				; test if boundary mappings are there
     (every #'mapped-p
	    (find-cells #'(lambda (cell)
			    (and (= (dimension cell) 1)
				 (not (eq (patch-of-cell cell mesh) ball))))
			mesh)))
     #+(or)
     (fl.plot:plot (refine mesh)))

  )

;;; (test-triangulate)
(fl.tests:adjoin-test 'test-triangulate)
