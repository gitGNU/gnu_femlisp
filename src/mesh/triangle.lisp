
;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; triangle.lisp - Interface to Shewchuk's Triangle
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

(file-documentation
 "This file provides an interface to the 2D triangulation software
@program{triangle} written by Jonathan R. Shewchuk.")

(defvar *triangle-pathname*
  (or (aand fl.start::*triangle-path* (probe-file (pathname it)))
      (fl.port:find-executable "triangle")
      (probe-file #p"femlisp:external;triangle;triangle"))
  "Pathname of the @program{triangle} binary.")

(defvar *meshes-pathname*
  (or (aand fl.start::*meshes-directory* (pathname it))
      (translate-logical-pathname #p"femlisp:meshes;"))
  "Pathname of the directory for @femlisp{} meshes.")

(defun mesh-file (kind)
  "Returns the pathname for the mesh component @arg{kind} which may be
:poly, :edge, :node, :element."
  (make-pathname
   :name (concatenate
	  'string "mesh."
	  (ecase kind
	    (:poly "poly") (:poly-backup "poly.bak")
	    (:element "1.ele") (:edge "1.edge") (:node "1.node")))
   :directory (pathname-directory *meshes-pathname*)))

(defun write-1d-skeleton-to-poly-file (patch mesh vertices)
  "Writes out the 1d-skeleton for @arg{mesh} to a @code{.poly}-file.
Returns a hash-table mapping the cells mesh to ids, a vector mapping ids to
vertices and a vector mapping ids to lines.  @arg{vertices} is a list of
vertices defining the correct ordering of vertices."
  (let ((cell->id (make-hash-table)))
    ;; if the file is already there, we move it
    (when (probe-file (mesh-file :poly))
      (rename-file (mesh-file :poly) (mesh-file :poly-backup)))
    ;; write output
    (with-open-file (stream (mesh-file :poly) :direction :output
			    :if-exists :supersede)
      ;; nodes
      (format stream "~D 2 0 0~%" (length vertices))
      (loop for count from 1 and vtx in vertices do
	    (let ((pos (vertex-position vtx)))
	      (format stream "~D ~F ~F~%"
		      count (aref pos 0) (aref pos 1))
	      (setf (gethash vtx cell->id) count)))
      ;; edges
      (format stream "~D 0~%" (nr-of-cells mesh 1))
      (let ((count 0))
	(doskel (line mesh :dimension 1)
	  (let ((sides (boundary line)))
	    (format stream "~D ~D ~D~%"
		    (incf count)
		    (gethash (aref sides 1) cell->id)
		    (gethash (aref sides 0) cell->id)))))
      ;; holes
      (let ((holes (and (typep patch '<boundary-cell>)
			(holes patch))))
	(format stream "~D~%" (length holes))
	(loop for count from 1 and hole in holes do
	      (format stream "~D ~F ~F~%" count (aref hole 0) (aref hole 1))))
      )
    ;; return the translation table
    cell->id))

;;; (write-boundary-to-poly-file (n-cube-domain 2) (constantly 10))

(defun call-triangle ()
  "Calls Shewchuk's triangle program."
  (unless *triangle-pathname*
    (error "Cannot call TRIANGLE.  Please install it, e.g. by issuing the
command 'make triangle' from the Femlisp main directory."))
  (fl.port:run-program
   *triangle-pathname*
   (list (concatenate 'string "-YqeB" (dbg-when :triangle "V"))
	 (concatenate 'string (namestring *meshes-pathname*)
		      (pathname-name (mesh-file :poly))))
   :wait t :output (dbg-when :triangle *trace-output*)))

(defun read-next-uncommented-line (stream)
  "Reads the next line from @arg{stream} which does not begin with
@code{#}."
  (loop for string = (read-line stream)
	unless string do (error "Unexpected EOF.")
	until (char-not-equal #\# (elt string 0))
	finally
	(return (read-from-string (concatenate 'string "(" string ")")))))

(defun read-triangle-mesh-for (mesh interior-patch)
  "Extends @arg{mesh} by the output from Triangle for meshing
@arg{interior-patch}."
  (declare (optimize (debug 3)))
  (let (nodes ; will be an array, its size is not yet known
	(corners->cell (make-hash-table :test 'equalp)))
    ;; insert known part (the boundary)
    (doskel (cell mesh)
      (setf (gethash (corners cell) corners->cell) cell))
    ;; read and insert interior vertices
    (with-open-file (stream (mesh-file :node))
      (destructuring-bind (nr-nodes dim nr-attributes nr-markers)
	  (read-next-uncommented-line stream)
	(assert (and (= 2 dim) (zerop nr-attributes) (zerop nr-markers)))
	(setq nodes (make-array (1+ nr-nodes)))
	(loop
	 repeat nr-nodes
	 for (node-id xc yc attribute) = (read-next-uncommented-line stream)
	 for pos = (map 'double-vec (rcurry #'coerce 'double-float)
			(vector xc yc))
	 for corners = (list pos) do
	 (setf (aref nodes node-id) pos)
	 (let ((found (gethash corners corners->cell)))
	   (dbg :triangle "~:[New~;Old~] node ~D: ~A" found node-id pos)
	   (unless found
	     (let ((vtx (make-vertex pos)))
	       (setf (skel-ref mesh vtx) (list 'PATCH interior-patch))
	       (setf (gethash corners corners->cell) vtx)))))))
    ;; read and insert interior edges (should be done automatically)
    (dbg-when :triangle
      (doskel (edge mesh :dimension 1)
	(dbg :triangle "Old edge: ~A"
	     (loop for i from 1 downto 0 collect
		   (vertex-position (aref (boundary edge) i))))))
    (with-open-file (stream (mesh-file :edge))
      (destructuring-bind (nr-edges nr-markers)
	  (read-next-uncommented-line stream)
	(assert (zerop nr-markers))
	(loop
	 repeat nr-edges
	 for (edge-id node-1 node-2) = (read-next-uncommented-line stream)
	 for node-ids = (sort (list node-1 node-2) #'<)
	 for corners = (mapcar (curry #'aref nodes) node-ids) do
	 (let ((found (gethash corners corners->cell)))
	   (dbg :triangle "~:[New~;Old~] edge ~D: ~A" found edge-id node-ids)
	   (unless found
	     (insert-cell-from-corners mesh corners->cell (simplex-class 1) corners
				       (list 'PATCH interior-patch)
				       :create-subcells nil))))))
    ;; read and insert triangles
    (with-open-file (stream (mesh-file :element))
      (destructuring-bind (nr-triangles nodes-per-triangle attributes)
	  (read-next-uncommented-line stream)
	(assert (and (= 3 nodes-per-triangle) (zerop attributes)))
	(loop
	 repeat nr-triangles
	 for (triangle-id node-1 node-2 node-3) = (read-next-uncommented-line stream)
	 for node-ids = (sort (list node-1 node-2 node-3) #'<)
	 for corners = (mapcar (curry #'aref nodes) node-ids)
	 unless (gethash corners corners->cell) do
	 (dbg :triangle "New cell ~D: ~A" triangle-id node-ids)
	 (insert-cell-from-corners
	  mesh corners->cell (simplex-class 2) corners
	  (list 'PATCH interior-patch) :create-subcells nil))))
    mesh))

(defmethod triangulate-2d (mesh (program (eql :triangle))
			   &key vertex-table)
  "Calls @program{triangle} for each 2D patch of the domain.  @arg{mesh}
should already be generated for 0 and 1-dimensional patches.  The mesh
should not cover any 2-dimensinal patches."
  (assert (zerop (nr-of-cells mesh 2))) ; might be relaxed later on
  (let ((domain (domain mesh)))
    (doskel (patch domain :dimension 2)
      (let* ((patch-boundary-mesh
	      (skeleton
	       (let ((patches (subcells patch)))
		 (find-cells #'(lambda (cell) (find (patch-of-cell cell mesh) patches))
			     mesh))))
	     (vertices (sort (cells-of-dim patch-boundary-mesh 0)
			     #'(lambda (v1 v2)
				 (< (gethash v1 vertex-table) (gethash v2 vertex-table))))))
	(write-1d-skeleton-to-poly-file patch patch-boundary-mesh vertices))
      (call-triangle)
      (read-triangle-mesh-for mesh patch))))

;;;; Testing: (test-triangle)
(defun test-triangle ()
  (dbg-off)
  (dbg-on :triangulate)
  (check (triangulate (n-cube-domain 2)))
  (nr-of-cells (triangulate (n-cube-domain 2) :meshsize 0.2) 1)
  (check (triangulate (n-ball-domain 2) :meshsize 0.5))
  (dbg-off :triangle)
  (when *triangle-pathname*
    (triangulate (n-ball-domain 2) :patch->resolution
		 #'(lambda (patch)
		     (if (< (norm (midpoint patch)) 0.8) 1 20))))
  )

(fl.tests:adjoin-test 'test-triangle)