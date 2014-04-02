;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; tetgen.lisp - Femlisp interface to Hang Si's Tetgen
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;
;;; Copyright (C) 2010 Nicolas Neuss, KIT Karlsruhe
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
;;; MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
;;; IN NO EVENT SHALL THE AUTHOR, THE KARLSRUHE INSTITUTE OF TECHNOLOGY OR
;;; OTHER CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
;;; SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
;;; LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
;;; DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
;;; THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
;;; (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
;;; OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(in-package :fl.mesh)

(file-documentation
 "This file provides an interface to the 3D triangulation software
@program{tetgen} written by Hang Si (WIAS Berlin).")

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Tetgen basics
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defvar *tetgen-pathname*
  (or (aand fl.start::*tetgen-path* (probe-file (pathname it)))
      (fl.port:find-executable "tetgen")
      (fl.start:femlisp-pathname "external/tetgen/tetgen"))
  "Pathname of the @program{tetgen} binary.")

(defun call-tetgen (&key meshsize &allow-other-keys)
  "Calls Hang Si's Tetgen."
  (unless *tetgen-pathname*
    (error "Cannot call Tetgen.  Please install it, e.g. by issuing the
command 'make tetgen' from the Femlisp main directory."))
  (dbg :tetgen "Deleting previous Tetgen output...")
  (loop for type in '(:element :face :edge :node)
     do (awhen (probe-file (mesh-file type))
          (delete-file it)))
  (dbg :tetgen "Calling Tetgen...")
  (fl.port:run-program
   *tetgen-pathname*
   (remove nil
           (list (concatenate 'string "-YqeB" (dbg-when :tetgen "V"))
                 (when meshsize (format nil "-a~A"
                                        (* (/ (sqrt 2) 12) (expt meshsize 3))))
                 (concatenate 'string (namestring (meshes-pathname))
                              (pathname-name (mesh-file :poly)))))
   :wait t :output (dbg-when :tetgen *trace-output*))
  (assert (probe-file (mesh-file :element)) ()
          "Tetgen did not succeed"))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Read Tetgen meshes
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun read-nodes (stream nr-nodes &optional (dim 3))
  (lret ((result (make-array (1+ nr-nodes) :initial-element nil)))
    (loop repeat nr-nodes do
       (destructuring-bind (node-id xc yc zc)
           (read-next-numbers-line stream)
         (setf (aref result node-id)
               (make-vertex (ecase dim
                              (2 (double-vec xc yc))
                              (3 (double-vec xc yc zc))
                              )))))))

(defun insert-simplex-in-table (node-ids table)
  (ensure (gethash node-ids table)
          (make-simplex
           (map 'vector
                (_ (insert-simplex-in-table (remove _ node-ids) table))
                node-ids))))

(defun read-tetgen-simplices (stream vertices &optional (n 0))
  (declare (optimize debug))
  (lret ((ids->cell (make-hash-table :test 'equalp)))
    ;; populate with vertices
    (loop+ ((vertex vertices) i)
       when vertex do
       (setf (gethash (list i) ids->cell) vertex))
    (loop for i from 0
       and line = (read-next-numbers-line stream)
       while line
       do
       (let ((node-ids (safe-sort (nthcdr n line) #'<)))
         (dbg :triangulate "Insert simplex: ~A~%" node-ids)
         (insert-simplex-in-table node-ids ids->cell)))))

(defun read-tetgen-mesh ()
  "Reads a mesh from Tetgen output."
  (declare (optimize debug))
  (let ((nodes
         (with-open-file (stream (mesh-file :node))
           (destructuring-bind (nr-nodes dim nr-attributes nr-markers)
               (read-next-numbers-line stream)
             (dbg :tetgen "nr-nodes=~D dim=~D nr-attributes=~D nr-markers=~D"
                  nr-nodes dim nr-attributes nr-markers)
             (assert (and (= 3 dim) (zerop nr-attributes) (zerop nr-markers)))
             (read-nodes stream nr-nodes)))))
    (with-open-file (stream (mesh-file :element))
      (destructuring-bind (nr-elements nodes-per-cell attributes)
	  (read-next-numbers-line stream)
        (dbg :tetgen "nr-elements=~D nodes-per-cell=~D nr-attributes=~D"
                  nr-elements nodes-per-cell attributes)
	(assert (and (= 4 nodes-per-cell) (zerop attributes)))
        (let ((simplices
               (hash-table-values (read-tetgen-simplices stream nodes 1))))
          (make-instance '<skeleton> :dimension 3 :cells simplices))))))

(defun read-tetgen-mesh-for (interior-patch &optional mesh)
  "Reads a mesh from Tetgen output."
  (declare (optimize debug))
  (let ((dim (dimension interior-patch)))
    (ensure mesh (make-instance '<skeleton>f :dimension dim))
    (let ((nodes
           (with-open-file (stream (mesh-file :node))
             (destructuring-bind (nr-nodes file-dim nr-attributes nr-markers)
                 (read-next-numbers-line stream)
               (dbg :tetgen "nr-nodes=~D file-dim=~D nr-attributes=~D nr-markers=~D"
                    nr-nodes file-dim nr-attributes nr-markers)
               (assert (and (= file-dim dim) (zerop nr-attributes) (zerop nr-markers)))
               (read-nodes stream nr-nodes dim))))
          (corners->cell (make-hash-table :test 'equalp)))
      ;; insert known part (the boundary)
      (doskel (cell mesh)
        (setf (gethash (corners cell) corners->cell) cell))
      (with-open-file (stream (mesh-file :element))
        (destructuring-bind (nr-elements nodes-per-cell attributes)
            (read-next-numbers-line stream)
          (dbg :tetgen "nr-elements=~D nodes-per-cell=~D nr-attributes=~D"
               nr-elements nodes-per-cell attributes)
          (assert (and (= (1+ dim) nodes-per-cell) (zerop attributes)))
          (loop
             repeat nr-elements
             for (cell-id . node-ids) = (read-next-numbers-line stream)
             for sorted-node-ids = (sort node-ids #'<)
             for corners = (mapcar (_ (vertex-position (aref nodes _))) sorted-node-ids)
             unless (gethash corners corners->cell)
             do
             (dbg :triangle "New cell ~D: ~A" cell-id sorted-node-ids)
             (insert-cell-from-corners
              mesh corners->cell (simplex-class dim) corners
              (list 'PATCH interior-patch) :create-subcells t)))))))

;;; (fl.plot:plot (read-tetgen-mesh) :program :vtk)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Write 2D-skeleton in poly-format
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun write-2d-skeleton-to-poly-file (patch mesh vertices)
  "Writes out the 2d-skeleton for @arg{mesh} to a @code{.poly}-file.
Returns a hash-table mapping the cells mesh to ids, a vector mapping ids to
vertices and a vector mapping ids to lines.  @arg{vertices} is a list of
vertices defining the correct ordering of vertices."
  (let ((cell->id (make-hash-table)))
    ;; if the file is already there, we move it
    (when (probe-file (mesh-file :poly))
      (awhen (probe-file (mesh-file :poly-backup))
        (delete-file it))
      (rename-file (mesh-file :poly) (mesh-file :poly-backup)))
    ;; write output
    (with-open-file (stream (mesh-file :poly) :direction :output
			    :if-exists :supersede)
      ;; nodes
      (format stream "~D 3 0 0~%" (length vertices))
      (loop for count from 1 and vtx in vertices do
           (let ((pos (vertex-position vtx)))
             (format stream "~D~{ ~F~}~%" count
                     (coerce pos 'list)))
           (setf (gethash vtx cell->id) count))
      ;; facets
      (format stream "~D 1~%" (nr-of-cells mesh 2))
      (let ((count 0))
	(doskel (cell mesh :dimension 2)
          (unless (simplex-p cell)
            (error "Cannot handle non-triangular boundary meshes."))
          (format stream "1 0 ~D~%3~{ ~D~}~%"
                  (incf count)
                  (mapcar (rcurry #'gethash cell->id) (vertices cell)))))
      ;; holes
      (let ((holes (and (typep patch '<boundary-cell>)
			(holes patch))))
	(format stream "~D~%" (length holes))
	(loop for count from 1 and hole in holes do
	      (format stream "~D~{ ~F~}~%"
                      count (coerce hole 'list))))
      )
    ;; return the translation table
    cell->id))

#+(or)
(let* ((mesh (skeleton (n-simplex 3)))
       (patch (first (cells-of-highest-dim mesh)))
       (vertices (vertices patch)))
  (write-2d-skeleton-to-poly-file patch mesh vertices))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Tetrahedration of 3D-domains
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmethod extend-triangulation (mesh (dim (eql 3)) (embedded-dim (eql 3))
                                 (program (eql :tetgen))
                                 &key meshsize &allow-other-keys)
  "Calls @program{tetgen} for each 3D patch of the domain.  @arg{mesh}
should already be generated for 0,1 and 2-dimensional patches."
  (let ((domain (domain mesh)))
    (doskel (patch domain :dimension 3)
      (let* ((patch-boundary-mesh
	      (skeleton
	       (let ((patches (subcells patch)))
		 (find-cells (_ (find (patch-of-cell _ mesh) patches))
			     mesh)))))
	(write-2d-skeleton-to-poly-file patch patch-boundary-mesh
                                        (consistent-vertex-numbering patch-boundary-mesh))
      (call-tetgen :meshsize meshsize)
      (read-tetgen-mesh-for patch mesh)
      ))
    
    ;; introduce blending maps if necessary
    (doskel (cell mesh :dimension '(2 . 3))
      (dovec (side (boundary cell))
        (when (and (mapped-p side)
                   (not (mapped-p cell)))
          (change-to-blending cell))))
    ))

;;; (write-boundary-to-poly-file (n-cube-domain 2) (constantly 10))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Test Tetgen on a sample domain with a curved boundary
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun sample-domain ()
  "Generates a sample domain which is a modified unit-cube whose upper side
bends upwards in the x_0-direction."
  (lret ((domain
          (make-domain (copy-skeleton (skeleton (n-cube 3))))))
    (flet ((osc (x) (1+ (* .2 (sin (* pi x)))))
           (dosc (x) (* .2 pi (cos (* pi x)))))
      (doskel (cell domain)
        (let ((midpoint (midpoint cell)))
          (when (= (aref midpoint 0) 0.5)
            (cond
              ((= (dimension cell) 3)
               (change-class cell 'fl.mesh::<boundary-cell>
                             :dimension 3 :midpoint #d(0.5 0.5 0.5)))
              ((= (aref midpoint 2) 1.0)
               (case (dimension cell)
                 (1 (change-class
                     cell 'fl.mesh::<mapped-1-simplex>
                     :mapping
                     (make-instance
                      '<special-function>
                      :domain-dimension 1 :image-dimension 3
                      :evaluator
                      (_ (double-vec (aref _ 0) (aref midpoint 1) (osc (aref _ 0))))
                      :gradient
                      (_ (make-real-matrix
                          `((1.0)
                            (0.0)
                            (,(dosc (aref _ 0)))))))))
                 (2 (change-class
                     cell 'fl.mesh::<mapped-1-1-product-cell>
                     :mapping
                     (make-instance
                      '<special-function>
                      :domain-dimension 2 :image-dimension 3
                      :evaluator
                      (_ (double-vec (aref _ 0) (aref _ 1) (osc (aref _ 0))))
                      :gradient
                      (_ (make-real-matrix
                          `((1.0 0.0)
                            (0.0 1.0)
                            (,(dosc (aref _ 0)) 0.0)))))))))
              ((and (= (dimension cell) 2)
                    (plusp (aref midpoint 2)))
               (change-class
                cell 'fl.mesh::<mapped-1-1-product-cell>
                :mapping
                (make-instance
                 '<special-function>
                 :domain-dimension 2 :image-dimension 3
                 :evaluator
                 (_ (double-vec (aref _ 0) (aref midpoint 1)
                                (* (aref _ 1) (osc (aref _ 0)))))
                 :gradient
                 (_ (make-real-matrix
                     `((1.0 0.0)
                       (0.0 0.0)
                       (,(* (aref _ 1) (dosc (aref _ 0))) ,(osc (aref _ 0))))))))))))))
    (check domain)))

(defun sample-mesh (n &optional (check-p t))
  (lret ((mesh (triangulate (sample-domain) :meshsize 0.5 :check-p check-p)))
    (when check-p (check mesh))
    (loop repeat n do
         (setf mesh (refine mesh))
         (when check-p (check mesh)))
    (change-class mesh '<hierarchical-mesh>)))

;;;; Testing:
(defun test-tetgen ()
  (dbg-on :triangulate)
  (dbg-on :tetgen)
  (dbg-off)
  (time (check (sample-mesh 0)))
  (time (check (sample-mesh 1)))
  #+(or)
  (fl.application:storing
    (fl.plot:plot (sample-mesh 0 nil)))

  )

;;; (test-tetgen)

(fl.tests:adjoin-test 'test-triangle)