;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; triangle.lisp - Interface to Shewchuk's Triangle
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;
;;; Copyright (C) 2004 Nicolas Neuss, University of Heidelberg.
;;; Copyright (C) 2010 Nicolas Neuss, KIT Karlsruhe.
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
;;; IN NO EVENT SHALL THE AUTHOR, THE UNIVERSITY OF HEIDELBERG, THE KIT,
;;; OR OTHER CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
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
 "This file provides an interface to the 2D triangulation software
@program{Triangle} written by Jonathan R. Shewchuk.")

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Triangle basics
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun meshes-pathname ()
  "Returns the pathname of the directory for @femlisp{} meshes."
  (or (aand fl.start::*meshes-directory* (pathname it))
      (fl.port::getenv "FEMLISP_MESHES")
      (fl.start:femlisp-pathname "meshes/")))

(defvar *meshfile-basename*
  "mesh"
  "Basename for mesh decomposition.")

(defvar *triangle-pathname*
  (or fl.start::*triangle-path*
      (fl.port:find-executable "triangle")
      (probe-file (fl.start:femlisp-pathname "external/triangle/triangle")))
  "Pathname of the @program{triangle} binary.")

(defun mesh-file (kind)
  "Returns the pathname for the mesh component @arg{kind} which may be
:poly, :edge, :node, :element and :face."
  (make-pathname
   :name (format nil "~A.~A" *meshfile-basename*
                 (ecase kind
                   (:poly "poly")
                   (:poly-backup "poly.bak")
                   (:element "1.ele")
                   (:face "1.face")
                   (:edge "1.edge")
                   (:node "1.node")))
   :directory (pathname-directory (meshes-pathname))))

(defun call-triangle (&key meshsize &allow-other-keys)
  "Calls Shewchuk's triangle program."
  (unless *triangle-pathname*
    (error "Cannot call Triangle.  Please install it, e.g. by issuing the
command 'make triangle' from the Femlisp main directory."))
  (dbg :triangle "Deleting previous Triangle output...")
  (loop for type in '(:element :face :edge :node)
     do (awhen (probe-file (mesh-file type))
          (delete-file it)))
  (dbg :triangle "Calling Triangle...")
  (fl.port:run-program
   *triangle-pathname*
   (remove nil
           (list (concatenate 'string "-YqeB" (dbg-when :triangle "V"))
                 (when meshsize (format nil "-a~A"
                                        (* (/ (sqrt 3) 8) (square meshsize))))
                 (concatenate 'string (namestring (meshes-pathname))
                              (pathname-name (mesh-file :poly)))))
   :wait t :output (dbg-when :triangle *trace-output*))
  (assert (probe-file (mesh-file :element)) ()
          "Triangle did not succeed"))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Read Triangle meshes
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun read-numbers-from-string (string)
  #-cl-ppcre (read-from-string (concatenate 'string "(" string ")"))
  #+cl-ppcre  ; safer
  (mapcar #'read-from-string
          (cl-ppcre:all-matches-as-strings
           "[+-]?(\\d+\\.\\d+|\\d+\\.|\\.\\d+|\\d+)([eE][+-]?\\d+)?"
           string)))

(defun read-next-uncommented-line (stream)
  "Reads the next nonempty line from @arg{stream} which does not begin with
@code{#}.  NIL indicates EOF."
  (loop for string = (read-line stream nil)
	until (or (null string)
                  (and (not (zerop (length string)))
                       (char-not-equal #\# (elt string 0))))
	finally (return string)))

(defun read-next-numbers-line (stream)
  "Reads the numbers from the next nonempty and uncommented line from
@arg{stream}."
  (aand (read-next-uncommented-line stream)
        (or (read-numbers-from-string it)
            it)))

(defun read-triangle-mesh-for (interior-patch &optional mesh)
  "Builds or extends @arg{mesh} using the output from Triangle for meshing
@arg{interior-patch}."
  (ensure mesh (make-instance '<skeleton> :dimension 2))
  (let (nodes ; will be an array, its size is not yet known
        (corners->cell (make-hash-table :test 'equalp)))
    ;; insert known part (the boundary)
    (doskel (cell mesh)
      (setf (gethash (corners cell) corners->cell) cell))
    ;; read and insert interior vertices
    (with-open-file (stream (mesh-file :node))
      (destructuring-bind (nr-nodes dim nr-attributes nr-markers)
	  (read-next-numbers-line stream)
	(assert (and (= 2 dim) (zerop nr-attributes) (zerop nr-markers)))
	(setq nodes (make-array (1+ nr-nodes)))
	(loop
	 repeat nr-nodes
	 for (node-id xc yc nil) = (read-next-numbers-line stream)
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
	  (read-next-numbers-line stream)
	(assert (zerop nr-markers))
	(loop
	 repeat nr-edges
	 for (edge-id node-1 node-2) = (read-next-numbers-line stream)
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
	  (read-next-numbers-line stream)
	(assert (and (= 3 nodes-per-triangle) (zerop attributes)))
	(loop
           repeat nr-triangles
           for (triangle-id node-1 node-2 node-3) = (read-next-numbers-line stream)
           for node-ids = (sort (list node-1 node-2 node-3) #'<)
           for corners = (mapcar (curry #'aref nodes) node-ids)
           unless (gethash corners corners->cell) do
             (dbg :triangle "New cell ~D: ~A" triangle-id node-ids)
             (insert-cell-from-corners
              mesh corners->cell (simplex-class 2) corners
              (list 'PATCH interior-patch) :create-subcells nil))))
    mesh))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Write 1D skeleton in poly-format
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun write-1d-skeleton-to-poly-file (patch mesh vertices)
  "Writes out the 1d-skeleton for @arg{mesh} to a @code{.poly}-file.
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

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; consistent numbering of vertices
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass priority-queue ()
  ((tables :initform (make-dll))
   (objects :initform (make-hash-table)))
  (:documentation "tables is a dll of hash-tables, at the start is the
  highest priority, at the end the lowest.  objects is a hash-table from
  objects in the pq pointing to the dll-item where the object is."))

(defun pq-get (pq obj)
  "Lookup @arg{obj} from priority queue."
  (with-slots (tables objects) pq
    (whereas ((dli (gethash obj objects)))
      (let ((ht (dli-object dli)))
        (gethash obj ht)))))

(defun (setf pq-get) (value pq obj)
  "Set the value @arg{obj} from priority queue."
  (with-slots (tables objects) pq
    (whereas ((dli (gethash obj objects)))
      (let ((ht (dli-object dli)))
        (setf (gethash obj ht) value)))))

(defun pq-insert (pq obj &optional value)
  "Always insert at highest priority."
  (with-slots (tables objects) pq
    (assert (not (gethash obj objects)))
    (when (dll-empty-p tables)
      (dll-front-insert (make-hash-table) tables))
    (setf (gethash obj (dll-peek-first tables)) value)
    (setf (gethash obj objects) (dll-first tables))))

(defun drop-empty-limits (pq)
  (with-slots (tables) pq
    (when (aand (dll-peek-first tables) (dic-empty-p it))
      (dll-pop-first tables))
    (when (aand (dll-peek-last tables) (dic-empty-p it))
      (dll-pop-last tables))))

(defun pq-remove (pq obj)
  "Remove @arg{obj} from priority queue."
  (with-slots (tables objects) pq
    (whereas ((dli (gethash obj objects)))
      (let ((ht (dli-object dli)))
        (remhash obj ht))
      (remhash obj objects)
      (drop-empty-limits pq))))

(defun pq-pop (pq)
  (with-slots (tables) pq
    (unless (dll-empty-p tables)
      (maphash (lambda (key value)
                 (pq-remove pq key)
                 (return-from pq-pop (values key value)))
               (dll-peek-first tables)))))

(defun increase-priority (pq obj)
  (with-slots (tables objects) pq
    (whereas ((dli (gethash obj objects)))
      (when (eq dli (dll-first tables))
        (dll-front-insert (make-hash-table) tables))
      (let ((pred (dli-pred dli)))
        (setf (gethash obj (dli-object pred))
              (gethash obj (dli-object dli)))
        (remhash obj (dli-object dli))
        (setf (gethash obj objects) pred))
      (drop-empty-limits pq))))

(defun decrease-priority (pq obj)
  (with-slots (tables objects) pq
    (whereas ((dli (gethash obj objects)))
      (when (eq dli (dll-last tables))
        (dll-rear-insert (make-hash-table) tables))
      (let ((succ (dli-succ dli)))
        (setf (gethash obj (dli-object succ))
              (gethash obj (dli-object dli)))
        (remhash obj (dli-object dli))
        (setf (gethash obj objects) succ))
      (drop-empty-limits pq))))

(defmethod show ((pq priority-queue) &key &allow-other-keys)
  (with-slots (tables) pq
    (let ((count 0))
      (dll-for-each (_ (format t "Priority = ~A:~%" (incf count))
                       (display-ht _))
                    tables))))

#+(or)
(let ((pq (make-instance 'priority-queue)))
  (pq-insert pq 1 1)
  (pq-insert pq 2 4)
  (pq-insert pq 3 9)
  (increase-priority pq 2)
  (increase-priority pq 3)
  (increase-priority pq 3)
  (show pq)
  (pq-pop pq)
  (show pq)
  (pq-pop pq)
  (show pq))

(defun consistent-vertex-numbering (skel)
  "Returns a vertex table with consistently numbered vertices - if
possible."
  (let ((pq (make-instance 'priority-queue))
        (vertices ()))
    ;; initialize pq
    (doskel (vtx skel :dimension 0)
      (pq-insert pq vtx ()))
    ;; insert connections
    (doskel (cell skel :dimension 1)
      (let ((conn (cons :connection (vertices cell))))
        (loop for vtx in (rest conn) and k from 0 do
             (push conn (pq-get pq vtx))
             (when (plusp k)
               (decrease-priority pq vtx)))))
    ;; work
    (loop while
         (multiple-value-bind (vtx conns) (pq-pop pq)
           (when vtx
             (loop for conn in conns do
                  (assert (eq (second conn) vtx))
                  (setf (cdr conn) (cddr conn))  ; drop second
                  (awhen (second conn) (increase-priority pq it)))
             (push vtx vertices))))
    ;; result
    (reverse vertices)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Triangulation of 2D-domains
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmethod extend-triangulation (mesh (dim (eql 2)) (embedded-dim (eql 2))
                                 (program (eql :triangle))
                                 &key meshsize &allow-other-keys)
  "Calls @program{triangle} for each 2D patch of the domain.  @arg{mesh}
should already be generated for 0 and 1-dimensional patches."
  (let ((domain (domain mesh)))
    (doskel (patch domain :dimension 2)
      (let* ((patch-boundary-mesh
	      (skeleton
	       (let ((patches (subcells patch)))
		 (find-cells (_ (find (patch-of-cell _ mesh) patches))
			     mesh))))
	     (vertices (consistent-vertex-numbering patch-boundary-mesh)))
	(write-1d-skeleton-to-poly-file patch patch-boundary-mesh vertices))
      (call-triangle :meshsize meshsize)
      (read-triangle-mesh-for patch mesh))
    
    ;; introduce blending for cells meeting a curved boundary
    (doskel (cell mesh :dimension :highest)
      (dovec (side (boundary cell))
        (when (and (mapped-p side)
                   (not (mapped-p cell)))
          (change-to-blending cell))))
    ))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Triangulation of 2D-patches embedded in 3D 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun map-skeleton (func skel)
  (let ((table (make-hash-table))
        (result (make-analog skel)))
    (doskel (cell skel :direction :up)
      (multiple-value-bind (object properties)
          (funcall func cell)
        (let ((copy
               (cond
                 ((cellp object) object)
                 ((vertex-p cell) (make-vertex object))
                 (t (apply #'make-instance
                           (funcall (if object #'mapped-cell-class #'unmapped-cell-class)
                                    (class-of cell))
                           :boundary (vector-map (rcurry #'gethash table) (boundary cell))
                           (and object (list :mapping object)))))))
          (insert-cell! result copy properties)
          (setf (gethash cell table) copy))))
    (values result table)))

(defun patch-boundary-mesh (patch mesh)
  (let ((patches (subcells patch)))
    (subskeleton mesh (lambda (cell) (find (patch-of-cell cell mesh) patches)))))

(defun refpatch-boundary-mesh (patch boundary-mesh)
  (let ((patches (subcells patch))
        (reference-patches (subcells (reference-cell patch))))
    (map-skeleton
     (lambda (cell)
       (let* ((patch (patch-of-cell cell boundary-mesh))
              (pos (position patch patches))
              (refpatch (aref reference-patches pos)))
         (values
          (cond ((vertex-p refpatch) (vertex-position refpatch))
                ((vertex-p cell)
                 (l2g refpatch (get-cell-property cell boundary-mesh 'LOCAL-COORD))))
          (list 'SOURCE-CELL cell))))
     boundary-mesh)))

(defun triangulate-patch (patch boundary-mesh &key meshsize)
  (write-1d-skeleton-to-poly-file patch boundary-mesh
                                  (consistent-vertex-numbering boundary-mesh))
  (call-triangle :meshsize meshsize)
  (lret ((result (make-instance '<skeleton> :dimension 2)))
    (skel-add! result boundary-mesh) 
    (read-triangle-mesh-for patch result)))

(defun patch-mesh (patch refpatch-mesh)
  (map-skeleton
   (lambda (cell)
     (values (or (and (<= (dimension cell) 1)
                      (get-cell-property cell refpatch-mesh 'SOURCE-CELL))
                 (if (vertex-p cell)
                     (local->global patch (vertex-position cell))
                     (compose (cell-mapping patch) (cell-mapping cell))))
             (list 'PATCH patch)))
   refpatch-mesh))

(defmethod extend-triangulation (mesh (dim (eql 2)) (embedded-dim (eql 3))
                                 (program (eql :triangle))
                                 &key meshsize &allow-other-keys)
  "Triangle can only treat 2D domains, so we triangulate 2D reference
patches in the 2D-embedded-in-3D case."
  (declare (optimize debug))
  (let ((domain (domain mesh)))
    (doskel (patch domain :dimension 2)
      (let* ((boundary-mesh (patch-boundary-mesh patch mesh))
             (refpatch-boundary-mesh (refpatch-boundary-mesh patch boundary-mesh))
             (refpatch-mesh (triangulate-patch
                             (reference-cell patch) refpatch-boundary-mesh
                             :meshsize meshsize))
             (patch-mesh (patch-mesh patch refpatch-mesh)))
        
        (skel-add! mesh patch-mesh)))
    ;; no blending here because it would probably do the wrong thing in 3D
    ))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Tests
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(in-suite mesh-suite)

(test triangle
  (finishes
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
  (let* ((vtx (make-vertex #d(1.0 0.0)))
	 (circle (make-line vtx vtx :mapping (circle-function 1.0 #d(0.0 0.0) (* 2 pi))))
	 (ball (make-instance '<boundary-cell>
			      :dimension 2
			      :boundary (vector circle)
			      :midpoint #d(0.0 0.0)))
	 (domain (change-class (skeleton ball) '<domain>)))
    (triangulate domain :meshsize 0.3))
  ))
