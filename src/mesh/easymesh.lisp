;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; easymesh.lisp - Interface to easymesh
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

(in-package :cl-user)
(defpackage "FL.EASYMESH"
  (:use "COMMON-LISP" "FL.MACROS" "FL.UTILITIES" "FL.DEBUG"
	"FL.MATLISP" "FL.ALGEBRA" "FL.FUNCTION"
	"FL.MESH")
  (:export )
  (:documentation "This package provides an interface to the 2D
triangulation software @program{easymesh}."))

(in-package "FL.EASYMESH")

(defvar *easymesh-pathname*
  (or (aand cl-user::*easymesh-path* (probe-file (pathname it)))
      (fl.port:find-executable "easymesh")
      (probe-file #p"femlisp:external;easymesh"))
  "Pathname of the @program{easymesh} binary.")

(defvar *meshes-pathname*
  (or (aand cl-user::*meshes-directory* (pathname it))
      (translate-logical-pathname #p"femlisp:meshes;"))
  "Pathname of the directory for @femlisp{} meshes.")

(defvar *meshes-directory*
  (namestring *meshes-pathname*)
  "Namestring of @var{*meshes-pathname*}.")

(defun write-easymesh-boundary (domain mesh-size)
  "Writes out the boundary for @arg{domain}.  Returns a hash-table mapping
objects to ids."
  (let ((pathname (make-pathname :name "mesh.d"
				 :directory (pathname-directory *meshes-pathname*)))
	(boundary (domain-boundary domain))
	(object->id (make-hash-table)))
    ;; wait until the file is moved by the graphics program
    (loop while (probe-file pathname) do (sleep 0.1))
    ;; write output
    (with-open-file (stream pathname :direction :output :if-exists :supersede)
      (let ((global-count -1))
	(let ((count -1))
	  (format stream "~D~%" (nr-of-cells boundary 0))
	  (doskel (vtx boundary :dimension 0)
	    (let ((pos (vertex-position vtx)))
	      (format stream "~D: ~F ~F ~F ~D~%"
		      (incf count) (aref pos 0) (aref pos 1)
		      (if (functionp mesh-size) (funcall mesh-size pos) mesh-size)
		      (incf global-count))
	      (setf (gethash vtx vertices) t)
	      (setf (gethash vtx object->id) global-count))))
	(let ((components (collect-components boundary)))
	  ...
	(let ((count -1))
	  (flet ((collect-component (vertex)
		   (if (remhash vertex vertices)
		       ()
		       (let ((next-edge (pop (gethash vertex vertex->edge))))
			 (cons next-edge
			       (collect-component
				(let ((corners (boundary next-edge)))
				  (aref corners (if (eq vertex (aref corners 0)) 1 0)))))))))
	    
			       (aref corners 1)
			       (aref corners 0)
			       (get-arbitrary-key-from-hash-table vertices)
			       (loop for vertex = 
				     (components ()))
			       (loop for vertex = (pop vertices))))))))

(defun collect-components (boundary)
  (let ((vertices (copy-hash-table (etable boundary 0)))
	(edges (copy-hash-table (etable boundary 1)))
	(vertex->edge (make-hash-table)))
    (doskel (edge boundary :dimension 1)
      (let ((corners (boundary edge)))
	(push edge (gethash (aref corners 0) vertex->edge))
	(push edge (gethash (aref corners 1) vertex->edge))))
    (loop
     while (plusp (hash-table-count vertices)) collecting
     (loop for vertex = (get-arbitrary-key-from-hash-table vertices) then
	   (aref (boundary edge)
		 (if (nth-value 1 (gethash (aref (boundary edge) 0) vertices))
		     0 1))
	   for (edge-1 edge-2) = (gethash vertex vertex->edge)
	   for edge = (if (nth-value 1 (gethash edge-1 edges)) edge-1 edge-2)
	   do (format t "~%vertex=~A, edge=~A~%~A~%~A~%"
		      vertex edge (hash-table-keys vertices)
		      (hash-table-keys edges))
	   while (and (remhash vertex vertices) (remhash edge edges))
	   collect vertex collect edge))))
	
				     (format stream "~D~%" nr-edges)
				     (doskel (edge boundary :dimension 1)
				       (let ((sides (boundary edge)))
					 (format stream "~D: ~D ~D ~D~%" (incf count)
						 (gethash (aref sides 1) object->id)
						 (gethash (aref sides 0) object->id)
						 (incf global-count))
					 (setf (gethash edge object->id) global-count))))))
    object->id))

(write-easymesh-boundary (n-cube-domain 2) (constantly 0.25))
	
      (apply #'graphic-write-data stream object program rest))
  ;; wait until the file is there
  (loop until (probe-file pathname) do (sleep 0.1))
  ;; send script commands to plot program
  (apply #'send-graphic-commands object program rest)

  )

(defun read-easymesh-mesh ()
  )

(defmethod triangulate (domain &key (program :easymesh) &allow-other-keys)
  
  )