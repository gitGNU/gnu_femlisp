;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; coarsen.lisp
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

(defun coarsen (h-mesh &key test (table (make-hash-table)))
  "Coarsen all cells in h-mesh whose children are marked by the predicate
TEST or are contained in the hash-table TABLE."
  (when test
    (doskel (cell h-mesh)
      (when (funcall test cell)
	(setf (gethash cell table) t))))
  ;; mark all subcells
  (doskel (cell h-mesh :direction :down)
    (when (gethash cell table)
      (loop for side across (boundary cell) do
	    (setf (gethash side table) t))))
  ;; now coarsen while keeping the mesh consistent
  (loop for level from (toplevel h-mesh) downto 0
	for skel = (cells-on-level h-mesh level) do
	(loop for dim from (dimension skel) down to 0
	      (doskel (cell skel :dimension dim)
		(when (gethash cell table)
		  (whereas ((children (children cell h-mesh)))
		    (when (every (rcurry gethash table) children)
		      (setf (children cell h-mesh) nil)
		      ;; remove all subcells from meshes
		      (loop for child across children do
			    (remove-cell skel child)
			    (remove-cell h-mesh child))))
		  (remhash cell table)))
	      ;; leave all refined boundaries of remaining refined cells
	      (doskel (cell skel :dimension dim)
		(loop for side across (boundary cell) do
		      (remhash side table)))))
  (assert (zerop (hash-table-count table)) ()
	  "Error: Unexpectedly there are cells left for coarsening."))
