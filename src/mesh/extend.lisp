;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;  extend.lisp
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

(defun extend (mesh &key (test (constantly t)))
  "Extends a mesh on all extensible cells for which test ---if provided---
yields T."
  (let ((domain (domain mesh)))
    (when (extensible-p domain)
      (let* ((level-0 (if (hierarchical-mesh-p mesh)
			  (cells-on-level mesh 0)
			  mesh))
	     (extensible-cells
	      (find-cells
	       #'(lambda (cell)
		   (and (get-patch-property cell mesh 'EXTENSION)
			(funcall test cell)))
	       level-0))
	     (changed-cells
	      (loop for cell in extensible-cells nconcing
		    (multiple-value-call #'nconc
		      (funcall (get-patch-property cell mesh 'EXTENSION)
			       cell mesh)))))
	(when (hierarchical-mesh-p mesh)
	  ;; insert an entry also in level-0 skeleton
	  (dolist (cell changed-cells)
	    (setf (skel-ref level-0 cell)
		  (skel-ref mesh cell))
	    ;; ensure that patch still agrees
	    (awhen (children cell mesh)
	      (loop for child across it do
		    (setf (patch-of-cell child mesh)
			  (patch-of-cell cell mesh))))
	    ))))))


(defun standard-extender (original-cell replacement)
  "Extension function replacing an original-cell with a replacement."
  #'(lambda (cell mesh)
      (let* ((shift (vec- (midpoint cell) (midpoint original-cell)))
	     (shifted-replacement
	      (shift-skeleton replacement shift
			      :properties '(PATCH EXTENSION))))
	(multiple-value-bind (skel-1 skel-2 overlap)
	    (skel-add! mesh shifted-replacement
		       :override '(PATCH EXTENSION)
		       :active-skel-1
		       (if (hierarchical-mesh-p mesh) (cells-on-level mesh 0) mesh))
	  (declare (ignore skel-1))
	  (values
	   (find-cells #'(lambda (cell) (not (gethash cell overlap)))
		       skel-2)
	   (hash-table-values overlap))))))

(defun cube-extender (domain-cube direction)
  "Makes domain-cube ---which should be a cube in a domain--- extensible in
the given direction."
  (let* ((replacement (skeleton domain-cube))
	 (cube (car (cells-of-highest-dim replacement)))
	 (old-ext (aref (boundary cube) (1+ (* 2 direction))))
	 (new-ext (aref (boundary cube) (* 2 direction)))
	 (extension (standard-extender old-ext replacement)))
    ;; setup patches of replacement
    (loop for subcell-1 across (subcells cube)
	  and subcell-2 across (subcells domain-cube) do
	  (setf (get-cell-property subcell-1 replacement 'PATCH) subcell-2))
    ;; correct patch for old boundary to be in the interior
    (setf (get-cell-property old-ext replacement 'PATCH)
	  domain-cube)
    ;; setup extensions
    (setf (get-cell-property old-ext replacement 'EXTENSION) nil)
    (setf (get-cell-property new-ext replacement 'EXTENSION) extension)
    ;; and return extension
    (values extension replacement)))


;;;; Testing

(defun test-extend ()
  (let* ((dim 1) (direction 0)
	 (direction (min direction (1- dim)))
	 (domain (n-cube-domain dim))
	 (domain-cube (car (cells-of-highest-dim domain)))
	 (domain-ext (aref (boundary domain-cube) (* 2 direction)))
	 (extension (cube-extender domain-cube direction)))
    ;; make domain extensible
    (setf (get-cell-property domain-ext domain 'EXTENSION) extension)
    (ensure-secondary-information domain)
    ;; make mesh
    (let ((mesh (make-hierarchical-mesh-from-domain domain)))
      (refine mesh)
      (extend mesh)
      #+(or)(plot:plot mesh)))
  )

(tests:adjoin-femlisp-test 'test-extend)
  
