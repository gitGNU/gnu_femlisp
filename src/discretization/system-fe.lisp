;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; system-fe.lisp - generally useful stuff for system discretizations
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

(in-package :fl.discretization)

(defmethod essential-boundary-constraints
    ((problem <pde-problem>) (ansatz-space <ansatz-space>)
     &key level (where :surface) interface)
  (assert (or (not (eq where :surface)) interface))
  (let* ((problem (problem ansatz-space))
	 (h-mesh (hierarchical-mesh ansatz-space))
	 (level-skel (if level (cells-on-level h-mesh level) h-mesh))
	 (fe-class (fe-class ansatz-space))
	 (constraints-P (make-ansatz-space-automorphism ansatz-space))
	 (constraints-Q (make-ansatz-space-automorphism ansatz-space))
	 (constraints-rhs (make-ansatz-space-vector ansatz-space)))
    (doskel (cell level-skel)
      (when (ecase where
	      (:refined (refined-p cell h-mesh))
	      (:surface (or (not (refined-p cell h-mesh))
			    (member-of-skeleton? cell interface)))
	      (:all t))
	(let* ((coeffs (coefficients-of-cell cell h-mesh problem))
	       (constraint-id (constraint-identifier problem cell))
	       (constraint-function (getf coeffs constraint-id))
	       (cell-key (cell-key cell h-mesh)))
	  ;; old version
	  #+(or)
	  (when constraint-function
	    (loop with fe = (get-fe fe-class cell)
		  for dof in (fe-dofs fe)
		  when (interior-dof? dof) do
		  (let ((k (dof-in-vblock-index dof))
			(comp (dof-component dof))
			(ci (list :local (dof-coord dof)
				  :global (local->global cell (dof-gcoord dof)))))
		    ;; Very simple component-wise constraints for Lagrange
		    ;; fe. This should be extended to constraints in the
		    ;; form P, Q, r.  It could also be accelerated by not
		    ;; evaluating every bc several times.
		    (multiple-value-bind (constraint-flag constraint-val)
			(evaluate constraint-function ci)
		      (when (aref constraint-flag comp)
			(setf (mref (mref constraints-P cell-key cell-key) k k) 1.0)
			(minject (aref constraint-val comp) (vref constraints-rhs cell-key)
				 k 0))))))
	  #-(or) ; new version
	  (when constraint-function
	    (let ((fe (get-fe fe-class cell)))
	      (do-dof (dof fe)
		(when (interior-dof? dof)
		  (let ((k (dof-in-vblock-index dof))
			(comp (dof-component dof))
			(ci (list :local (dof-coord dof)
				  :global (local->global cell (dof-gcoord dof)))))
		    ;; Simple component-wise constraints for Lagrange fe.
		    ;; It could be accelerated by not evaluating every bc
		    ;; several times in the case of equal dofs for several
		    ;; components.
		    (multiple-value-bind (constraint-P constraint-val)
			(evaluate constraint-function ci)
		      (whereas ((cp-val (aref constraint-P comp)))
			(let ((cp-val   ; compatibility: T means componentwise
			       (if (eq cp-val T) (eye 1) cp-val)))
			  ;; if values are given blockwise, the following
			  ;; entries must be given no constraint
			  (assert (loop for i from 1 below (nrows cp-val)
				     never (or (aref constraint-val (+ comp i))
					       (aref constraint-P (+ comp i)))))
			  ;; this is not yet correct for blockwise
			  ;; assignments and order p>2, because on edges a
			  ;; nonconsecutive submatrix should be set
			  (minject cp-val (mref constraints-P cell-key cell-key) k k)
			  (minject (aref constraint-val comp)
				   (vref constraints-rhs cell-key) k 0)
			  )))))))))))
    (values constraints-P constraints-Q constraints-rhs)))

(defun dof-submatrix-indices (dof vecfe n)
  "Find the indices of dofs corresponding to @arg{n} consecutive components
with the same dof location."
  )
