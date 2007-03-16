;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; coeffplot.lisp - plotting of coefficient functions
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

(in-package :fl.plot)

(defmethod graphic-commands ((asv <pde-problem>) (program (eql :dx))
			     &rest rest)
  (apply 'fl.graphic::dx-commands-data rest))

(defmethod plot ((problem <pde-problem>) &rest rest &key mesh (refinements 0)
		 (depth 0) (key #'identity) parametric coefficient
		 (rank 0) shape &allow-other-keys)
  "Plots a coefficient function for the problem on the given mesh.  Does
handle coefficients depending on finite element functions."
  (ensure mesh (uniformly-refined-mesh
		(domain problem) refinements :parametric parametric))
  (apply #'graphic-output problem :dx
	 :dimension (dimension mesh)
	 :cells (plot-cells mesh)
	 :rank rank :shape shape
	 :cell->values
	 (lambda (cell)
	   (whereas ((coeff-func (get-coefficient (coefficients-of-cell cell mesh problem)
						  coefficient)))
	     (let* ((sample-points
		     (refcell-refinement-vertex-positions cell depth))
		    (geometry (fe-cell-geometry cell sample-points))
		    (fe-paras (loop for obj in (required-fe-functions (list coeff-func))
				    for asv = (get-property problem
							    (if (symbolp obj) obj (car obj)))
				    for fe = (get-fe (ansatz-space asv) cell)
				    collect obj
				    collect (cons fe (get-local-from-global-vec cell asv))))
		    (fe (when fe-paras
			  (assert (<= (length fe-paras) 2))
			  (prog1 (car (second fe-paras))
			    (setf (second fe-paras) (cdr (second fe-paras))))))
		    (dummy (unless fe
			     (make-array (length sample-points) :initial-element nil))))
	       (map 'vector
		    (lambda (global Dphi Dphi^-1 shape-vals shape-grads)
		      (let* ((gradients (and shape-grads
					     (map 'vector (rcurry #'m* Dphi^-1) shape-grads)))
			     (coeff-input (construct-coeff-input
					   cell global Dphi shape-vals gradients fe-paras)))
			(funcall key (evaluate coeff-func coeff-input))))
		    (getf geometry :global-coords)
		    (getf geometry :gradients)
		    (getf geometry :gradient-inverses)
		    (or dummy (local-evaluation-matrix fe depth))
		    (or dummy (local-evaluation-matrix fe depth :gradient))))))
	 rest))


;;; Testing:

(defun test-coeffplot ()
  (let* ((dim 1)
	 (domain (n-cell-domain dim))
	 (problem
	  (make-instance
	   '<pde-problem> :domain domain :components '(u) :patch->coefficients
	   `((:d-dimensional (,(f[x]->coefficient
				'MY-COEFFICIENT (lambda (x) (aref x 0)))))))))
    (plot problem :refinements 2 :coefficient 'MY-COEFFICIENT))
  )

;;; (test-coeffplot)
(fl.tests:adjoin-test 'test-coeffplot)