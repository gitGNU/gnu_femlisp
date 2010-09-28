;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; function-plot.lisp - plotting of coefficient functions
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

(defmethod graphic-commands ((f <function>) (program (eql :dx)) &rest rest)
  (apply #'fl.graphic::dx-commands-data rest))

(defmethod plot ((f <function>) &rest rest &key cells mesh domain parametric
		 (depth 0) (key #'identity) (refinements 0) &allow-other-keys)
  "Plots @arg{function} on the given cell list @arg{cells}.  If @arg{cells}
is empty, the highest-dimensional cells of @arg{mesh} are used.  If this is
NIL, then a temporary mesh on @arg{domain} is creatend and refined up to
level @arg{refinements}.  Each cell is additionally refined @arg{depth}
times."
  (ensure cells
	  (cells-of-highest-dim
	   (or mesh (uniformly-refined-mesh domain refinements :parametric parametric))))
  (apply #'call-next-method f
	 :dimension (domain-dimension f)
	 :rank (if (= (image-dimension f) 1) 1 2)
	 :shape (image-dimension f)
	 :cells cells
	 :cell->values
	 #'(lambda (cell)
	     (map 'vector
		  #'(lambda (lcoords)
		      (funcall key (evaluate f (local->global cell lcoords))))
		  (refcell-refinement-vertex-positions cell depth)))
	 rest))

;;; Testing

(defun test-function-plot ()
  (let ((domain (n-cube-domain 2)))
    (plot (make-instance
	   '<special-function> :evaluator
	   #'(lambda (x)
	       (let ((x (aref x 0))
		     (y (aref x 1)))
		 (if (or (> (abs (- 0.5 x)) 0.4)
			 (> (abs (- 0.5 y)) 0.4))
		     0.0
		     (let ((x (/ (- x 0.1) 0.8))
			   (y (/ (- y 0.1) 0.8)))
		       #I(x*(1-x)*y*(1-y)*sin(10*pi*x)*sin(8*pi*y))))))
	   :domain-dimension (dimension domain)
	   :image-dimension 1)
	  :domain domain
	  :refinements 4 :depth 2))
  )

;;; (test-function-plot)
(fl.tests::adjoin-test 'test-function-plot)