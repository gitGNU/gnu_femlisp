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

(in-package :plot)

(defmethod graphic-commands ((problem <problem>) (program (eql :dx))
			     &key (foreground :white) &allow-other-keys)
  (let ((axis-color (ecase foreground (:black "black")(:white "white")))
	(graph-color (ecase foreground (:black "black")(:white "yellow"))))
    (case (dimension (domain problem))
      (1 (list
	  "data = Options(data, \"mark\", \"circle\");"
	  (format nil "data = Color(data,~S);" graph-color)
	  (format nil "image = Plot(data, colors=~S);" axis-color)
	  ;;"xyplot = Plot(data);"
	  ;;"camera = AutoCamera(xyplot);"
	  ;;"image = Render (xyplot, camera);"
	  ))
      (2 (list
	  "colored = AutoColor(data);"
	  "surface = Isosurface(data, number=20);"
	  "image = Collect(surface,colored);"
	  ;;"camera = AutoCamera(image);"
	  ;;"image = Render(image, camera);"
	  ))
      (3 (list
	  "connections = ShowConnections(data);"
	  "tubes = Tube(connections, 0.01);"
	  "image = AutoColor(tubes);"
	  ;;"camera = AutoCamera(tubes, \"off diagonal\");"
	  ;;"image = Render(tubes, camera);"
	  )))))

(defmethod plot ((problem <problem>) &rest rest
		 &key (depth 0) (key #'identity) refinements parametric coefficient &allow-other-keys)
  "Plots a coefficient function for problem.  It generates a temporary
mesh, refines it as often as given in the keyword parameter refinements and
plots the coefficient function on this mesh where each cell is resolved as
given by the additional parameter depth."
  (let* ((mesh (uniformly-refined-mesh (domain problem)	refinements
				       :parametric parametric))
	 (cells (cells-of-highest-dim mesh)))
    (apply #'graphic-output problem :dx
	   :dimension (dimension mesh)
	   :cells cells
	   :cell->values 
	   #'(lambda (cell)
	       (let* ((coeff-func (getf (coefficients-of-cell cell mesh problem)
					coefficient))
		      (local-vertices (refcell-refinement-vertices
				       (reference-cell cell) depth))
		      (values (make-double-vec (length local-vertices))))
		 (dotimes (i (length local-vertices))
		   (let* ((lcoords (vertex-position (aref local-vertices i)))
			  (ci (list :local lcoords :global (local->global cell lcoords))))
		      (setf (aref values i)
			    (funcall key (evaluate coeff-func ci)))))
		 values))
	   rest)))

;;; Testing: (test-coeffplot)

(defun test-coeffplot ()
  (let* ((dim 1) (domain (n-cell-domain dim))
	 (my-problem
	  (make-instance
	   '<problem> :domain domain :patch->coefficients
	   (constantly (list 'MY-COEFFICIENT
			     (function->coefficient #'(lambda (x) (aref x 0))))))))
    (plot my-problem :refinements 0 :coefficient 'MY-COEFFICIENT))
  )

(tests::adjoin-femlisp-test 'test-coeffplot)