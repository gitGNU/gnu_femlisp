;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; feplot.lisp - plotting of fe functions
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

(defmemo local-evaluation-matrix (fe depth)
  "Returns an evaluation matrix for the fe with the specified local
refinement depth."
  (let* ((vertices (refcell-refinement-vertices (reference-cell fe) depth))
	 (result (make-real-matrix (length vertices) (nr-of-dofs fe))))
    (loop for vertex across vertices and i from 0 do
      (loop for shape in (fe-basis fe) and j from 0 do
	    (setf (mref result i j)
		  (evaluate shape (vertex-position vertex)))))
    result))

(defmethod graphic-commands ((asv <ansatz-space-vector>) (program (eql :dx))
			     &key tubes (foreground :white) &allow-other-keys)
  (let ((axis-color (ecase foreground (:black "black")(:white "white")))
	(graph-color (ecase foreground (:black "black")(:white "yellow")))
	(dim (dimension (mesh asv))))
    (case dim
      (1 (list
	  "data = Options(data, \"mark\", \"circle\");"
	  (format nil "data = Color(data,~S);" graph-color)
	  (format nil "image = Plot(data, colors=~S);" axis-color)
	  ;;"image = Render (image, camera);"
	  ))
      (2 (list
	  "colored = AutoColor(data);"
	  "surface = Isosurface(data, number=20);"
	  "image = Collect(surface,colored);"
	  ;;(format nil "camera = AutoCamera(image, background=~S);" background-string)
	  ;;"image = Render(image, camera);"
	  ))
      (t (list
	  "connections = ShowConnections(data);"
	  (if tubes
	      (format nil "tubes = Tube(connections,~F);"
		      (if (numberp tubes)
			  tubes
			  (case dim ((1 2) 0.01) (t 0.01))))
	      "tubes = connections;")
	  "image = AutoColor(tubes);"
	  ;;(format nil "camera = AutoCamera(image, \"off diagonal\", background=~S);" background-string)
	  ;;"image = Render(image, camera);"
	  )))))

(defmethod plot ((asv <ansatz-space-vector>) &rest rest
		 &key depth component index key transformation
		 &allow-other-keys)
  "Plots a certain component of asv.  Index is useful when asv consists of
several vectors."
  (let* ((fe-class (fe-class (ansatz-space asv)))
	 (order (discretization-order fe-class))
	 (mesh (mesh asv))
	 (dim (manifold-dimension mesh)))
    (when (and (arrayp order) component)
      (setq order (aref order component)))
    (unless depth
      ;; even better might be a dependence on the mesh width
      (setq depth (case order
		    ((0 1) 0)
		    ((2) 1)
		    ((3 4) 2)
		    (t 3))))
    (apply #'graphic-output asv :dx
	   :depth depth
	   :cells (plot-cells mesh)
	   :dimension (plot-dimension dim)
	   :transformation (or transformation (plot-transformation dim))
	   :cell->values
	   #'(lambda (cell)
	       (let* ((fe (get-fe fe-class (reference-cell cell)))
		      (local-vec (get-local-from-global-vec cell fe asv))
		      (m (nr-of-refinement-vertices cell depth))
		      (components (etypecase fe
				    (<fe> (vector fe))
				    (<vector-fe> (components fe)))))
		 (unless (typep fe '<vector-fe>)
		   (setq component 0)
		   (setq local-vec (vector local-vec)))
		 (when (= (multiplicity asv) 1)
		   (setq index 0))
		 (assert (or (and key (not component))
			     (and component (not key))))
		 (let ((all (loop for k below (length components) collect
				  (m* (local-evaluation-matrix (aref components k) depth)
				      (aref local-vec k)))))
		   ;; we have to extract a one component value array
		   (if key
		       (let ((result (make-real-matrix m 1)))
			 (dotimes (i m result)
			   (setf (vref result i)
				 (funcall key (map 'vector (rcurry #'mref i index) all)))))
		       (matrix-slice (elt all component) :from-col index :ncols 1)))))
	   (sans rest :depth))))

#| Test of 1d plotting with dx
dx -script

data = Import("image-data-1d.dx");
data = Options(data, "mark", "circle");
xyplot = Plot(data);
camera = AutoCamera(xyplot);
image = Render (xyplot, camera);
where = SuperviseWindow("femlisp-output", size=[480,480], visibility=2);
Display (image,where=where);
|#
