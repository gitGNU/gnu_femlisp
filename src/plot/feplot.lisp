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
	    (setf (matrix-ref result i j)
		  (evaluate shape (vertex-position vertex)))))
    result))

(defmethod graphic-commands ((asv <ansatz-space-vector>) (program (eql :dx))
			     &key (foreground :white) &allow-other-keys)
  (let ((axis-color (ecase foreground (:black "black")(:white "white")))
	(graph-color (ecase foreground (:black "black")(:white "yellow"))))
    (case (dimension (mesh asv))
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
      (3 (list
	  "connections = ShowConnections(data);"
	  "tubes = Tube(connections, 0.01);"
	  "image = AutoColor(tubes);"
	  ;;(format nil "camera = AutoCamera(image, \"off diagonal\", background=~S);" background-string)
	  ;;"image = Render(image, camera);"
	  )))))

(defmethod plot ((asv <ansatz-space-vector>) &rest rest
		 &key depth component (index 0) (key #'identity)
		 &allow-other-keys)
  "Plots a certain component of asv.  Index is useful when asv consists of
several vectors."
  (let* ((fe-class (fe-class (ansatz-space asv)))
	 (order (discretization-order fe-class))
	 (mesh (mesh asv)))
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
	   :dimension (dimension mesh)
	   :cells (surface-cells-of-highest-dim mesh)
	   :depth depth
	   :cell->values
	   #'(lambda (cell)
	       (let* ((fe (get-fe fe-class (reference-cell cell)))
		      (local-vec (get-local-from-global-vec cell fe asv)))
		 (when (typep fe '<vector-fe>)
		   (setq fe (aref (components fe) component))
		   (setq local-vec (aref local-vec component)))
		 (map-matrix key (m* (local-evaluation-matrix fe depth)
				     (matrix-slice local-vec :from-col index :ncols 1)))))
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
