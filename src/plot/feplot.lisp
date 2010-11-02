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

(in-package :fl.plot)

(with-memoization (:id 'local-evaluation-matrix :test 'equal)
  (defun local-evaluation-matrix (fe depth &optional (type :solution))
    "Returns a local evaluation matrix for the refcell-refinement-vertices of
the given depth."
    (memoizing-let ((fe fe) (depth depth) (type type))
      (funcall (ecase type
		 (:solution 'ip-values)
		 (:gradient 'ip-gradients))
	       fe (refcell-refinement-vertex-positions
		   (reference-cell fe) depth)))))

(defmethod graphic-commands ((asv <ansatz-space-vector>) (program (eql :dx))
			     &key (foreground :white) dimension shape range
                             &allow-other-keys)
  (let ((axis-color (ecase foreground (:black "black")(:white "white")))
	(graph-color (ecase foreground (:black "black")(:white "white"))))
    (case dimension
      (1 (list
	  "data = Options(data, \"mark\", \"circle\");"
	  ;; workaround for a DX bug, see bug-2.dx 
	  (if (eq foreground :white)
	      (format nil "data = Color(data);")  ; ??
	      (format nil "data = Color(data,~S);" graph-color))
	  (format nil "image = Plot(data, colors=~S);" axis-color)
	  ;;"xyplot = Plot(data);"
	  ;;"camera = AutoCamera(xyplot);"
	  ;;"image = Render (xyplot, camera);"
	  ))
      (2 (cons
	  (if range
	      (format nil "colored = AutoColor(data, min=~A, max=~A);"
		      (first range) (second range))
	      "colored = AutoColor(data);")
	  (if (and shape (= shape 2))
	      (list
	       "samples = Sample(data, 400);"
	       (if range
		   (format nil "glyphs = AutoGlyph(samples, type=\"arrow\",min=~A,max=~A);"
			   (first range) (second range))
		   "glyphs = AutoGlyph(samples, type=\"arrow\");")
	       "image = Collect(colored,glyphs);")
	      (list
	       "surface = Isosurface(data, number=20);"
	       "image = Collect(surface,colored);"))))
      (3 (list
	  "connections = ShowConnections(data);"
	  "tubes = Tube(connections, 0.01);"
	  (if range
	      (format nil "image = AutoColor(tubes, min=~A, max=~A);"
		      (first range) (second range))
	      "image = AutoColor(tubes);")
	  ;;"camera = AutoCamera(tubes, \"off diagonal\");"
	  ;;"image = Render(tubes, camera);"
	  )))))

(defgeneric plot-cells (skel &key &allow-other-keys)
  (:method ((skel <skeleton>) &key &allow-other-keys)
    (cells-of-dim skel 1))
  (:method ((h-mesh <hierarchical-mesh>) &key part &allow-other-keys)
    (let ((dim (min 3 (if part
                          (dimension-of-part h-mesh part)
                          (dimension h-mesh)))))
      (find-cells
       (lambda (cell)
         (and (= (dimension cell) dim)
              (or (null part)
                  (eq (fl.mesh::get-patch-property cell h-mesh :part) part))))
       h-mesh :where :surface))))

(defmethod plot ((asv <ansatz-space-vector>) &rest rest
		 &key cells depth (index 0) (component 0) key transformation
		 (rank 0) shape part &allow-other-keys)
  "Plots a certain component of the ansatz space vector @arg{asv},
e.g. pressure for Navier-Stokes equations.  @arg{index} chooses one of
several solutions if @arg{asv} has multiplicity larger 1."
  (let* ((as (ansatz-space asv))
	 (mesh (mesh as))
	 (problem (problem as))
	 (dim (or (let ((fl.mesh::*check-well-defined-embedded-dimension* t))
                    (embedded-dimension mesh))
                  (dimension-of-part mesh part))))
    (ensure depth
	    #'(lambda (cell)
		(let* ((components (components-of-cell cell mesh problem))
		       (pos (component-position components component))
		       (order (discretization-order (get-fe as cell))))
		  (when (arrayp order)
		    (setq order (aref order pos)))
			(case order
			  ((0 1) 0)
			  ((2) 1)
			  ((3 4) 2)
			  (t 3)))))
    ;; (ensure rank (if (and component index) 0 1))
    (apply #'call-next-method asv
	   :depth depth
	   :cells (or cells (plot-cells mesh :part part))
	   :dimension (plot-dimension dim)
	   :transformation (or transformation (plot-transformation dim))
	   :rank rank :shape shape
	   :cell->values
	   #'(lambda (cell)
	       (let ((fe (get-fe as cell))
		     (local-vec (get-local-from-global-vec cell asv))
		     (depth (if (numberp depth) depth (funcall depth cell)))
		     ;;(m (nr-of-refinement-vertices cell depth))
		     )
		 (let ((values
			(multiple-evaluate-local-fe
			 local-vec (local-evaluation-matrix fe depth)))
		       (components (components-of-cell cell mesh problem)))
		   (multiple-value-bind (pos length)
		       (component-position components component)
		     (vector-map (lambda (vals)
				   (let ((value (extract-from vals pos length nil index)))
				     (if key (funcall key value) value)))
				 values)))))
	   ;; other paramters
	   rest
	   )))

#| Test of 1d plotting with dx
dx -script

data = Import("output.dx");
data = Options(data, "mark", "circle");
xyplot = Plot(data);
camera = AutoCamera(xyplot);
image = Render (xyplot, camera);
where = SuperviseWindow("femlisp-output", size=[480,480], visibility=2);
Display (image,where=where);
|#

#| Test of 2d vector plotting with dx
dx -script

data = Import("output.dx");
colored = AutoColor(data);
samples = Sample(data, 400);
glyphs = AutoGlyph(samples, type="arrow");
image = Collect(colored,glyphs);
image = Options(image, "rendering mode", "hardware");
where=SuperviseWindow("femlisp-image",size=[480,480],visibility=2);
where=SuperviseWindow("femlisp-image",size=[480,480],visibility=1);
camera = AutoCamera(image, direction="front", background="black", resolution=480, aspect=1.0);
Display (image, camera, where=where);

content = ReadImageWindow(where);
WriteImage(content, "test", "tiff");
|#
