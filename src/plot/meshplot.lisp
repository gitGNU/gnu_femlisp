;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; meshplot.lisp - plotting of meshes
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

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Helper routines
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun plot-dimension (dim)
  (min dim 3))

(defun plot-cells (skel)
  (let ((dim (dimension skel)))
    (if (typep skel '<hierarchical-mesh>)
	(surface-cells-of-dim skel (min 3 dim))
	(cells-of-dim skel (min 3 dim)))))

(defun plot-transformation (dim)
  "The default plot transformation is simply a projection on R^3."
  (when (> dim 3)
    (make-instance
     '<linear-function> :A (eye 3 dim))))

(defmethod graphic-commands ((skel <skeleton>) (program (eql :dx))
			     &key tubes (glyphs t) background
			     &allow-other-keys)
  (let ((dim (dimension skel)))
    (list
     "connections = ShowConnections(data);"
     (if tubes
	 (format nil "tubes = Tube(connections,~F);"
		 (if (numberp tubes)
		     tubes
		     (case dim ((1 2) 0.01) (t 0.01))))
	 "tubes = connections;")
     (when (eq background :white)
	 "tubes = Color(tubes, \"black\");")
     (when glyphs
       (format nil "glyphs = Glyph(data~A);"
	       (case dim (1 ",scale=0.01") (2 ",scale=0.01") (t ""))))
     (if glyphs
	 "image = Collect(tubes,glyphs);"
	 "image = tubes;")
     #+(or)
     (format nil "camera = AutoCamera(all~A);"
	     (case dim ((1 2) "") (t ", \"off diagonal\"")))
     )))

(defmethod plot ((skel <skeleton>) &rest rest &key transformation depth &allow-other-keys)
  (when depth
    (error "The depth option is not yet allowed for mesh plotting."))
  (let ((dim (manifold-dimension skel)))
    (apply #'graphic-output skel :dx
	   :cells (plot-cells skel)
	   :dimension (plot-dimension dim)
	   :transformation (or transformation (plot-transformation dim))
	   rest)))

(defmethod plot ((cell <cell>) &rest rest)
  (apply #'plot (skeleton cell) rest))

#| Test of 1d plotting with dx
dx -script

data = Import("output-0.dx");
data = Options(data, "mark", "circle");
xyplot = Plot(data, corners={[0.0,-0.1],[1.0,0.1]});
camera = AutoCamera(xyplot);
image = Render (xyplot, camera);
Display (image);

data = Import("output-1.dx");
connections = ShowConnections(data);
positions = ShowPositions(data);
image = Collect(positions,connections);
image = Options(image, "cache", 0);
camera = AutoCamera(image,direction="front");
image = Render(image, camera);
Display (image);
|#

(defun test-meshplot ()
  (plot (n-cube 4) :transformation
	(make-instance '<linear-function>
		       :A #m((1.0 0.0   0.0  0.0)
			     (0.0 0.8  -0.6  0.0)
			     (0.0 0.56  0.57 0.6))))
  )

(fl.tests::adjoin-test 'test-meshplot)