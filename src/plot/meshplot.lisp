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

(in-package :plot)

(defmethod graphic-commands ((skel <skeleton>) (program (eql :dx))
			     &key (tubes t) (glyphs t) background
			     &allow-other-keys)
  (let ((dim (dimension skel)))
    (list
     "connections = ShowConnections(data);"
     (if tubes
	 (format nil "tubes = Tube(connections,~F);"
		 (if (numberp tubes)
		     tubes
		     (ecase dim ((1 2) 0.01) (3 0.01))))
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
	     (ecase dim ((1 2) "") (3 ", \"off diagonal\"")))
     )))

(defmethod plot ((skel <skeleton>) &rest rest &key depth &allow-other-keys)
  (when depth
    (error "The depth option is not yet allowed for mesh plotting."))
  (apply #'graphic-output skel :dx
	 :cells (if (typep skel '<hierarchical-mesh>)
		    (surface-cells-of-highest-dim skel)
		    (cells-of-highest-dim skel))
	 :dimension (manifold-dimension skel)
	 rest))

(defmethod plot ((cell <cell>) &key &allow-other-keys)
  (plot (skeleton cell)))
#| Test of 1d plotting with dx
dx -script

data = Import("image-mesh-1d.dx");
data = Options(data, "mark", "circle");
xyplot = Plot(data, corners={[0.0,-0.1],[1.0,0.1]});
camera = AutoCamera(xyplot);
image = Render (xyplot, camera);
Display (image);

connections = ShowConnections(data);
positions = ShowPositions(data);
colored = AutoColor(data);
surface = Isosurface(data, number=20);
surface = RubberSheet(surface);
image = Collect(surface,colored);
image = Options(image, "cache", 0);
camera = AutoCamera(image,direction="front");
image = Render(image, camera);
Display (image);
|#
