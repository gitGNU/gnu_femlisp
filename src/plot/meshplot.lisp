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

(defun plot-transformation (dim)
  "The default plot transformation is simply a projection on R^3."
  (when (> dim 3)
    (make-instance
     '<linear-function> :A (eye 3 dim))))

(defmethod plot ((skel <skeleton>) &rest rest
		 &key transformation &allow-other-keys)
  (let ((dim (embedded-dimension skel)))
    (apply #'call-next-method skel
	   :dimension (plot-dimension dim)
	   :transformation (or transformation (plot-transformation dim))
	   rest)))

(defmethod plot ((cell <cell>) &rest rest)
  (apply #'plot (skeleton cell) rest))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Plot with DX
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmethod graphic-commands ((skel <skeleton>) (program (eql :dx))
			     &rest rest)
  (apply #'fl.graphic::dx-commands-mesh rest))


(defmethod graphic-write-data (stream (skel <skeleton>) (program (eql :dx))
			       &key (cells nil cells-p) part transformation)
  "Plots a mesh. If provided, @arg{cells} should be 1-cells."
  (unless cells-p
    (setq cells (1d-surface-cells skel part)))
  (unless cells (return-from graphic-write-data))
  (assert (every #'(lambda (cell) (= (dimension cell) 1)) cells))
  (let* ((position-indices (compute-position-indices cells 0))
	 (position-array (position-array cells position-indices 0 transformation)))
    ;; write positions
    (dbg :plot "Positions:~%~A" position-array)
    (write-positions stream position-array 'dx-position-header)
    ;; write connections
    (format stream "object 2 class array type int rank 1 shape 2 items ~D data follows~%"
	    (length cells))
    (let* ((bdry (boundary (n-cube 1)))
	   (from (aref bdry 1)) (to (aref bdry 0)))
      (dolist (cell cells)
	(format stream "~D ~D~%"
		(gethash (cons cell from) position-indices)
		(gethash (cons cell to) position-indices))))
    (format stream "~
attribute \"ref\" string \"positions\"~%~
attribute \"element type\" string \"lines\"~%~
object \"grid\" class field~%~
component \"positions\" value 1~%~
component \"connections\" value 2~%~
end~%")))

#| Test of 1d plotting with dx
dx -script

data = Import("output.dx");
data = Options(data, "mark", "circle");
xyplot = Plot(data, corners={[0.0,-0.1],[1.0,0.1]});
camera = AutoCamera(xyplot);
image = Render (xyplot, camera);
Display (image);

data = Import("output.dx");
connections = ShowConnections(data);
positions = ShowPositions(data);
image = Collect(positions,connections);
image = Options(image, "cache", 0);
camera = AutoCamera(image,direction="front");
image = Render(image, camera);
Display (image);
|#

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Plot with VTK
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmethod graphic-write-data (stream (skel <skeleton>) (program (eql :vtk))
			       &key (cells nil cells-p) part transformation)
  "Plots a mesh. If provided, @arg{cells} should be 1-cells."
  (unless cells-p
    (setq cells (1d-surface-cells skel part)))
  (unless cells (return-from graphic-write-data))
  (assert (every #'(lambda (cell) (= (dimension cell) 1)) cells))
  ;;; Not yet changed from DX output
  (let* ((position-indices (compute-position-indices cells 0))
	 (position-array (position-array cells position-indices 0 transformation)))
    ;; write positions
    (dbg :plot "Positions:~%~A" position-array)
    (write-positions stream position-array 'vtk-position-header)
    ;; write connections
    (let ((n (length cells)))
      (format stream "CELLS ~D ~D~%" n (* 3 n))
      (let* ((bdry (boundary (n-cube 1)))
             (from (aref bdry 1)) (to (aref bdry 0)))
        (dolist (cell cells)
          (format stream "2 ~D ~D~%"
                  (gethash (cons cell from) position-indices)
                  (gethash (cons cell to) position-indices)))))
    (let ((n (length position-array)))
      (format stream "POINT_DATA ~D~%SCALARS scalar_data float 1~%LOOKUP_TABLE default~%~{~A~%~}"
            n (make-list n :initial-element 0)))))

;;;; Testing

(defun test-meshplot ()
  (plot (n-cube 4) :transformation
	(make-instance '<linear-function>
		       :A #m((1.0 0.0   0.0  0.0)
			     (0.0 0.8  -0.6  0.0)
			     (0.0 0.56  0.57 0.6))))
  (plot (skeleton-boundary (n-ball-domain 3)))
  )

(fl.tests::adjoin-test 'test-meshplot)