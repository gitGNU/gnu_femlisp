;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; plot-dx.lisp - DX plotting of meshes and finite element functions
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

(defun dx-position-header (dim nr-positions)
   (format nil "object 1 class array type float rank 1 shape ~D items ~D data follows~%"
           dim nr-positions))

(defparameter *dx-bug-transformation*
  (make-instance '<linear-function>
		 :A (let ((eps 1.0d-5))
		      (make-real-matrix `((,(cos eps) ,(- (sin eps)))
					  (,(sin eps) ,(cos eps))))))
  "Transformation which drops the ugly line occuring for me due to some bug
either in DX or Mesa.")

(defparameter *shapes*
  '((2 (1) 2 "lines")
    (3 (2) 3 "triangles")
    (4 (1 1) 4 "quads")
    (5 (3) 4 "tetrahedra")
    (6 (2 1) 8 "cubes")
    (7 (1 2) 8 "cubes")
    (8 (1 1 1) 8 "cubes"))
  "The shapes which are handled by DX plotting.")

(defmethod graphic-write-data (stream object (program (eql :dx))
			       &key cells cell->values
			       (rank 0) shape
			       (depth 0) transformation)
  "Rather general plotting procedure writing data to a file.  Plots in DX format
for 1D, 2D, and 3D.  Can plot data either discontinuous or continuous at the
cell boundaries when coefficient functions are involved.  cell->values has to be
either nil --meaning no data-- or some function mapping cells to a list of
corner values."
  (declare (ignore object))
  (unless cells (return-from graphic-write-data))
  (let* ((position-indices (compute-position-indices cells depth))
	 (position-array (position-array cells position-indices depth transformation))
	 (values (and cell->values
		      (compute-all-position-values
		       cells position-indices depth cell->values)))
         (groups ()))
	 
    (write-positions stream position-array 'dx-position-header)
    
    (loop for (object shape nr-vertices element-type) in *shapes*
       for current-cells = (filter shape cells :key #'factor-dimensions :test #'equalp)
       for connections = (connections current-cells position-indices depth)
       when connections do
         (assert (same-p connections :key #'length))
         (push (list object element-type) groups)
         (format stream "object ~D class array type int rank 1 shape ~D items ~D data follows~%"
                 object nr-vertices (length connections))
         (loop for connection in connections do
              (format stream "~{~D~^ ~}~%" (rest connection)))
         (format stream "attribute \"element type\" string \"~A\"~%" element-type)
         (format stream "attribute \"ref\" string \"positions\"~%"))
    ;; we want groups from lines to cubes
    (setf groups (reverse groups))
    
    ;; write data
    (when values
      (format stream "object 9 class array type float rank ~D" rank)
      (when shape (format stream " shape ~D" shape))
      (format stream " items ~D data follows~%" (length values))
      (dovec (value values)
	(dovec (num value)
	  (format stream "~10,8,,,,,'eG " (coerce num 'single-float)))
	(terpri stream))
      (format stream "attribute \"dep\" string \"positions\"~%"))

    ;; epilogue
    (loop for (object element-type) in groups do
         (format stream "object \"~A-part\" class field~%" element-type)
         (format stream "component \"positions\" value 1~%")
         (format stream "component \"connections\" value ~D~%" object)
         (when values
           (format stream "component \"data\" value 9~%")))
    (when groups
      (format stream "object \"all-grids\" class group~%")
      (mapc (_ (format stream "member \"~A-part\" value \"~A-part\"~%" _ _))
            (mapcar #'second groups)))
    (format stream "end~%")))

