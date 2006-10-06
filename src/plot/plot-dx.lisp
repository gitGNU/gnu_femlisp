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

(defun write-positions (stream position-array &optional (header t))
  "Write a header and all positions to the stream."
  (let ((dim (length (elt position-array 0))))
    (assert (every #'(lambda (pos) (= dim (length pos))) position-array))
    ;; write positions
    (when header
      (format stream "object 1 class array type float rank 1 shape ~D items ~D data follows~%"
	      dim (length position-array)))
    (loop for position across position-array do
	 (loop for coord across position do
	      (format stream "~G " (float coord 1.0))
	    finally (terpri stream)))))

(defparameter *dx-bug-transformation*
  (make-instance '<linear-function>
		 :A (let ((eps 1.0d-5))
		      (make-real-matrix `((,(cos eps) ,(- (sin eps)))
					  (,(sin eps) ,(cos eps))))))
  "Transformation which drops the ugly line occuring for me due to some bug
either in DX or Mesa.")

(defmethod graphic-write-data (stream object (program (eql :dx))
			       &key cells cell->values
			       (rank 0) shape
			       (depth 0) transformation)
  "Rather general plotting procedure writing data to a file.  Plots in
Gnuplot format in 1D, to DX format in 2D and 3D.  Can plot data either
discontinuous or continuous at the cell boundaries when coefficient
functions are involved.  cell->values has to be either nil --meaning no
data-- or some function mapping cells to a list of corner values."
  (declare (ignore object))
  (unless cells (return-from graphic-write-data))
  (let* ((dim (dimension (car cells)))
	 (position-indices (compute-position-indices cells depth))
	 (position-array (position-array cells position-indices depth transformation))
	 (values (and cell->values
		      (compute-all-position-values
		       cells position-indices depth cell->values)))
	 (simplices (remove-if-not #'simplex-p cells))
	 (other-cells (remove-if #'simplex-p cells)))
	 
    (write-positions stream position-array)
    
    ;; write simplices
    (when simplices
      (let ((connections (connections simplices position-indices depth)))
	(format stream "object 2 class array type int rank 1 shape ~D items ~D data follows~%"
		(1+ dim) (length connections))
	(loop for connection in connections do
	      (format stream "~{~D~^ ~}~%" connection))
	(format stream "attribute \"element type\" string \"~A\"~%"
		(case dim (1 "lines") (2 "triangles") (3 "tetrahedra")))
	(format stream "attribute \"ref\" string \"positions\"~%")))
       
    ;; write other-cells (these are written as cubes)
    (when other-cells
      (let ((connections (connections other-cells position-indices depth)))
	(format stream "object 3 class array type int rank 1 shape ~D items ~D data follows~%"
		(expt 2 dim) (length connections))
	(loop for connection in connections do
	      (format stream "~{~D~^ ~}~%" connection))
	(format stream "attribute \"element type\" string \"~A\"~%"
		(ecase dim (2 "quads") (3 "cubes")))
	(format stream "attribute \"ref\" string \"positions\"~%")))
    
    ;; write data
    (when values
      (format stream "object 4 class array type float rank ~D" rank)
      (when shape (format stream " shape ~D" shape))
      (format stream " items ~D data follows~%" (length values))
      (loop for value across values do
	    (setq value (etypecase value
			  (number (list value))
			  (vector (coerce value 'list))))
	    (loop for num in value do
		  (format stream "~10,8,,,,,'eG " (coerce num 'single-float))
		  finally (terpri stream)))
      (format stream "attribute \"dep\" string \"positions\"~%"))

    ;; epilogue
    (when simplices
      (format stream "object \"simplex-part\" class field~%")
      (format stream "component \"positions\" value 1~%")
      (format stream "component \"connections\" value 2~%")
      (when values
	(format stream "component \"data\" value 4~%")))
    (when other-cells
      (format stream "object \"quad-part\" class field~%")
      (format stream "component \"positions\" value 1~%")
      (format stream "component \"connections\" value 3~%")
      (when values
	(format stream "component \"data\" value 4~%")))
    (when (and simplices other-cells)
      (format stream "object \"both-grids\" class group~%")
      (format stream "member \"simplex-part\" value \"simplex-part\"~%")
      (format stream "member \"quad-part\" value \"quad-part\"~%"))
    (format stream "end~%")))

