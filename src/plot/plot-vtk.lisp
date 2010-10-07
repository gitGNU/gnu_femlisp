;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; plot-vtk.lisp - VTK plotting of meshes and finite element functions
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;
;;; Copyright (C) 2010 Nicolas Neuss, Karlsruhe Institute of Technology.
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
;;; MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
;;; IN NO EVENT SHALL THE AUTHOR, THE KARLSRUHE INSTITUTE OF TECHNOLOGY
;;; OR OTHER CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
;;; SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
;;; LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
;;; DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
;;; THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
;;; (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
;;; OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(in-package :fl.plot)

(defun vtk-position-header (dim n)
  (declare (ignore dim))
  (format nil "POINTS ~D float~%" n))

(defun write-vtk-positions (stream position-array)
  (format nil "POINTS ~D float~%" (length position-array))
  (loop for position across position-array do
	 (loop for i below 3 do
	      (format stream "~G "
                      (if (< i (length position))
                          (float (aref position i) 1.0)
                          0.0))
	    finally (terpri stream))))

(defmethod graphic-write-data (stream object (program (eql :vtk))
			       &key cells cell->values shape
			       (depth 0) transformation
                               &allow-other-keys)
  (declare (optimize debug))
  (declare (ignore object))
  (unless cells (return-from graphic-write-data))
  (let* ((position-indices (compute-position-indices cells depth))
	 (position-array (position-array cells position-indices depth transformation))
	 (values (and cell->values
		      (compute-all-position-values
		       cells position-indices depth cell->values))))

    ;; write datafile header
    (format stream "# vtk DataFile Version 2.0~%femlisp-vtk~%ASCII~%DATASET UNSTRUCTURED_GRID~%~%")

    ;; write point data
    (write-vtk-positions stream position-array)
    (terpri stream)

    (let ((connections (connections cells position-indices depth)))
      ;; corner info header
      (format stream "CELLS ~D ~D~%"
              (length connections)
              (reduce #'+ connections :key #'length))
      ;; corner info
      (loop for (cell . conn) in connections do
           (typecase cell
             (fl.mesh::<1-1-product-cell>
              (destructuring-bind (a b c d) conn
                (format stream "4 ~D ~D ~D ~D~%" a b d c)))
             (fl.mesh::<1-1-1-product-cell>
              (destructuring-bind (a b c d e f g h) conn
                (format stream "8 ~D ~D ~D ~D ~D ~D ~D ~D~%" a b d c e f h g)))
             (t (format stream "~D ~{~D~^ ~}~%" (length conn) conn))))
      
      ;; write cell type info
      (format stream "~%CELL_TYPES ~D~%" (length connections))
      (loop for conn in connections do
           (format stream "~D~%"
                   (etypecase (first conn)
                     (fl.mesh::<1-simplex> 13)
                     (fl.mesh::<2-simplex> (error "NYI"))
                     (fl.mesh::<3-simplex> 10)
                     (fl.mesh::<1-1-product-cell> 9)
                     (fl.mesh::<1-1-1-product-cell> 12)))))
    ;; write data
    (when values
      (format stream "POINT_DATA ~A~%SCALARS scalar_data float ~A~%LOOKUP_TABLE default~%"
              (reduce #'+ values :key #'nrows)
              (or shape 1))
      (dovec (value values)
	(dovec (num value)
	  (format stream "~10,8,,,,,'eG " (coerce num 'single-float)))
	(terpri stream)))
    ))

