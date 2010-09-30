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

(defmethod graphic-write-data (stream object (program (eql :vtk))
			       &key cells cell->values shape
			       (depth 0) transformation
                               &allow-other-keys)
  (declare (optimize debug))
  (declare (ignore object))
  (unless cells (return-from graphic-write-data))
  (let* ((dim (dimension (car cells)))
	 (position-indices (compute-position-indices cells depth))
	 (position-array (position-array cells position-indices depth transformation))
	 (values (and cell->values
		      (compute-all-position-values
		       cells position-indices depth cell->values))))

    ;; write datafile header
    (format stream "# vtk DataFile Version 2.0~%femlisp-vtk~%ASCII~%DATASET UNSTRUCTURED_GRID~%~%")

    ;; write point data
    (write-positions stream position-array 'vtk-position-header)
    (terpri stream)

    (assert (and (= dim 3) (every #'simplex-p cells)) ()
            "Only tetrahedra are implemented because I do not know other VTK cell types")

    (let ((connections (connections cells position-indices depth)))
      ;; corner info header
      (format stream "CELLS ~D ~D~%"
              (length connections)
              (reduce #'+ connections :key (_ (1+ (length _)))))
      ;; corner info
      (loop for conn in connections do
           (format stream "~D ~{~D~^ ~}~%" (length conn) conn))
      
      ;; write cell type info
      (format stream "~%CELL_TYPES ~D~%~{~D~%~}~%"
              (length connections)
              (mapcar (constantly 10) connections)))
    ;; Quad: 9, Hexahedron: 12, Wedge=Dreieck x Line: 13
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

