;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; asaplot.lisp - plotting of ansatz-space-automorphism graphs
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


(defmethod graphic-commands ((asa <ansatz-space-automorphism>) (program (eql :dx))
			     &key &allow-other-keys)
  (list
   "connections = ShowConnections(data);"
   "tubes = Tube(connections, 0.01);"
   "glyphs = Glyph(data, scale=0.01);"
   "image = Collect(tubes,glyphs);"
   #+(or)(case (dimension (mesh asa))
     ((1 2) "camera = AutoCamera(image);")
     (3 "camera = AutoCamera(image, \"off diagonal\");")
     (t (error "Not implemented.  Maybe by projecting on 3D.")))
   ))

#|
dx -script

data = Import("image-data-1d.dx");
connections = ShowConnections(data);
glyphs = Glyph(data, scale=0.05);
tubes = Tube(connections, 0.02);
all = Collect(tubes,glyphs);
camera = AutoCamera(all);
image = Render(all, camera);
Display (image);
|#

(defmethod graphic-write-data (stream (asa <ansatz-space-automorphism>) (program (eql :dx))
				   &key dimension &allow-other-keys)
  (let ((node-indices (make-hash-table)))
    ;; write positions
    (format stream "object 1 class array type float rank 1 shape ~D items ~D data follows~%"
	    dimension (nr-nonempty-rows asa))
    (let ((count 0))
      (for-each-row-key
       #'(lambda (cell)
	   (setf (gethash cell node-indices) count)
	   (incf count)
	   (loop for coord across (midpoint cell) do
		 (format stream "~G " (float coord 1.0))
		 finally (terpri stream)))
       asa))
    ;; write connections
    (format stream "object 2 class array type int rank 1 shape ~D items ~D data follows~%"
	    2 (nr-of-entries asa))
    (for-each-key
     #'(lambda (row-key col-key)
	 (format stream "~D ~D "
		 (gethash row-key node-indices)
		 (gethash col-key node-indices)))
     asa)
    (format stream "~%attribute \"element type\" string \"lines\"~%")
    (format stream "attribute \"ref\" string \"positions\"~%"))
  ;; write epilogue
  (format stream "object \"line-part\" class field~%")
  (format stream "component \"positions\" value 1~%")
  (format stream "component \"connections\" value 2~%")
  (format stream "end~%"))

(defmethod plot ((asa <ansatz-space-automorphism>) &rest rest)
  (apply #'graphic-output asa :dx :dimension (dimension (mesh asa)) rest))

