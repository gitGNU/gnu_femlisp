;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; plot-gnuplot.lisp - Plotting with Gnuplot
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
;;; Default graphic commands
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmethod graphic-commands (object (program (eql :gnuplot)) &rest paras)
  "Default gnuplot plotting command."
  (list (format nil "plot ~S title \"data\" with linespoints"
		(namestring (apply #'graphic-file-name object program paras)))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Plotting of ordinary graphs (list of 2D-points)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun ensure-polygons (polygons)
  (if (listp (car polygons))
      polygons
      (list polygons)))

(defmethod graphic-write-data (stream (polygons list) (program (eql :gnuplot))
				      &key &allow-other-keys)
  "Can handle either a single list of points, or a single list of the form
\(string . points) or a list of such lists which is interpreted as multiple
graphs."
  (dovec (polygon (ensure-polygons polygons))
    (dovec (point (if (stringp (car polygon)) (cdr polygon) polygon))
      (format stream "~G ~G~%" (elt point 0) (elt point 1)))
    (format stream "~%~%")))

(defmethod graphic-commands ((polygons list) (program (eql :gnuplot))
			     &key (linewidth 1) &allow-other-keys)
  (let ((gnuplot-file (namestring (truename (graphic-file-name t :gnuplot)))))
    (list
     (apply #'concatenate 'string
	    "plot "
	    (loop for polygon in (ensure-polygons polygons)
		  and i from 0 
		  unless (zerop i) collect ", " end
		  collect
		  (format nil "~S index ~D title ~S with lines lw ~D"
			  gnuplot-file i
			  (if (stringp (car polygon))
			      (car polygon)
			      (format nil "data-~D" i))
			  linewidth))))))

(defmethod plot ((polygons list) &rest rest &key &allow-other-keys)
  (apply #'graphic-output polygons :gnuplot rest))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Plotting of finite element functions (1d)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


(defmethod graphic-write-data (stream object (program (eql :gnuplot))
				      &key cells cell->values (depth 0))
  "Writes data in gnuplot format to a stream."
  (declare (ignore object))
  (assert (= 1 (dimension (car cells))))
  (let* ((position-indices (compute-position-indices cells depth))
	 (position-array (position-array cells position-indices depth))
	 (values (and cell->values
		      (compute-all-position-values
		       cells position-indices depth cell->values))))
    (flet ((index->xpos (index) (aref (aref position-array index) 0)))
      (dolist (connection (sort (connections cells position-indices depth)
				#'<= :key (compose #'index->xpos #'car)))
	(dolist (index connection)
	  (format stream "~G ~G~%" (index->xpos index)
		  (if values (aref values index) 0.0)))))))

;;;; Testing: (test-plot-gnuplot)
(defun test-plot-gnuplot ()
  (let ((graph '(("graph-1" #(1.0 2.0) #(3.0 4.0)))))
    (plot graph :debug t))
  (let ((graph '("graph-2" #(1.0 3.0) #(4.0 2.0))))
    (plot graph :debug t))
  (plot
   (list
    (cons
     "{/Symbol \\266W}"
     (loop for phi from 0 upto (* 2 pi) by (* 0.005 pi)
	   collect (vector (cos phi) (sin phi))))
    (cons
     "{/Symbol \\266W^e}"
     (loop for phi from 0 upto (* 2 pi) by (* 0.005 pi)
	   for r = #I"1 + 0.1*sin(phi)^^2*sin(40*phi)"
	   collect (vector (* r (cos phi)) (* r (sin phi))))))
   :border nil :tics nil
   :terminal #-(or)"postscript eps enhanced color" #+(or) "epslatex color"
   :left -1.7 :right 1.7 :top 1.15 :bottom -1.15 :linewidth 3)
  )

(fl.tests:adjoin-test 'test-plot-gnuplot)

