;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; gnuplot.lisp - Interfacing to gnuplot
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

(in-package :graphics)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Routines for establishing the communication line
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defvar *gnuplot-stream* nil
  "Stream to gnuplot.  Should perhaps be coalesced with the CLOCC
version.")

(defun ensure-gnuplot-stream ()
  (setq *gnuplot-stream*
	(if (and *gnuplot-stream* (open-stream-p *gnuplot-stream*))
	    *gnuplot-stream*
	    (ext::process-input
	     (ext::run-program *gnuplot-path* '() :input :stream
			       :output nil :wait nil)))))

(defmethod graphic-stream ((program (eql :gnuplot)))
  (ensure-gnuplot-stream))


#+(or)  ; (ensure-gnuplot-stream)
(progn
  (close *gnuplot-stream*)
  (setq *gnuplot-stream* nil))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Gnuplot enhancements
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun gnuplot-tics (filename)
  (with-open-file (stream filename :direction :input)
    (loop for line = (read-line stream nil 'eof)
	  until (eq line 'eof)
	  collect (read-from-string line))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Plotting a fe function with gnuplot
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defparameter *gnuplot-toggle* 0
  "Depending on *gnuplot-toggle* we determine the filename.  This is a
trick to avoid writing the file while our gnuplot job is still reading.")

(defmethod graphic-file-name (object (program (eql :gnuplot))
			      &key &allow-other-keys)
  (format nil "output~D.dat" *gnuplot-toggle*))

(defmethod send-graphic-commands (stream object (program (eql :gnuplot)) &rest paras
					 &key left right top bottom
					 (border t) (tics t) (terminal "x11")
					 (output "gnuplot.out") &allow-other-keys)
  (if (and left right)
      (format stream "set xrange [~a:~a]~%" left right)
      (format stream "set autoscale x~%"))
  (if (and bottom top)
      (format stream "set yrange [~a:~a]~%" bottom top)
      (format stream "set autoscale y~%"))
  (if border
      (format stream "set border~%")
      (format stream "set noborder~%"))
  (if tics
      (format stream "set xtics~%set ytics~%")
      (format stream "set noxtics~%set noytics~%"))
  (format stream "set size 1.0,1.0; set terminal ~A~%" terminal)
  (format stream "set output ~S~%" (concatenate 'string *images-directory* output))
  (loop for script-command in (apply #'graphic-commands object program paras) do
	(format stream "~A~%" script-command))
  (force-output stream))

(defmethod graphic-output :after (object (program (eql :gnuplot)) &key &allow-other-keys)
  (setq *gnuplot-toggle* (if (zerop *gnuplot-toggle*) 1 0)))

#+(or) ; inactive
(defmethod graphic-end (object (program (eql :gnuplot)) &key &allow-other-keys)
  (setq *gnuplot-toggle* (if (zerop *gnuplot-toggle*) 1 0)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Old code, maybe still useful
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

#|
(defun gnuplot-action (action &key (plot t) left right bottom top &allow-other-keys)
  "Performs a certain gnuplot action."
  (let ((stream (if (eq plot :file)
		    (open *gnuplot-file* :direction :output :if-exists :supersede)
		    (ensure-gnuplot-stream))))
    (unwind-protect
	 (progn
	   (if (and left right)
	       (format stream "set xrange [~a:~a]~%" left right)
	       (format stream "set autoscale x~%"))
	   (if (and bottom top)
	       (format stream "set yrange [~a:~a]~%" bottom top)
	       (format stream "set autoscale y~%"))
	   (funcall action stream))
      (ecase plot
	((t :plot) (format *standard-output* "Done plotting.~%"))
	(:wait
	 (fresh-line *terminal-io*)
	 (princ "Press <enter> to continue..." *terminal-io*)
	 (force-output *terminal-io*) (read-line *terminal-io* nil nil))
	(:print (format *standard-output* "Sent the plot to `~a'.~%" cllib::*gnuplot-printer*)
		(format stream "set output~%"))
	(:file (format *standard-output* "Wrote `~a'.~%Type \"load '~a'\" at the gnuplot prompt.~%"
		       (extensions:unix-namestring *gnuplot-file*)
		       (extensions:unix-namestring *gnuplot-file*))))
      (force-output stream))))

(defun gnuplot-vector-field (vec-field &rest rest &key (title "vf")
			     (left -1.0) (right 1.0) (bottom -0.8) (top 0.8)
			     n nx ny scale &allow-other-keys)
  "Plot the given vec-field function."
  (apply #'gnuplot-action
	 #'(lambda (stream)
	     (format stream "plot '-' using 1:2:3:4 title \"~a\" with vector~%~%" title)
	     (loop with dx = (/ (- right left) (or nx n))
		   and dy =  (/ (- top bottom) (or ny n))
		   for x = left then (+ x dx)
		   until (> x right) do
		   (loop for y = bottom then (+ y dy)
			 until (> y top) do
			 (let* ((vf (funcall vec-field (vector x y)))
				(vec (etypecase scale
				       (number (map (type-of vf) #'(lambda (x) (* scale x)) vf))
				       (function (funcall scale vf)))))
			   (format stream "~f ~f ~f ~f~%"
				   (- x (* 0.5 (aref vec 0))) (- y (* 0.5 (aref vec 1)))
				   (aref vec 0) (aref vec 1))))
		   finally (format stream "e~%")))
	 :left left :right right :bottom bottom :top top
	 rest))

(defun gnuplot-polygon (points &rest rest &key (title "data") &allow-other-keys)
  "Plot all points in list connected by lines."
  (apply #'gnuplot-action
	 #'(lambda (stream)
	     (format stream "plot '-' using 1:2 title \"~a\" with lines~%~%" title)
	     (dolist (point points)
	       (format stream "~f ~f~%" (vec-ref point 0) (vec-ref point 1)))
	     (format stream "e~%"))
	 rest))

(defun gnuplot-polygons (polygon-list &rest rest)
  "Plot all points in list connected by lines."
  (apply #'gnuplot-action
	 #'(lambda (stream)
	     (loop initially (format stream "plot ")
		   for polygon in polygon-list
		   for beginning = t then nil do
		   (unless beginning (format stream ", "))
		   (format stream "'-' title \"~a\" with lines" (car polygon))
		   finally (format stream "~%~%"))
	     (dolist (polygon polygon-list)
	       (dolist (point (cdr polygon))
		 (format stream "~f ~f~%" (vec-ref point 0) (vec-ref point 1)))
	       (format stream "e~%")))
	 rest))

(defun gnuplot-array (mat &rest rest &key (h 1.0) (title "data") &allow-other-keys)
  "Plot all points in the array mat connected by lines."
  (apply #'gnuplot-action
	 #'(lambda (stream)
	     (format stream "set hidden3d~%")
	     (format stream "splot '-' using 1:2:3 title \"~a\" with lines~%~%" title)
             (dotimes (row (nrows mat))
               (dotimes (col (ncols mat))
                 (format stream "~f ~f ~f~%" (* row h) (* col h) (mat-ref mat row col)))
               (terpri stream))
	     (format stream "e~%")
	     (format stream "set nohidden3d~%"))
	 rest))
|#
