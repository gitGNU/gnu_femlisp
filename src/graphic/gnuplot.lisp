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
	    (whereas ((process
		       (ext::run-program *gnuplot-path* '() :input :stream
					 :output nil :wait nil)))
	      (ext::process-input process))))
  (unless *gnuplot-stream*
    (format *error-output*
	    "Could not open stream to Gnuplot.  Please ensure that the
special variable CL-USER::*GNUPLOT-PATH* has a reasonable value.  Usually,
this is set in femlisp:src;femlisp-config.lisp."))
  *gnuplot-stream*)

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

