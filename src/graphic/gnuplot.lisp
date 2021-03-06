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

(in-package :fl.graphic)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Routines for establishing the communication line
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defvar *gnuplot-name*
  "gnuplot"
  "Pathname of the @program{Gnuplot} binary.")

(defvar *gnuplot-process* nil
  "The current @program{Gnuplot} process.")

(defun ensure-gnuplot-process ()
  (when *gnuplot-process* ; and: (eq (fl.port:process-status *gnuplot-process*) :running) ?
    (return-from ensure-gnuplot-process *gnuplot-process*))
  (setq *gnuplot-process*
	(when *gnuplot-name*
	  (fl.port:run-program-report-errors
	   *gnuplot-name* '() :wait nil
	   :input :stream :output :stream
	   :directory (images-pathname))))
  (unless *gnuplot-process*
    (format *error-output* "~&ENSURE-GNUPLOT-PROCESS: could not start GNUPLOT.~%"))
  *gnuplot-process*)

(defmethod graphic-input-stream ((program (eql :gnuplot)))
  (fl.port:process-input (ensure-gnuplot-process)))

(defun gnuplot-input-stream ()
  (whereas ((process (ensure-gnuplot-process)))
    (fl.port:process-input process)))
(defun gnuplot-output-stream ()
  (whereas ((process (ensure-gnuplot-process)))
    (fl.port:process-output process)))

#+(or)  ; (ensure-gnuplot-process)
(when *gnuplot-process*
  (fl.port:process-close *gnuplot-process*)
  (setq *gnuplot-process* nil))

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

(defmethod graphic-file-name (object (program (eql :gnuplot)) &key &allow-other-keys)
  "Returns the output file for @program{Gnuplot}."
  (declare (ignore object))
  (make-pathname :name "output.gnuplot"
                 :directory (pathname-directory (images-pathname))))

(defun wait-for-gnuplot ()
  "Does not work - somehow Gnuplot does not print anything to its output
stream in contrast to DX."
  (whereas ((stream (gnuplot-output-stream)))
    (dbg :graphic "Waiting for Gnuplot:")
    (loop for line = (prog1 (read-line stream nil nil)
                       (dbg :graphic "Waiting... Read: ~A~%" line))
       until (search "Femlisp request processed" line)
       finally (break))))

(defmethod send-graphic-commands (object (program (eql :gnuplot)) &rest paras
				  &key left right top bottom
				  (border t) (tics t) (terminal "x11")
				  (output "gnuplot.ps") &allow-other-keys)
  (let ((stream (if (dbg-p :graphic)
		    (make-broadcast-stream (gnuplot-input-stream) *trace-output*)
		    (gnuplot-input-stream))))
;;    (format stream "set size 1.0,1.0;~%")
    (format stream "set size square;~%")
    (if (and left right)
	(format stream "set xrange [~a:~a];~%" left right)
	(format stream "set autoscale x;~%"))
    (if (and bottom top)
	(format stream "set yrange [~a:~a];~%" bottom top)
	(format stream "set autoscale y;~%"))
    (if border
	(format stream "set border;~%")
	(format stream "set noborder;~%"))
    (if tics
	(format stream "set xtics~%set ytics;~%")
	(format stream "set noxtics~%set noytics;~%"))
    (format stream "set terminal ~A;~%" terminal)
    (let ((output-file
           (concatenate 'string (namestring (images-pathname)) output)))
      (format stream "set output ~S;~%" output-file))
    (loop for script-command in (apply #'graphic-commands object program paras) do
         (format stream "~A;~%" script-command))
    (format stream "print \"\\nFemlisp request processed\\n\";~%")
    (finish-output stream)
    ;; (wait-for-gnuplot)          ; does not work maybe due to readline
    ))

