;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; dx.lisp - Interfacing with IBM's Data Explorer
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

(defvar *dx-pathname*
  (or (whereas ((dx-name cl-user::*dx-path*))
	(probe-file (pathname dx-name)))
      (fl.port:find-executable "dx"))
  "Pathname to DX.")

(defvar *dx-stream* nil
  "The current Data Explorer output stream.")

(defparameter *dx-toggle* 0
  "First, the filename is toggled between outputs for making writing to
files and reading from that file via dx simultaneously possible.  This
works, but it is not really safe.  Better would be waiting for dx to finish
output.  Furthermore, depending on *dx-toggle* the cache option is switched
on/off.  This is a trick to make dx redraw the picture.")

(defun ensure-dx-stream ()
  (setq *dx-stream*
	(if (and *dx-stream* (open-stream-p *dx-stream*))
	    *dx-stream*
	    (when *dx-pathname*
	      (whereas ((process
			 (fl.port:run-program
			  *dx-pathname*
			  `("-script" "-cache off" "-log on"
			    ,(concatenate 'string "-directory " *images-directory*))
			  :input :stream :output nil :wait nil)))
		(fl.port:process-input process)))))
  (unless *dx-stream*
    (format *error-output* "~&ENSURE-DX-STREAM: could not open stream.~%"))
  *dx-stream*)

(defmethod graphic-stream ((program (eql :dx)))
  (ensure-dx-stream))

#+(or)  ; (ensure-dx-stream)
(progn
  (close *dx-stream*)
  (setq *dx-stream* nil))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Communication with graphic servers
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmethod graphic-file-name (object (program (eql :dx)) &key &allow-other-keys)
  (declare (ignore object))
  (format nil "output-~D.dx" *dx-toggle*))

(defmethod send-graphic-commands (stream object (program (eql :dx)) &rest paras
				  &key (plot t) dimension (background :black)
				  (resolution 480) (width 480) (height 480)
				  format filename &allow-other-keys)
  (let* ((file-name (apply #'graphic-file-name object program paras))
	 (long-file-name (concatenate 'string (namestring *images-pathname*) file-name)))
    (format stream "data = Import(~S);~%" long-file-name)
    (format stream "data = Options(data, \"cache\", ~D);~%" *dx-toggle*)
    (loop for script-command in (apply #'graphic-commands object program paras)
	  when script-command do
	  (format stream "~A~%" script-command))
    (format stream "camera = AutoCamera(image, direction=~S, background=~S, resolution=~D, aspect=1.0);~%"
	    (if (= dimension 3) "off-diagonal" "front")
	    (ecase background (:black "black") (:white "white"))
	    resolution)
    (format stream "image = Render(image, camera);~%")
    (ecase plot
      ((t :plot)
       (let ((control "where=SuperviseWindow(\"femlisp-image\",size=[~D,~D],visibility=~D);~%"))
	 (format stream control width height 2)
	 (format stream control width height 1))
       (format stream "Display (image, where=where);~%"))
      (:file
       (when format
	 (format stream "WriteImage (image,~S,~S);~%"
		 (concatenate 'string (namestring *images-pathname*)
			      (or filename file-name))
		 format))))
    ))

(defmethod graphic-output :after (object (program (eql :dx)) &key &allow-other-keys)
  (setq *dx-toggle* (if (zerop *dx-toggle*) 1 0)))



