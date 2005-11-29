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
  (or (aand fl.start::*dx-path* (probe-file (pathname it)))
      (fl.port:find-executable "dx"))
  "Pathname of the @program{DX} binary.")

(defvar *dx-process* nil
  "The current @program{dx} process.")

(defvar *dx-file*
 (make-pathname :name "output.dx"
		:directory (pathname-directory *images-pathname*))
  "The output file for @program{dx}.")

;;; something went wrong last time, because this file is still there
(when (probe-file *dx-file*)
  (delete-file *dx-file*)
  (warn "Deleted DX image file."))

(defun ensure-dx-process ()
  (when *dx-process*  ; (eq (fl.port:process-status *dx-process*) :running) ?
    (return-from ensure-dx-process *dx-process*))
  ;;; execute it within the images directory
  (setq *dx-process*
	(when *dx-pathname*
	  (fl.port:run-program
	   *dx-pathname* '("-script" "-cache" "off" "-log" "on") :wait nil
	   :input :stream :output (dbg-when :graphic *trace-output*)
	   :directory *images-pathname*)))
  (unless *dx-process*
    (warn "~&ENSURE-DX-PROCESS: could not start DX.~%"))
  *dx-process*)

(defmethod graphic-input-stream ((program (eql :dx)))
  (fl.port:process-input (ensure-dx-process)))

(defun dx-input-stream () (fl.port:process-input (ensure-dx-process)))
(defun dx-output-stream () (fl.port:process-output (ensure-dx-process)))

#+(or)  ; (ensure-dx-process)
(when *dx-process*
  (fl.port:process-close *dx-process*)
  (setq *dx-process* nil))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Communication with graphic servers
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defparameter *dx-toggle* 0
  "Depending on @var{*dx-toggle*} the cache option is switched on/off.
This is a trick to make @arg{dx} redraw the picture.")

(defmethod graphic-file-name (object (program (eql :dx)) &key &allow-other-keys)
  (declare (ignore object))
  *dx-file*)

(defmethod graphic-output :after (object (program (eql :dx)) &key &allow-other-keys)
  (declare (ignore object))
  (setq *dx-toggle* (if (zerop *dx-toggle*) 1 0)))

(defmethod send-graphic-commands (object (program (eql :dx)) &rest paras
				  &key (plot t) dimension (background :black)
				  (resolution 480) (width 480) (height 480)
				  format filename &allow-other-keys)
  (whereas ((stream (dx-input-stream)))
    (when (dbg-p :graphic)
      (setq stream (make-broadcast-stream stream *trace-output*)))
    (dbg :graphic "~A" *dx-file*)
    (let ((dx-file (namestring (probe-file (graphic-file-name object :dx)))))
      (format stream "data = Import(~S);~%" dx-file)
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
		   (concatenate 'string (namestring *images-pathname*) filename)
		   format))))
      (format stream "System(\"mv -f ~A ~A\");~%"
	      dx-file (concatenate 'string dx-file ".bak"))
      (force-output stream)
      )))

