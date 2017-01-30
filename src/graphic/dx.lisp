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

(defvar *dx-name*
  "dx"
  "Name of the @program{DX} binary.")

(defvar *dx-process* nil
  "The current @program{dx} process.")

(defun ensure-dx-process ()
  (when *dx-process*  ; (eq (fl.port:process-status *dx-process*) :running) ?
    (return-from ensure-dx-process *dx-process*))
  ;; execute it within the images directory
  (setq *dx-process*
        (when *dx-name*
          (fl.port:run-program-report-errors
           *dx-name*
	   ;; -processors 1 is apparently needed with opendx-4.4.4
           ;; but unfortunately breaks CLISP
           '("-script" "-cache" "off" "-log" "on"
             "-processors" "1"
             #+os-windows "-native")
           :wait nil
           :input :stream :output :stream
           :directory (images-pathname))))
  (unless *dx-process*
    (warn "~&ENSURE-DX-PROCESS: could not start DX.~%"))
  ;; there is a bug in DX which leaves the DX window around
  ;; when not all children are killed which we do here
  ;; when exiting SBCL
  (fl.port:add-exit-hook 'dx-close)
  *dx-process*)

;;; #+sbcl (pushnew 'dx-close sb-ext:*exit-hooks*)

(defmethod graphic-input-stream ((program (eql :dx)))
  (fl.port:process-input (ensure-dx-process)))

(defun dx-input-stream ()
  (whereas ((process (ensure-dx-process)))
    (fl.port:process-input process)))
(defun dx-output-stream ()
  (whereas ((process (ensure-dx-process)))
    (fl.port:process-output process)))

(defun dx-close ()
  (whereas ((process *dx-process*))
    (fl.port::kill-process process)
    (fl.port:process-close process)
    (setq *dx-process* nil)))

;; (ensure-dx-process)
#+(or)  
(when *dx-process* (dx-close))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Communication with graphic servers
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defparameter *dx-toggle* 0
  "Depending on @var{*dx-toggle*} the cache option is switched on/off.
This is a trick to make @arg{dx} redraw the picture.")

(defmethod graphic-file-name (object (program (eql :dx)) &key &allow-other-keys)
  "Returns the output file for @program{dx}."
  (declare (ignore object))
  (make-pathname :name "output.dx"
                 :directory (pathname-directory (images-pathname))))

(defmethod graphic-output :after (object (program (eql :dx)) &key &allow-other-keys)
  (declare (ignore object))
  (setq *dx-toggle* (if (zerop *dx-toggle*) 1 0)))

(defun wait-for-dx ()
  (whereas ((stream (dx-output-stream)))
    (loop for line = (prog1 (read-line stream)
                       (when (dbg-p :graphic)
                         (format *trace-output* "Waiting... Read: ~A~%" line)
                         (force-output *trace-output*)))
       until (search "0:  ECHO:  Femlisp request processed" line))))

(defvar *dx-bug-workaround* nil
  "If T switches on hardware rendering.  This variant is problematic,
although it gets rid of black lines in dx pictures in some situations.
E.g. it does not do xy-graphs correctly and fails for @lisp{(plot (n-cube
1))}.  It can also kill the Xwindows interface on some computers.")

(defvar *show-dx-window* nil
  "Show the DX window when something is plotted.  This may be useful on
Laptops when the window is hidden.")

(defun handshake-with-dx ()
  (whereas ((stream (dx-input-stream)))
    (format stream "Echo (\"Femlisp request processed\");~%")
    (force-output stream)
    (wait-for-dx)))

(defmethod send-graphic-commands (object (program (eql :dx)) &rest paras
				  &key dimension (background :black)
				  (resolution 480) (width 480) (height 480)
				  format filename (window "femlisp-image") &allow-other-keys)
  (whereas ((stream (dx-input-stream)))
    (when (dbg-p :graphic)
      (setq stream (make-broadcast-stream stream *trace-output*)))
    (let ((dx-file (graphic-file-name object :dx :system-p t)))
      (dbg :graphic "File: ~A~%" dx-file)
      (format stream "data = Import(~S);~%" dx-file)
      (format stream "data = Options(data, \"cache\", ~D);~%" *dx-toggle*)
      (loop for script-command in (apply #'graphic-commands object program paras)
	    when script-command do
	    (format stream "~A~%" script-command))
      (format stream "camera = AutoCamera(image, direction=~S, background=~S, resolution=~D, aspect=1.0);~%"
	      (if (= dimension 3) "off-diagonal" "front")
	      (ecase background (:black "black") (:white "white"))
	      resolution)
      (if *dx-bug-workaround*
	  (format stream "image = Options(image, \"rendering mode\", \"hardware\");~%")
	  (format stream "image = Render(image,camera);"))
      (when (find :screen *output-types*)
        (let ((control "where=SuperviseWindow(~S,size=[~D,~D],visibility=~D);~%"))
          (when *show-dx-window*
            (format stream control window width height 2))
          (format stream control window width height 1))
        (if *dx-bug-workaround*
            ;; corresponding to the above problematic variant
            (format stream "Display (image, camera, where=where);~%")
            (format stream "Display (image, where=where);~%")))
      (awhen (find-if (_ (and (listp _) (eql (first _) :images)))
                      *output-types*)
        (whereas ((pathname (funcall (second it))))
          (flet ((file (type)
                   (namestring (merge-pathnames pathname (make-pathname :type type)))))
            (format stream "WriteImage(image, ~S, \"tiff\");" (file nil))
            ;; unfortunately writing jpg directly does not work (contrary to the DX Reference Manual):
            ;; "ImageMagick supported format gamma=2.2 compression=JPEG quality=90 resize=50"
            (handshake-with-dx)
            (awhen (probe-file "/usr/bin/convert")
              (fl.port:run-program it (list "-quality" "90" (file "tiff") (file "jpg"))
                                   :wait t))
            )))
      (when filename
	(ensure format "tiff")
	(cond
	  (*dx-bug-workaround*
	   (format stream "content = ReadImageWindow(where);~%")
	   (format stream "WriteImage (content,~S,~S);~%"
		   (concatenate 'string (namestring (images-pathname)) filename)
		   format))
	  (t
	   (format stream "WriteImage (image,~S,~S);~%"
		   (concatenate 'string (namestring (images-pathname)) filename)
		   format))))
      (handshake-with-dx)
      )))
