;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; graphics.lisp
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

(defparameter *default-graphic-program* :dx
  "Default graphics program.")

(defparameter *output-types*
  (list :screen)
  "Types of output.  Possible entries in this list are:
   - :screen
   - (:images successive-pathname-producer)")

(defun images-pathname ()
  "Pathname of the directory for @femlisp{} images."
  (or (aand (or (fl.port:getenv "FEMLISP_IMAGES")
		fl.start:*images-directory*)
	    (pathname it))
      (fl.start:femlisp-pathname "images/")))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Public interface
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defgeneric graphic-input-stream (program)
  (:documentation "Return the input stream for the graphic program."))
(defgeneric graphic-output-stream (program)
  (:documentation "Return the output stream for the graphic program."))

(defgeneric graphic-file-name (object program &rest rest &key &allow-other-keys)
  (:documentation "Return a filename for the data of this plot.")
  (:method :around (object program &key system-p &allow-other-keys)
  "This around-method handles the system-p flag: if true,
@function{graphic-file-name} returns a filename suitable for use in
external programs."
  (declare (ignore object program))
  (let ((result (call-next-method)))
    (if system-p
        (fl.port::system-namestring (probe-file result))
        result))))

(defgeneric graphic-write-data (stream object program &rest rest &key &allow-other-keys)
  (:documentation "Write the data file for @arg{program} to
@arg{stream}."))

(defgeneric graphic-commands (object program &rest rest &key &allow-other-keys)
  (:documentation "Returns commands for plotting to be sent to the graphics
program."))

(defgeneric send-graphic-commands (object program &rest rest &key &allow-other-keys)
  (:documentation "Routine for sending commands to the graphics server."))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; General graphics output
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defgeneric graphic-output (object program &key &allow-other-keys)
  (:documentation "Writes the graphic commands for @arg{object}
as needed by the graphic program @arg{program}.")
  (:method (object program &rest rest)
      "Calls the generic graphic interface in appropriate order."
    (let* ((pathname (apply #'graphic-file-name object program rest)))
      ;; write output to a standard file
      (with-open-file (stream pathname :direction :output :if-exists :supersede)
        (apply #'graphic-write-data stream object program rest))
      ;; wait until the file is there
      (loop until (probe-file pathname) do (sleep 0.01))
      ;; send script commands to plot program
      (apply #'send-graphic-commands object program rest))))

