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

(in-package :graphics)

(defparameter *default-graphic-program* :dx
  "Default graphics program.")

(defparameter *images-pathname*
  (or (and cl-user::*images-directory*
	   (pathname cl-user::*images-directory*))
      (translate-logical-pathname #p"femlisp:images;"))
  "Pathname of the images directory.")

(defvar *images-directory*
  (namestring *images-pathname*)
  "Namestring of *IMAGES-PATHNAME*.")

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Public interface
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defgeneric graphic-stream (program)
  (:documentation "Return the stream for the program."))

(defgeneric graphic-file-name (object program &rest rest &key &allow-other-keys)
  (:documentation "Return a filename for the data of this plot."))

(defgeneric graphic-write-data (stream object program &rest rest &key &allow-other-keys)
  (:documentation "Will usually be the default method depending on several
data elements in the rest parameters."))

(defgeneric graphic-commands (object program &rest rest &key &allow-other-keys)
  (:documentation "Returns commands for plotting to be sent to the graphics
server."))

(defgeneric send-graphic-commands (stream object program &rest rest &key &allow-other-keys)
  (:documentation "Routine for sending commands to the graphics server."))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; General graphics output
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmethod graphic-output (object program &rest rest
			   &key debug &allow-other-keys)
  "Calls the generic graphic interface in appropriate order."
  (let* ((filename (apply #'graphic-file-name object program rest))
	 (pathname (concatenate 'string "femlisp:images;" filename)))
    ;; write output to a standard file
    (with-open-file (stream pathname :direction :output :if-exists :supersede)
      (apply #'graphic-write-data stream object program rest))
    ;; send script commands to plot program
    (whereas ((stream (graphic-stream program)))
      (when debug (apply #'send-graphic-commands *trace-output* object program rest))
      (apply #'send-graphic-commands stream object program rest)
      (force-output stream)))
  ;; this is quite useful, similar to the behaviour of print
  object)

