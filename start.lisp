;;; -*- mode: lisp; fill-column: 64; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; start.lisp - Femlisp start file
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

(in-package "COMMON-LISP-USER")

(defpackage "FL.START"
  (:use "COMMON-LISP")
  (:export "*FEMLISP-VERSION*" "FEMLISP-HERALD" "FEMLISP-BANNER")
  (:documentation "This package contains some routines called
during initialization of Femlisp."))

(in-package :fl.start)

(defparameter *femlisp-version* "0.9.6")
(defparameter *process* nil
  "This variable should be set externally for identifying a certain process.")

(defun femlisp-version () *femlisp-version*)
(defun femlisp-herald () (format nil "    Femlisp/~a" (femlisp-version)))

(defun femlisp-banner ()
  (format
   t "~&~%*** Femlisp-~A ***

Copyright (C) 2003-2005
Nicolas Neuss, University of Heidelberg.

Femlisp comes with ABSOLUTELY NO WARRANTY, for details see the
file LICENSE in the Femlisp main directory.  This is free
software, and you are welcome to redistribute it under certain
conditions.

You can enter \"(femlisp-demo)\" to get a guided tour through
Femlisp, and enter \"(quit)\" to leave the program.~%~%"
*femlisp-version*))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Setup the logical host "FEMLISP"
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defvar *femlisp-pathname*
  (make-pathname :directory (pathname-directory *load-truename*))
  "The pathname for the Femlisp main directory.  This should be the
location of this file when it is loaded.")

(defvar *femlisp-directory* (namestring *femlisp-pathname*)
  "The namestring for @var{*femlisp-pathname*}.")

(let ((directory (pathname-directory *femlisp-pathname*)))
  (setf (logical-pathname-translations "FEMLISP")
	`(("**;*.*.*" 
	   ,(make-pathname :directory `(,@directory :wild-inferiors)
			   :name :wild :type :wild :version :wild)))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Configuration
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(load "femlisp:femlisp-config.lisp")

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Implementation-dependent
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

#+cmu (setq extensions:*gc-verbose* nil)
#+sbcl (require 'sb-posix)
#+sbcl (require 'sb-introspect)
#+sbcl (require 'asdf)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Ensure the presence of Common Lisp libraries
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; ASDF
#-(or asdf mk-defsystem) (load #p"femlisp:external;asdf")

;;; INFIX
#-infix (load "femlisp:external;infix.cl")

#+asdf 
(pushnew *femlisp-pathname* asdf::*central-registry*)
#+mk-defsystem
(pushnew *femlisp-pathname* mk::*central-registry*)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Compiling and loading of Femlisp
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; we want to work generally with double float numbers
(setq *READ-DEFAULT-FLOAT-FORMAT* 'double-float)

(load #p"femlisp:femlisp.asd")

#+asdf (asdf:operate 'asdf::load-op 'femlisp)
#-asdf (mk:oos 'femlisp 'compile)

(pushnew :femlisp *features*)

(let ((private (probe-file #p"femlisp:start-private.lisp")))
  (when private (load private)))

(femlisp-banner)

