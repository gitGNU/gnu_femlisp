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

(in-package :cl-user)

;;; Note: We expect that the logical host "femlisp:" is
;;; correctly set to the femlisp directory and that
;;; "femlisp:src;" is registered for ASDF.  CL-PPCRE and infix
;;; should also be loaded.  This can be ensured by loading
;;; "femlisp:femlisp-init.lisp" when starting up Lisp.

(load "femlisp:femlisp-config.lisp")

;;; we want to work generally with double float numbers
(setq *READ-DEFAULT-FLOAT-FORMAT* 'double-float)

#+asdf (asdf:operate 'asdf::load-op 'femlisp)
#-asdf (mk:oos 'femlisp 'compile)

(pushnew :femlisp *features*)

(defparameter *femlisp-version* "0.9.3")

(defun femlisp-version () *femlisp-version*)
(defun femlisp-herald () (format nil "    Femlisp/~a" (femlisp-version)))

(defun femlisp-banner ()
  (format
   t "~&~%*** Femlisp-~A ***

Copyright (C) 2003-2004
Nicolas Neuss, University of Heidelberg.

Femlisp comes with ABSOLUTELY NO WARRANTY, for details see the
file LICENSE in the Femlisp main directory.  This is free
software, and you are welcome to redistribute it under certain
conditions.

You can enter \"(femlisp-demo)\" to get a guided tour through
Femlisp, and enter \"(quit)\" to leave the program.~%~%"
*femlisp-version*))

(femlisp-banner)
