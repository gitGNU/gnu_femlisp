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

(pushnew :femlisp *features*)

;;; Note: The code below expects that the logical host "femlisp:" is
;;; correctly set to the femlisp directory.  Furthermore, "cl:matlisp"
;;; should point to the matlisp directory and cl:utilities to a directory
;;; containing infix.cl.

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Initialization of external utilities
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(load "femlisp:src;femlisp-config.lisp")

#-matlisp (setq debug:*debug-readtable* (copy-readtable *readtable*))
#-matlisp (load "cl:matlisp;start.lisp")
(unexport 'MATLISP::REAL :MATLISP)
#-infix (load "cl:utilities;infix.cl")

(load "cl:cl-ppcre;load.lisp")

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; register femlisp for defsystem
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; Register 'femlisp;src' for the defsystem utility
(let ((name "femlisp:src;"))
  (setq mk::*central-registry*
	(cons (translate-logical-pathname name)
	      (remove name mk::*central-registry* :test #'equal
		      :key #'(lambda (obj) (and (typep obj 'pathname)
						(pathname-name obj)))))))

;;(mk:oos 'femlisp 'compile :force :all :verbose t :compile-during-load t)
(mk:oos 'femlisp 'compile :minimal-load t :verbose t :compile-during-load t)

(eval-when (:load-toplevel :compile-toplevel :execute)
(defparameter *femlisp-version* "0.8.3")
(defun femlisp-version () *femlisp-version*)
(defun femlisp-herald () (format nil "    Femlisp/~a" (femlisp-version)))
#+cmu (setf (getf ext:*herald-items* :femlisp)
	    (list (femlisp-herald))))

(defun femlisp-banner ()
  (format
   t "~&~%*** Femlisp-0.8.5 ***

Copyright (C) 2003-2004
Nicolas Neuss, University of Heidelberg.

Femlisp comes with ABSOLUTELY NO WARRANTY, for details see the
file LICENSE in the Femlisp main directory.  This is free
software, and you are welcome to redistribute it under certain
conditions.

Type (demo) or (femlisp-demo:demo) to get a guided tour through
Femlisp.~%~%"))

(femlisp-banner)

