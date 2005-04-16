;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; port.lisp - portability issues between Lisp implementations
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

;;; This file serves the same purpose as the PORT module in CLOCC and is
;;; inspired by this module.  It will be dropped when CLOCC/port is easily
;;; installable in all CL implementations we are interested in or if the
;;; maintenance of this file should become too difficult.

(defpackage "FL.PORT"
  (:use "COMMON-LISP")
  (:import-from
   #+cmu "C-CALL" #+sbcl "SB-ALIEN" #+clisp "FFI"
   "INT" "DOUBLE")
  #+(or cmu sbcl)
  (:import-from
   #+cmu "SYSTEM" #+sbcl "SB-SYS"
   "WITHOUT-GCING")
  (:import-from
   #+(or cmu clisp) "EXT" #+sbcl "SB-EXT"
   "QUIT")
  (:export
   ;; UNIX environment
   "FIND-EXECUTABLE" "GETENV" "UNIX-CHDIR"
   
   ;; process communication
   "RUN-PROGRAM" "PROCESS-INPUT" "PROCESS-OUTPUT"
   "PROCESS-CLOSE" "PROCESS-STATUS"

   ;; load alien code
   "LOAD-SHARED-OBJECT" "DEFINE-ALIEN-ROUTINE"
   "INT" "DOUBLE" "LOAD-FOREIGN"
   "VECTOR-SAP" "WITHOUT-GCING")
  (:export "SAVE-FEMLISP-CORE-AND-DIE" "FEMLISP-RESTART")
  (:documentation "This package should contain the implementation-dependent
parts of Femlisp with the exception of the MOP."))

(in-package :fl.port)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; UNIX environment access
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun find-executable (name)
  #+cmu (probe-file (pathname (concatenate 'string "path:" name)))
  #+sbcl (sb-ext:find-executable-in-search-path name)
  #-(or cmu sbcl)
  (error "Unknown Lisp implementation.")
  )

(defun getenv (var)
  "Return the value of the environment variable."
  #+allegro (sys::getenv (string var))
  #+clisp (ext:getenv (string var))
  #+(or cmu scl)
  (cdr (assoc (string var) ext:*environment-list* :test #'equalp
              :key #'string))
  #+gcl (si:getenv (string var))
  #+lispworks (lw:environment-variable (string var))
  #+mcl (ccl::getenv var)
  #+sbcl (sb-ext:posix-getenv var)
  #-(or allegro clisp cmu gcl lispworks mcl sbcl scl)
  (error "Unknown Lisp implementation.")
  )
  
(defun unix-chdir (path)
  "Change the directory to @arg{path}."
  #+cmu (unix:unix-chdir path)
  #+sbcl (sb-posix:chdir path)
  #+clisp (ext:cd path)
  #-(or cmu sbcl clisp) nil
  )

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Process communication
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun run-program (program args &rest opts)
  "Interface to run-program."
  #+clisp (apply #'ext:run-program program :arguments args opts)
  #+cmu (apply #'ext:run-program program args opts)
  #+sbcl (apply #'sb-ext:run-program program args opts)
  #-(or clisp cmu sbcl) (error "Unknown Lisp implementation")
  )

(defun process-input (process)
  "Interface to process-input."
  #+cmu (ext:process-input process)
  #+sbcl (sb-ext:process-input process)
  #-(or cmu sbcl) (error "Unknown Lisp implementation")
  )

(defun process-output (process)
  "Interface to process-input."
  #+cmu (ext:process-output process)
  #+sbcl (sb-ext:process-output process)
  #-(or cmu sbcl) (error "Unknown Lisp implementation")
  )

(defun process-close (process)
  "Interface to process-input."
  #+cmu (ext:process-close process)
  #+sbcl (sb-ext:process-close process)
  #-(or cmu sbcl) (error "Unknown Lisp implementation")
  )

(defun process-status (process)
  "Interface to process-input."
  #+cmu (ext:process-status process)
  #+sbcl (sb-ext:process-status process)
  #-(or cmu sbcl) (error "Unknown Lisp implementation")
  )

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Foreign libraries
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun load-shared-object (file)
  #+cmu (sys::load-object-file file)
  #+sbcl (sb-alien:load-shared-object file))

(defmacro define-alien-routine (&rest args)
  #+cmu `(c-call::def-alien-routine ,@args)
  #+sbcl `(sb-alien:define-alien-routine ,@args)
  #+clisp
  (destructuring-bind (name return-type &rest args)
      args
    `(def-call-out (intern name)
	 (cons :arguments args)
       (list :return-type return-type))))

(declaim (inline vector-sap))
(defun vector-sap (ptr)
  #+cmu (system:vector-sap ptr)
  #+sbcl (sb-sys:vector-sap ptr)
  #+clisp (ffi::foreign-pointer ptr)
  )

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Saving a core and restarting
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun save-femlisp-core-and-die (core-file-name)
  #+cmu
  (progn
    (ext:save-lisp core-file-name :print-herald nil)
    (ext:quit))
  #+sbcl
  (sb-ext:save-lisp-and-die core-file-name)
  #-(or cmu sbcl) (error "Unknown Lisp implementation"))

(defun femlisp-restart ()
  #+cmu (progn (ext::print-herald))
  (fl.start:femlisp-banner))

  