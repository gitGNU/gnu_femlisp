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

(defpackage "FL.PORT"
  (:use "COMMON-LISP" "FL.UTILITIES")
  #+(or cmu sbcl)
  (:import-from
   #+cmu "C-CALL" #+sbcl "SB-ALIEN"
   "INT" "DOUBLE")
  #+(or cmu sbcl)
  (:import-from
   #+cmu "SYSTEM" #+sbcl "SB-SYS"
   "WITHOUT-GCING")
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

(file-documentation
 "This file serves the same purpose as the PORT module in CLOCC and is
inspired by this module.  It will be dropped when CLOCC/port is easily
installable in all CL implementations we are interested in or if the
maintenance of this file should become too difficult.")

;;;; required modules
#+allegro (require :osi)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Utility
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun portability-warning (function &optional (error t))
  (let ((message (format nil "The function ~A has not yet been written for
your Lisp.  If you want full functionality of Femlisp you should provide it
in the file @path{femlisp:src;basic;port.lisp}." function)))
    (if error
	(error message)
	(warn message))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Quitting
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun quit ()
  #+cmu (ext:quit)
  #+sbcl (sb-ext:quit)
  #+allegro (excl:exit)
  #-(or allegro cmu sbcl)
  (portability-warning 'quit)
  )

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; UNIX environment access
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun find-executable (name)
  #-(or allegro cmu sbcl)
  (portability-warning 'find-executable nil)
  #+allegro (excl.osi:find-in-path name)
  #+cmu (probe-file (pathname (concatenate 'string "path:" name)))
  #+sbcl (sb-ext:find-executable-in-search-path name)
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
  (portability-warning 'getenv)
  )
  
(defun unix-chdir (path)
  "Change the directory to @arg{path}."
  #+allegro (excl.osi:chdir path)
  #+cmu (unix:unix-chdir path)
  #+sbcl (sb-posix:chdir path)
  #+clisp (ext:cd path)
  #-(or cmu sbcl clisp)
  (portability-warning 'unix-chdir)
  )

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Process communication
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun run-program (program args &key wait input output error-output directory)
  "Interface to run-program."
  #+allegro
  (let ((program (namestring program)))
    (cond
      (wait (excl:run-shell-command (apply #'vector program program args) :wait t
				    :directory directory))
      (t
       ;; for now, process-close expects that only input is a newly generated stream
       (assert (and (eq input :stream)
		    (not (eq output :stream))
		    (not (eq error-output :stream))))
       (multiple-value-bind (istream ostream estream proc-id)
	   (excl:run-shell-command
	    (apply #'vector program program args)
	    :wait nil :separate-streams t :input input :output output
	    :error-output error-output :directory directory)
	 (declare (ignore ostream estream))
	 (list :process istream output error-output proc-id)))))
  #+clisp (ext:run-program program :arguments args :wait wait :input input :output output)
  #+cmu (progn
	  (when directory (unix-chdir (namestring directory)))
	  (ext:run-program program args :wait wait :input input :output output))
  #+sbcl (progn
	   (when directory (unix-chdir (namestring directory)))
	   (sb-ext:run-program program args :wait wait :input input :output output))
  #-(or allegro clisp cmu sbcl)
  (portability-warning 'run-program)
  )

(defun process-input (process)
  "Interface to process-input."
  #+allegro (second process)
  #+cmu (ext:process-input process)
  #+sbcl (sb-ext:process-input process)
  #-(or allegro cmu sbcl)
  (portability-warning 'process-input)
  )

(defun process-output (process)
  "Interface to process-output."
  #+allegro (third process)
  #+cmu (ext:process-output process)
  #+sbcl (sb-ext:process-output process)
  #-(or allegro cmu sbcl)
  (portability-warning 'process-output)
  )

(defun process-close (process)
  "Interface to process-close."
  #+allegro (close (second process))
  #+cmu (ext:process-close process)
  #+sbcl (sb-ext:process-close process)
  #-(or allegro cmu sbcl)
  (portability-warning 'process-close)
  )

(defun process-status (process)
  "Interface to process-status."
  #+cmu (ext:process-status process)
  #+sbcl (sb-ext:process-status process)
  #-(or cmu sbcl)
  (portability-warning 'process-status)
  )

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Foreign libraries
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun load-shared-object (file)
  #+cmu (sys::load-object-file file)
  #+sbcl (sb-alien:load-shared-object file)
  #-(or cmu sbcl)
  (portability-warning 'load-shared-object))

(defmacro define-alien-routine (&rest args)
  #+cmu `(c-call::def-alien-routine ,@args)
  #+sbcl `(sb-alien:define-alien-routine ,@args)
  #+clisp
  (destructuring-bind (name return-type &rest args)
      args
    `(def-call-out (intern name)
	 (cons :arguments args)
       (list :return-type return-type)))
  #-(or cmu sbcl clisp)
  (portability-warning 'define-alien-routine))

(declaim (inline vector-sap))
(defun vector-sap (ptr)
  #+cmu (system:vector-sap ptr)
  #+sbcl (sb-sys:vector-sap ptr)
  #+clisp (ffi::foreign-pointer ptr)
  #-(or cmu sbcl clisp)
  (portability-warning 'vector-sap))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Saving a core and restarting
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun save-femlisp-core-and-die (core-file-name)
  #+allegro
  (progn
    (excl:dumplisp :name core-file-name)
    (exit 0))
  #+cmu
  (progn
    (ext:save-lisp core-file-name :print-herald nil)
    (ext:quit))
  #+sbcl
  (sb-ext:save-lisp-and-die core-file-name)
  #+clisp
  (EXT:SAVEINITMEM corefilename)
  #-(or cmu sbcl clisp)  ; do nothing
  (portability-warning 'save-femlisp-core-and-die))

(defun femlisp-restart ()
  #+cmu (progn (ext::print-herald))
  (fl.start:femlisp-banner))

  