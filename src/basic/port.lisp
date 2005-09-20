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
  (:export
   ;; UNIX environment
   "FIND-EXECUTABLE" "GETENV" "UNIX-CHDIR"
   
   ;; process communication
   "RUN-PROGRAM" "PROCESS-INPUT" "PROCESS-OUTPUT"
   "PROCESS-CLOSE" "PROCESS-STATUS"

   ;; load alien code
   "LOAD-FOREIGN-LIBRARY" "DEF-FUNCTION"
   "INT" "DOUBLE" "LOAD-FOREIGN"
   "VECTOR-SAP" "FOREIGN-CALL-WRAPPER"

   ;; weak pointers and GC
   "MAKE-WEAK-POINTER" "WEAK-POINTER-VALUE"
   "FINALIZE" "GC"
   )
  (:export "SAVE-FEMLISP-CORE-AND-DIE" "FEMLISP-RESTART")
  (:documentation "This package should contain the implementation-dependent
parts of Femlisp with the exception of the MOP."))

(in-package :fl.port)

(file-documentation
 "This file serves the same purpose as the PORT module in CLOCC and is
inspired by this module.  It will be dropped when CLOCC/port is easily
installable in all CL implementations we are interested in or if the
maintenance of this file should become too difficult.")

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Utility
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defparameter *portability-problem-handling* :error
  "Determines which action should be taken if we encounter a portability
problem for a Lisp.")

(defun portability-warning (function &rest args)
  (let ((message (format nil "The function~& ~A~% called with
parameters~&~A~%has not yet been written for your Lisp.  If you want full
functionality of Femlisp you should provide it in the file
@path{femlisp:src;basic;port.lisp}." function args)))
    (ecase *portability-problem-handling*
      (:error (error message))
      (:warn (warn message)))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Quitting
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun quit ()
  #+cmu (ext:quit)
  #+ecl (si:quit)
  #+gcl (lisp:quit)
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
  (portability-warning 'find-executable name)
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
  #+ecl (si:getenv (string var))
  #+gcl (system:getenv (string var))
  #+lispworks (lw:environment-variable (string var))
  #+mcl (ccl::getenv var)
  #+sbcl (sb-ext:posix-getenv var)
  #-(or allegro clisp cmu ecl gcl lispworks mcl sbcl scl)
  (portability-warning 'getenv var)
  )
  
(defun unix-chdir (path)
  "Change the directory to @arg{path}."
  #+allegro (excl.osi:chdir path)
  #+cmu (unix:unix-chdir path)
  #+ecl (si:chdir path)
  #+gcl (system:chdir path)
  #+sbcl (sb-posix:chdir path)
  #+clisp (ext:cd path)
  #-(or cmu clisp ecl gcl sbcl)
  (portability-warning 'unix-chdir path)
  )

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Process communication
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun run-program (program args &key wait input output error-output directory)
  "Interface to run-program."
  #+(or clisp cmu sbcl) (declare (ignore error-output))
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
  ;; #+(or ecl gcl) (apply #'si:run-process prog args)
  #+sbcl (progn
	   (when directory (unix-chdir (namestring directory)))
	   (sb-ext:run-program program args :wait wait :input input :output output))
  #-(or allegro clisp cmu sbcl)
  (portability-warning 'run-program args wait input output error-output directory)
  )

(defun process-input (process)
  "Interface to process-input."
  #+allegro (second process)
  #+cmu (ext:process-input process)
  #+sbcl (sb-ext:process-input process)
  #-(or allegro cmu sbcl)
  (portability-warning 'process-input process)
  )

(defun process-output (process)
  "Interface to process-output."
  #+allegro (third process)
  #+cmu (ext:process-output process)
  #+sbcl (sb-ext:process-output process)
  #-(or allegro cmu sbcl)
  (portability-warning 'process-output process)
  )

(defun process-close (process)
  "Interface to process-close."
  #+allegro (close (second process))
  #+cmu (ext:process-close process)
  #+sbcl (sb-ext:process-close process)
  #-(or allegro cmu sbcl)
  (portability-warning 'process-close process)
  )

(defun process-status (process)
  "Interface to process-status."
  #+cmu (ext:process-status process)
  #+sbcl (sb-ext:process-status process)
  #-(or cmu sbcl)
  (portability-warning 'process-status process)
  )

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Memory
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun memory-usage ()
  "Returns an approixmation to the memory used at this moment."
  #+cmu (lisp::dynamic-usage)
  #+sbcl (sb-kernel::dynamic-usage)
  #-(or cmu sbcl)
  (portability-warning 'memory-usage))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Foreign libraries
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun load-foreign-library (file)
  #+allegro (cl:load file)
  #+cmu (sys::load-object-file file)
  #+sbcl (sb-alien:load-shared-object file)
  #+(and uffi (not (or allegro cmu sbcl)))
  (uffi:load-foreign-library file)
  #-(or allegro cmu sbcl uffi)
  (portability-warning 'load-foreign-library file))

#+(or cmu sbcl)
(defun convert-type (type)
  (cond ((null type) ())
	((atom type)
	 (if (member type '(:int :double))
	     (find-symbol (symbol-name type)
			  #+cmu "C-CALL" #+sbcl "SB-ALIEN")
	     type))
	(t (cons (convert-type (car type)) (convert-type (cdr type))))))

#+(or allegro cmu sbcl)
(defmacro simplified-def-function ((c-name lisp-name) args &rest keys)
  (declare (ignorable lisp-name))
  `(#+cmu c-call::def-alien-routine
    #+sbcl sb-alien::define-alien-routine
    #+allegro ff:def-foreign-call
    #+allegro (,lisp-name ,c-name)
    #+(or cmu sbcl) ,c-name
    #+(or cmu sbcl) ,(convert-type (getf keys :returning))
    #+allegro ,args
    #+(or cmu sbcl)
    ,@(mapcar #'(lambda (pair) (cons (car pair) (convert-type (cdr pair)))) args)
    #+allegro ,@keys))

(defmacro def-function (&rest args)
  #+(or allegro cmu sbcl) `(simplified-def-function ,@args)
  #+(and uffi (not (or allegro cmu sbcl)))
  `(uffi:def-function ,@args)
  #-(or allegro cmu sbcl uffi)
  (portability-warning 'define-alien-routine args))

(declaim (inline vector-sap))
(defun vector-sap (ptr)
  #+allegro ptr
  #+cmu (system:vector-sap ptr)
  #+sbcl (sb-sys:vector-sap ptr)
  #+clisp (ffi::foreign-pointer ptr)
  #-(or allegro cmu sbcl clisp)
  (portability-warning 'vector-sap ptr)
  )

(defmacro foreign-call-wrapper (&rest body)
  "Ensures a safe environment for a FF call, especially so that no GC
changes array pointers obtained by @function{vector-sap}."
  #+allegro `(progn ,@body)
  #+cmu `(system:without-gcing ,@body)
  #+sbcl `(sb-sys:without-gcing ,@body)
  #-(or allegro cmu sbcl)
  (portability-warning 'foreign-call-wrapper ptr)
  )

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Weak pointers
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun make-weak-pointer (obj)
  #+allegro
  (let ((result (excl:weak-vector 1)))
    (setf (aref result 0) obj)
    result)
  #+cmu (ext:make-weak-pointer obj)
  #+sbcl (sb-ext:make-weak-pointer obj)
  #-(or allegro cmu sbcl)
  (portability-warning 'make-weak-pointer obj)
  )

(defun weak-pointer-value (wp)
  #+allegro (aref wp 0)
  #+cmu (ext:weak-pointer-value wp)
  #+sbcl (sb-ext:weak-pointer-value wp)
  #-(or allegro cmu sbcl)
  (portability-warning 'weak-pointer-value wp))

(defun finalize (obj func)
  #+allegro (excl:schedule-finalization
	     obj #'(lambda (obj)
		     (declare (ignore obj))
		     (funcall func)))
  #+cmu (ext:finalize obj func)
  #+sbcl (sb-ext: finalize obj func)
  #-(or allegro cmu sbcl)
  (portability-warning 'finalize obj func))

(defun gc (&rest args)
  #+allegro (apply #'excl:gc args)
  #+cmu (apply #'ext:gc args)
  #+sbcl (apply #'sb-ext:gc args)
  #-(or allegro cmu sbcl)
  (portability-warning 'gc args))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Saving a core and restarting
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun save-femlisp-core-and-die (core-file-name)
  #+allegro
  (progn
    (excl:dumplisp :name core-file-name)
    (excl:exit 0))
  #+cmu
  (progn
    (ext:save-lisp core-file-name :print-herald nil)
    (ext:quit))
  #+sbcl
  (sb-ext:save-lisp-and-die core-file-name)
  #+clisp
  (EXT:SAVEINITMEM corefilename)
  #-(or cmu sbcl clisp)  ; do nothing
  (portability-warning 'save-femlisp-core-and-die core-file-name))

(defun femlisp-restart ()
  #+cmu (progn (ext::print-herald))
  (fl.start::femlisp-banner))

  