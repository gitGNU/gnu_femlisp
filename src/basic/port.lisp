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
  (:use "COMMON-LISP")
  (:export
   "PORTABILITY-WARNING"
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
parts of Femlisp with the exception of the MOP.  It serves a similar
purpose as the PORT module in CLOCC and is somewhat inspired by this
module.  It will be dropped when there is a portable and easily installable
alternative in all CL implementations we are interested in"))

(in-package :fl.port)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Utility
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defparameter *portability-problem-handling* :warn
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
;;;; UNIX environment access
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun find-executable (name)
  "Finds an executable in the current path."
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
  "Runs @arg{program} with arguments @arg{args}."
  #+(or clisp ecl cmu sbcl) (declare (ignore error-output))
  #+(or ecl gcl) (declare (ignore wait))
  #+allegro
  (let ((program (namestring program)))
    (cond
      (wait (excl:run-shell-command (apply #'vector program program args) :wait t
				    :directory directory :show-window :hide))
      (t
       #+mswindows (assert (and (eq input :stream) (eq output :stream)))
       (multiple-value-bind (stream null proc-id)
           (excl:run-shell-command
            #-mswindows (apply #'vector program program args)
            #+mswindows (apply #'concatenate 'string program
                               (loop for arg in args collect " " collect arg))
            :wait nil :input input :output output
            :error-output error-output :directory directory
            :show-window :normal) ; :hide
	 (declare (ignore null))
	 (list :allegro-process
	       (if (eq input :stream) stream input)
	       (if (eq output :stream) stream output)
	       error-output proc-id)))))
  #+clisp (ext:run-program program :arguments args :wait wait :input input :output output)
  #+(or cmu ecl gcl sbcl)
  (when directory (unix-chdir (namestring directory)))
  #+cmu (ext:run-program program args :wait wait :input input :output output)
  #+(or gcl) (si:run-process (namestring program) args)
  #+ecl (si:run-program (namestring program) args)
  #+sbcl (sb-ext:run-program program args :wait wait :input input :output output)
  #-(or allegro clisp cmu ecl gcl sbcl)
  (portability-warning 'run-program args wait input output error-output directory)
  )

(defun process-input (process)
  "Process-input for @arg{process}."
  #+allegro (second process)
  #+cmu (ext:process-input process)
  #+ecl process
  #+gcl (si::fp-output-stream process)
  #+sbcl (sb-ext:process-input process)
  #-(or allegro cmu ecl gcl sbcl)
  (portability-warning 'process-input process)
  )

(defun process-output (process)
  "Process-output for @arg{process}."
  #+allegro (third process)
  #+cmu (ext:process-output process)
  #+ecl process
  #+gcl (si::fp-input-stream process)
  #+sbcl (sb-ext:process-output process)
  #-(or allegro cmu ecl gcl sbcl)
  (portability-warning 'process-output process)
  )

(defun process-close (process)
  "Closes @arg{process}."
  #+allegro (close (second process))
  #+cmu (ext:process-close process)
  #+(or ecl gcl) (close process)
  #+sbcl (sb-ext:process-close process)
  #-(or allegro cmu ecl gcl sbcl)
  (portability-warning 'process-close process)
  )

(defun process-status (process)
  "Returns the status of @arg{process}."
  #+cmu (ext:process-status process)
  #+sbcl (sb-ext:process-status process)
  #-(or cmu sbcl)
  (portability-warning 'process-status process)
  )

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Memory
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun memory-usage ()
  "Returns an approximation to the memory used at this moment."
  #+cmu (lisp::dynamic-usage)
  #+sbcl (sb-kernel::dynamic-usage)
  #-(or cmu sbcl)
  (portability-warning 'memory-usage))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Foreign libraries
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun load-foreign-library (file)
  "Loads the foreign library @arg{file}."
  #+allegro (cl:load file)
  #+cmu (sys::load-object-file file)
  #+ecl (ffi:load-foreign-library
	 (if (stringp file) file (namestring file)))
  #+sbcl (sb-alien:load-shared-object file)
  #+(and uffi (not (or allegro cmu ecl sbcl)))
  (uffi:load-foreign-library file)
  #-(or allegro cmu ecl sbcl uffi)
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
  "Defines a foreign function.  See examples in
@path{alien;src;superlu.lisp}."
  #+(or allegro cmu sbcl) `(simplified-def-function ,@args)
  #+ecl `(ffi:def-function ,@args)
  #+(and uffi (not (or allegro cmu ecl sbcl)))
  `(uffi:def-function ,@args)
  #-(or allegro cmu ecl sbcl uffi)
  (portability-warning 'define-alien-routine args))

(declaim (inline vector-sap))
(defun vector-sap (ptr)
  "Returns an array pointer which can be used in a foreign call."
  #+allegro ptr
  #+cmu (system:vector-sap ptr)
  #+ecl (ffi:c-inline (ptr) (:object) :pointer-void "(#0)->array.self.ch" :one-liner t)
  #+sbcl (sb-sys:vector-sap ptr)
  #+clisp (ffi::foreign-pointer ptr)
  #-(or allegro cmu ecl sbcl clisp)
  (portability-warning 'vector-sap ptr)
  )

(defmacro foreign-call-wrapper (&rest body)
  "Ensures a safe environment for a foreign function call, especially so
that no GC changes array pointers obtained by @function{vector-sap}."
  #+(or allegro ecl) `(progn ,@body)
  #+cmu `(system:without-gcing ,@body)
  #+sbcl `(sb-sys:without-gcing ,@body)
  #-(or allegro cmu ecl sbcl)
  (portability-warning 'foreign-call-wrapper ptr)
  )

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Weak pointers
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun make-weak-pointer (obj)
  "Creates a weak pointer pointing to @arg{obj}."
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
  "Returns the value of the weak pointer @arg{wp}."
  #+allegro (aref wp 0)
  #+cmu (ext:weak-pointer-value wp)
  #+sbcl (sb-ext:weak-pointer-value wp)
  #-(or allegro cmu sbcl)
  (portability-warning 'weak-pointer-value wp))

(defun finalize (obj func)
  "Sets up @arg{func} as finalizer for @arg{obj}."
  #+allegro (excl:schedule-finalization
	     obj #'(lambda (obj)
		     (declare (ignore obj))
		     (funcall func)))
  #+cmu (ext:finalize obj func)
  #+sbcl (sb-ext:finalize obj func)
  #-(or allegro cmu sbcl)
  (portability-warning 'finalize obj func))

(defun gc (&rest args)
  #+allegro (apply #'excl:gc args)
  #+cmu (apply #'ext:gc args)
  #+sbcl (apply #'sb-ext:gc args)
  #-(or allegro cmu sbcl)
  (portability-warning 'gc args))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Quitting
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun quit ()
  "Quits Femlisp."
  #+cmu (ext:quit)
  #+ecl (si:quit)
  #+gcl (lisp:quit)
  #+sbcl (sb-ext:quit)
  #+allegro (excl:exit)
  #-(or allegro cmu sbcl)
  (portability-warning 'quit)
  )

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Saving a core and restarting
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun save-femlisp-core-and-die (core-file-name)
  "Saves Femlisp core and quits."
  (when (symbolp core-file-name)
    (setq core-file-name (string-downcase (symbol-name core-file-name))))
  #+allegro
  (progn
    (setq excl:*restart-init-function*
	  (let ((directory (pathname-directory fl.start::*femlisp-pathname*)))
	    #'(lambda ()
		(setf (logical-pathname-translations "FEMLISP")
		      `(("**;*.*.*"
			 ,(make-pathname :directory `(,@directory :wild-inferiors)
					 :name :wild :type :wild :version :wild))))
		(tpl:setq-default *package* (find-package :fl.application))
		(rplacd (assoc 'tpl::*saved-package*
			       tpl:*default-lisp-listener-bindings*)
			'common-lisp:*package*))))
    (excl:dumplisp :name core-file-name))
  #+clisp   (EXT:SAVEINITMEM corefilename)
  #+cmu     (ext:save-lisp core-file-name :print-herald nil)
  #+gcl     (si:save-system core-file-name)
  #+ecl     (quit)
  #+sbcl    (sb-ext:save-lisp-and-die core-file-name)
  #-(or allegro clisp cmu ecl gcl sbcl)  ; do nothing
  (portability-warning 'save-femlisp-core-and-die core-file-name)
  ;; we quit in any case
  (quit))

(defun femlisp-restart ()
  #+cmu (progn (ext::print-herald))
  (fl.start::femlisp-banner)
  (in-package :fl.application)
  )
  
