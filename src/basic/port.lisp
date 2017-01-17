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
  (:use "COMMON-LISP" "FL.DEBUG")
  (:export
   "PORTABILITY-WARNING"
   ;; UNIX environment
   "FIND-EXECUTABLE" "FIND-SHARED-LIBRARY"
   "HOSTNAME" "GETENV" "UNIX-CHDIR" "SYSTEM-NAMESTRING"
   "DYNAMIC-SPACE-SIZE"
   
   ;; runtime compilation
   "COMPILE-SILENTLY" "RUNTIME-COMPILE" "COMPILE-AND-EVAL"
   
   ;; process communication
   "RUN-PROGRAM" "PROCESS-INPUT" "PROCESS-OUTPUT" "PROCESS-ERROR"
   "PROCESS-CLOSE" "PROCESS-STATUS" "RUN-PROGRAM-OUTPUT"
   "KILL-PROCESS"
   
   ;; load alien code
   "LOAD-FOREIGN-LIBRARY" "DEF-FUNCTION" "SIMPLIFIED-DEF-FUNCTION"
   "CONVERT-TYPE"
   "VECTOR-SAP" "FOREIGN-CALL"

   ;; weak pointers and GC
   "MAKE-WEAK-POINTER" "WEAK-POINTER-VALUE"
   "FINALIZE" "GC"

   ;; exiting
   "ADD-EXIT-HOOK"
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
@path{src/basic/port.lisp}." function args)))
    (ecase *portability-problem-handling*
      (:error (error message))
      (:warn (warn message)))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Runtime compilation with warning suppression
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun compile-silently (name source)
  "Compiles @arg{source} silently."
  ;; ANSI version
  (let ((*compile-print* nil)
	(*compile-verbose* nil))
    (compile name source))
  ;; posted by C. Rhodes 5.6.07 in the SBCL mailing list; however, it does
  ;; not work here
  #+(or)
  (handler-bind ((warning #'muffle-warning)
		 #+sbcl (sb-ext:compiler-note #'muffle-warning))
    (compile name source))
  )

(defun runtime-compile (source)
  "Calls compile on the provided @arg{source}.  When :compile is activated
for debugging, the source code is printed."
  (let ((*print-circle* nil))
    (dbg :compile "Compiling source: ~%~S~%" source))
  (funcall (if (dbg-p :compile) #'compile #'fl.port:compile-silently)
	   nil source))

(defun compile-and-eval (source)
  "Compiles and evaluates the given @arg{source}.  This should be an ANSI
compatible way of ensuring method compilation."
  (dbg :compile "Compiling and evaluating: ~%~S~%" source)
  (funcall (funcall (if (dbg-p :compile) #'compile #'fl.port:compile-silently)
		    nil `(lambda () ,source))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; CL environment
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun dynamic-space-size ()
  "Available memory for calculations"
  #+sbcl (sb-ext:dynamic-space-size)
  #-sbcl nil
  )

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; UNIX environment access
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
(defun hostname ()
  "Returns the hostname."
  #+asdf3 (uiop:hostname)
  #-asdf3 (getenv "HOSTNAME"))

(defun find-executable (name)
  "Finds an executable in the current path."
  #-(or allegro cmu scl)
  (declare (ignorable name))
  #+os-unix (probe-file (concatenate 'string "/usr/bin/" name))
  #+allegro (excl.osi:find-in-path name)
  #+cmu (probe-file (pathname (concatenate 'string "path:" name)))
  #+scl (probe-file (concatenate 'string "file://path/" name))
  )

(defun find-shared-library (name)
  "Finds a shared library."
  #-(or allegro cmu scl)
  (declare (ignorable name))
  #+os-unix (probe-file (concatenate 'string "/usr/lib/" name ".so"))
  #+allegro (excl.osi:find-in-path name)
  #+cmu (probe-file (pathname (concatenate 'string "path:" name)))
  #+scl (probe-file (concatenate 'string "file://path/" name))
  )

(defun getenv (var)
  "Return the value of the environment variable."
  (asdf::getenv var))

;; TODO: use uiop:chdir
(defun unix-chdir (path)
  "Change the directory to @arg{path}."
  #+allegro (excl.osi:chdir path)
  #+ccl (ccl::%chdir path)
  #+clisp (ext:cd path)
  #+(or cmu scl) (unix:unix-chdir path)
  #+ecl (si:chdir path)
  #+gcl (system:chdir path)
  #+lispworks (harlequin-common-lisp:change-directory path)
  #+sbcl (sb-posix:chdir path)
  #-(or allegro ccl clisp cmu scl ecl gcl lispworks sbcl)
  (portability-warning 'unix-chdir path))

;; TODO: use uiop:native-namestring
(defun system-namestring (path)
  #+scl (ext:unix-namestring path)
  #-scl (namestring path)
  )

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Process communication
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;; TODO: use uiop:run-program for synchronous use, and
;; uiop/run-program::%run-program for the asynchronous case
(defun run-program (program args
                    &key wait directory
                    input (output nil output-p) error-output)
  "Runs @arg{program} with arguments @arg{args}."
  (declare (ignorable output-p))
  #+(or clisp ecl) (declare (ignore error-output))
  #+(or ecl gcl) (declare (ignore wait))
  ;; change directory if not possible in the call itself
  #+(or cmu scl clisp ecl gcl lispworks sbcl ccl)
  (when directory (unix-chdir (namestring directory)))
  #+(or allegro lispworks)
  (let ((program (namestring program)))
    (cond
      (wait #+allegro (excl:run-shell-command
		       (apply #'vector program program args) :wait t
		       :directory directory :show-window :hide)
	    #+lispworks (system:run-shell-command
			 (apply #'vector program program args) :wait t))
      (t
       #+mswindows (assert (and (eq input :stream) (eq output :stream)))
       (multiple-value-bind (istream ostream error-output proc-id)
           #+allegro (excl:run-shell-command
		      #-mswindows (apply #'vector program program args)
		      #+mswindows (apply #'concatenate 'string program
					 (loop for arg in args collect " " collect arg))
		      :wait nil :input input :output output
		      :separate-streams t
		      :error-output error-output :directory directory
		      :show-window :normal) ; :hide
	   #+lispworks (system:run-shell-command
			(apply #'vector program program args)
			:wait nil :input input :output output
			:error-output error-output)
	 (declare (ignore null))
	 (list :process
	       (if (eq input :stream) istream input)
	       (if (eq output :stream) ostream output)
	       error-output proc-id)))))
  #+clisp
  (let ((output (if output-p output :terminal)))
    (cond
      (wait (ext:run-program program :arguments args :wait wait :input input :output output))
      (t (multiple-value-bind (iostream ostream istream)
             (ext:run-program program :arguments args :wait wait :input input :output output)
           (list :process
                 (or istream input)
                 (or ostream output)
                 iostream
                 nil)))))
  #+ccl
  (ccl:run-program program args :wait wait
                   :input input :output output :error error-output)
  #+(or cmu scl)
  (ext:run-program program args :wait wait
                   :input input :output output :error error-output)
  #+(or gcl) (si:run-process (namestring program) args)
  #+ecl (si:run-program (namestring program) args)
  #+sbcl (sb-ext:run-program program args :wait wait
			     :input input :output output :error error-output)
  #-(or allegro ccl clisp cmu ecl gcl lispworks sbcl scl)
  (portability-warning 'run-program args wait input output error-output directory)
  )

(defun process-input (process)
  "Process-input for @arg{process}."
  #+(or allegro lispworks clisp) (second process)
  #+(or cmu scl) (ext:process-input process)
  #+ccl (ccl:external-process-input-stream process)
  #+ecl process
  #+gcl (si::fp-output-stream process)
  #+sbcl (sb-ext:process-input process)
  #-(or allegro ccl lispworks clisp cmu scl ecl gcl sbcl)
  (portability-warning 'process-input process)
  )

(defun process-output (process)
  "Process-output for @arg{process}."
  #+(or allegro lispworks clisp) (third process)
  #+(or cmu scl) (ext:process-output process)
  #+ccl (ccl:external-process-output-stream process)
  #+ecl process
  #+gcl (si::fp-output-stream process)
  #+sbcl (sb-ext:process-output process)
  #-(or allegro ccl lispworks clisp cmu scl ecl gcl sbcl)
  (portability-warning 'process-output process)
  )

(defun process-error (process)
  "Process-output for @arg{process}."
  #+(or allegro lispworks) (fourth process)
  #+ccl (ccl:external-process-error-stream process)
  #+(or cmu scl) (ext:process-error process)
  #+sbcl (sb-ext:process-error process)
  #-(or allegro ccl lispworks cmu scl sbcl)
  (portability-warning 'process-error process)
  )

(defun process-close (process)
  "Closes @arg{process}."
  #+(or allegro lispworks) (close (second process))
  #+(or clisp)
  (loop for k in '(1 2 3) do
       (let ((stream (nth k process)))
         (when stream (close stream))))
  #+ccl (ccl:signal-external-process process 9)
  #+cmu (ext:process-close process)
  #+scl (ext:process-kill process :sigkill)
  #+(or ecl gcl) (close process)
  #+sbcl (sb-ext:process-close process)
  #-(or allegro ccl lispworks clisp cmu scl ecl gcl sbcl)
  (portability-warning 'process-close process)
  )

(defun process-status (process)
  "Returns the status of @arg{process}."
  #+(or cmu scl) (ext:process-status process)
  #+ccl (ccl:external-process-status process)
  #+sbcl (sb-ext:process-status process)
  #-(or cmu ccl scl sbcl)
  (portability-warning 'process-status process)
  )

(defun kill-process (process)
  "Kills the given process"
  #+(and sbcl linux)
  (run-program
   "/usr/bin/pkill"
   (list "-P"  (format nil "~D" (sb-impl::process-pid process)))
   :wait t)
  #-(and sbcl linux)
  (portability-warning 'kill-process process)
  )

(defun run-program-output (program args)
  "Returns a list of all lines of the output of @function{RUN-PROGRAM}
when called with the arguments @arg{PROGRAM} and @arg{ARGS}."
  (let ((process (run-program (find-executable program) args
                              :wait nil :input :stream :output :stream
                              :error-output nil)))
    (unwind-protect
	 (with-open-stream (stream (process-output process))
	   (loop for line = (read-line stream nil) while line collect line))
      (process-close process))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Memory
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun memory-usage ()
  "Returns an approximation to the memory used at this moment."
  #+(or cmu scl) (lisp::dynamic-usage)
  #+sbcl (sb-kernel::dynamic-usage)
  #-(or cmu scl sbcl)
  (portability-warning 'memory-usage))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Foreign libraries
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;; TODO: just use cffi everywhere
(defun load-foreign-library (file)
  "Loads the foreign library @arg{file}."
  #+allegro (cl:load file)
  #+lispworks (fli:register-module file)
  #+cmu (sys::load-object-file file)
  #+scl (ext:load-dynamic-object file)
  #+ecl (ffi:load-foreign-library
	 (if (stringp file) file (namestring file)))
  #+sbcl (sb-alien:load-shared-object file)
  #+(and uffi (not (or allegro cmu scl sbcl ecl)))
  (uffi:load-foreign-library file)
  #+(and cffi (not (or allegro cmu scl sbcl ecl uffi)))
  (eval
   (let ((name (gensym)))
     `(progn
        (cffi:define-foreign-library ,name
          (:unix (:or ,file)))
        (cffi:use-foreign-library ,name))))
  #-(or allegro lispworks cmu scl sbcl ecl uffi cffi)
  (portability-warning 'load-foreign-library file))

(defun convert-type (type)
  #+(or allegro cmu scl lispworks sbcl)
  (cond ((null type) ())
        #+lispworks
        ((and (listp type) (eq (car type) '*))
         :lisp-simple-1d-array)
        ((atom type)
         #+(or allegro lispworks) type
         #+(or cmu scl sbcl)
         (if (member type '(:void :char :int :float :double))
             (find-symbol (symbol-name type)
                          #+(or cmu scl) "C-CALL" #+sbcl "SB-ALIEN")
           type)
         )
	(t (cons (convert-type (car type)) (convert-type (cdr type)))))
  #-(or allegro cmu scl lispworks sbcl)
  (portability-warning 'convert-type type)
  )

#+(or allegro lispworks cmu scl sbcl)
(defmacro simplified-def-function ((c-name lisp-name) args &rest keys)
  (declare (ignorable lisp-name))
  (let (fortran-flag)
    (declare (ignorable fortran-flag))
    `(#+(or cmu scl) c-call::def-alien-routine
        #+sbcl sb-alien::define-alien-routine
        #+allegro ff:def-foreign-call
        #+lispworks fli:define-foreign-function
        #+(or allegro) (,lisp-name ,c-name)
        #+(or cmu scl sbcl) ,c-name
        #+(or lispworks) (,lisp-name ,c-name :source)
        #+(or cmu scl sbcl) ,(convert-type (getf keys :returning))
        #+allegro
        ,(mapcar #'(lambda (arg)
                     (destructuring-bind (name type &optional what) arg
                       (case what
                         ((:copy :in-out) (setq fortran-flag t)))
                       (list name (subst :array '* type))))
                 args)
        #+(or cmu scl sbcl)
        ,@(mapcar #'(lambda (arg)
                      (destructuring-bind (name type . rest) arg
                        (list* name (convert-type type) rest)))
                  args)
        #+lispworks
        ,(mapcar #'(lambda (arg)
                     (destructuring-bind (name type &optional what) arg
                       (list name
                             (ecase what
                               ((nil)
                                (if (and (listp type) (eq (car type) '*))
                                    :lisp-simple-1d-array
                                    type))
                               (:in :lisp-simple-1d-array)
                               (:copy (list :reference-pass type))
                               (:in-out (list :reference-return type))))))
                 args)
        #+lispworks ,@(list :result-type (getf keys :returning) :language :ansi-c)
        #+allegro ,@(if fortran-flag '(:convention :fortran))
        #+allegro ,@keys)))

(defmacro def-function (&rest args)
  "Defines a foreign function.  See examples in
@path{alien;src;superlu.lisp}."
  #+(or allegro lispworks cmu scl sbcl) `(simplified-def-function ,@args)
  #+ecl `(ffi:def-function ,@args)
  #+(and uffi (not (or allegro cmu scl sbcl ecl)))
  `(uffi:def-function ,@args)
  #+(and cffi (not (or allegro cmu scl sbcl ecl uffi)))
  (destructuring-bind (name args &key returning &allow-other-keys)
      args
    `(cffi:defcfun ,name ,returning
       ,@(subst :pointer '* args)))
  #-(or allegro lispworks cmu scl sbcl ecl uffi cffi)
  (portability-warning 'define-alien-routine args))

  ;; ensure that vector-sap is working before activating UMFPACK
(declaim (inline vector-sap))
(defun vector-sap (ptr)
   "Returns an array pointer which can be used in a foreign call."
   #+(or allegro lispworks) ptr
   #+(or cmu scl) (system:vector-sap ptr)
   #+ecl nil  ; (ffi:c-inline (ptr) (:object) :pointer-void "(#0)->array.self.ch" :one-liner t)
   #+sbcl (sb-sys:vector-sap ptr)
   #+clisp (ffi::foreign-pointer ptr)
   #-(or allegro lispworks cmu scl sbcl ecl clisp)
   (portability-warning 'vector-sap ptr)
   )

#+lispworks
(defun convert-to-static (arg)
  (cond
    ((not (vectorp arg)) arg)
    ((simple-vector-p arg)
     (let ((n (length arg)))
       (etypecase (aref arg 0)
	 ((complex single-float)
	  (let ((result (make-array (* 2 (length arg)) :element-type 'single-float
				     :initial-element 0.0f0 :allocation :static)))
	    (dotimes (i n result)
	      (let ((x (aref arg i)))
		(declare (type (complex single-float) x))
		(setf (aref result (* 2 i)) (realpart x)
		      (aref result (1+ (* 2 i))) (imagpart x))))))
	 ((complex double-float)
	  (let ((result (make-array (* 2 (length arg)) :element-type 'double-float
				     :initial-element 0.0d0 :allocation :static)))
	    (dotimes (i n result)
	      (let ((x (aref arg i)))
		(declare (type (complex double-float) x))
		(setf (aref result (* 2 i)) (realpart x)
		      (aref result (1+ (* 2 i))) (imagpart x)))))))))
    ((system:staticp arg) arg)
    (t (make-array (length arg) :element-type (array-element-type arg)
		   :allocation :static :initial-contents arg))))

#+lispworks
(defun copy-from-static (arg vec)
  (when (and (vectorp arg) (not (eql arg vec)))
    (cond
      ((simple-vector-p arg)
       (let ((n (length arg)))
	 (etypecase vec
	   ((simple-array single-float (*))
	    (dotimes (i n)
	      (setf (aref arg i)
		    (complex (aref vec (* 2 i)) (aref vec (1+ (* 2 i)))))))
	   ((simple-array double-float (*))
	    (dotimes (i n)
	      (setf (aref arg i)
		    (complex (aref vec (* 2 i)) (aref vec (1+ (* 2 i))))))))))
      (t (fli:replace-foreign-array arg vec)))))

(defun foreign-convert (x)
  (typecase x
    (number x)
    (vector (fl.port:vector-sap x))))

(defun execute-thunk-with-pinned-objects (thunk pinned-objects)
  (if (null pinned-objects)
      (funcall thunk)
      (#+scl ext:with-pinned-object #+sbcl sb-sys:with-pinned-objects
             ((car pinned-objects))
             (execute-thunk-with-pinned-objects thunk (cdr pinned-objects)))))

(defun foreign-call (function &rest args)
  "Ensures a safe environment for a foreign function call, especially so
that no GC changes array pointers obtained by @function{vector-sap}."
  #+ecl (apply function args)
  #+(or allegro ccl) (apply function args)
  #+lispworks
  (let ((transformed-args (mapcar #'convert-to-static args)))
    (multiple-value-prog1
        (apply function transformed-args)
      (mapc #'copy-from-static args transformed-args)))
  #+cmu
  (system:without-gcing
    (apply function args))
  #+(or sbcl scl)
  (#+scl ext:with-pinned-object
   #+sbcl sb-sys:with-pinned-objects (args)
         (apply function (mapcar #'foreign-convert args)))
  ;;(execute-thunk-with-pinned-objects (lambda () (apply function args)) args)
  ;;(apply function args)
  #-(or allegro ccl lispworks cmu scl sbcl ecl)
  (portability-warning 'foreign-call-wrapper)
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
  #+(or cmu scl) (ext:make-weak-pointer obj)
  #+sbcl (sb-ext:make-weak-pointer obj)
  #-(or allegro cmu scl sbcl)
  (portability-warning 'make-weak-pointer obj)
  )

(defun weak-pointer-value (wp)
  "Returns the value of the weak pointer @arg{wp}."
  #+allegro (aref wp 0)
  #+(or cmu scl) (ext:weak-pointer-value wp)
  #+sbcl (sb-ext:weak-pointer-value wp)
  #-(or allegro cmu scl sbcl)
  (portability-warning 'weak-pointer-value wp))

(defun finalize (obj func)
  "Sets up @arg{func} as finalizer for @arg{obj}."
  #+allegro (excl:schedule-finalization
	     obj #'(lambda (obj)
		     (declare (ignore obj))
		     (funcall func)))
  #+(or cmu scl) (ext:finalize obj func)
  #+sbcl (sb-ext:finalize obj func)
  #-(or allegro cmu scl sbcl)
  (portability-warning 'finalize obj func))

(defun gc (&rest args)
  #+(or allegro lispworks cmu scl sbcl) (declare (ignore args))
  #+allegro (excl:gc)
  #+lispworks (hcl::mark-and-sweep 3)
  #+(or cmu scl) (ext:gc)
  #+sbcl (sb-ext:gc)
  #-(or allegro lispworks cmu scl sbcl)
  (portability-warning 'gc args))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Quitting
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun add-exit-hook (hook)
  #+sbcl    (pushnew hook sb-ext:*exit-hooks*)
  #+ccl     (pushnew hook ccl:*lisp-cleanup-functions*)
  #+allegro (pushnew hook sys:*exit-cleanup-forms* :test #'equal)
  #-(or sbcl ccl allegro)
  (portability-warning 'add-exit-hook))
  
(defun quit ()
  "Quits Femlisp."
  #+asdf3 (uiop:quit)
  #+(and lispworks (not asdf3)) (lispworks:quit)
  #-(or asdf3 lispworks) (portability-warning 'quit))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Saving a core and restarting
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun save-femlisp-core-and-die (&optional core-file-name)
  "Saves Femlisp core and quits."
  (unless core-file-name
    (setq core-file-name
          (fl.start:femlisp-pathname
           (format nil "bin/femlisp-~A-core"
                   (string-downcase (lisp-implementation-type))))))
  (format t "Saving ~A~%" core-file-name)
  (setq *package* (find-package :fl.application))
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
  #+clisp (EXT:SAVEINITMEM core-file-name)
  #+(or cmu scl)
  (ext:save-lisp
   core-file-name #+cmu :print-herald #+cmu t
   :init-function
   (lambda ()
     (fl.start::femlisp-banner)
     (setq *package* (find-package :fl.application))
     (lisp::%top-level)))
  #+gcl     (si:save-system core-file-name)
  #+ecl     (quit)
  #+lispworks (hcl:save-image
	       core-file-name :environment nil :restart-function
	       (lambda ()
		 (fl.start::femlisp-banner)
		 (setq *package* (find-package :fl.application))))
  #+sbcl (sb-ext:save-lisp-and-die
          core-file-name :purify t :executable t :toplevel
          (lambda ()
            (fl.start::femlisp-restart)
            (sb-impl::toplevel-init)))
  #-(or allegro lispworks clisp cmu scl sbcl ecl gcl)  ; do nothing
  (portability-warning 'save-femlisp-core-and-die core-file-name)
  ;; we quit in any case
  #+(or)(quit)
  )

#+(or)
(defun save-femlisp-core-and-die (&optional core-file-name (restart-f 'fl.start::femlisp-restart))
  "Saves Femlisp core and quits."
  (unless core-file-name
    (setq core-file-name
          (fl.start:femlisp-pathname
           (format nil "bin/femlisp-~A-core"
                   (string-downcase (lisp-implementation-type))))))
  (format t "Saving ~A~%" core-file-name)
  (setf *package* (find-package :fl.application))
  ;; #+allegro
  ;; (progn
  ;;   (tpl:setq-default *package* (find-package :fl.application))
  ;;   (rplacd (assoc 'tpl::*saved-package*
  ;;                  tpl:*default-lisp-listener-bindings*)
  ;;           'common-lisp:*package*))
  #+asdf3 (push restart-f uiop:*image-restore-hook*)
  #+asdf3 (uiop:dump-image core-file-name :executable t)
  #-asdf3 (portability-warning 'save-femlisp-core-and-die core-file-name)
  #+asdf3 (quit))


;;;; Testing
(defun test-port ()
  (print
   (macroexpand-1
    '(SIMPLIFIED-DEF-FUNCTION ("dnrm2_" DNRM2-)
      ((N :INT :COPY) (X (* :DOUBLE) :IN) (INCX :INT :COPY)) :RETURNING
      :DOUBLE)))
  )
