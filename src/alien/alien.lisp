;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; alien.lisp - interface to external programs
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;
;;; Copyright (C) 2004 Nicolas Neuss, University of Heidelberg.
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

(defpackage "FL.ALIEN"
  (:use "COMMON-LISP" "FL.MACROS" "FL.PORT")
  (:documentation "This package loads some foreign libraries and does set
up a Lisp interface for it."))

(in-package "FL.ALIEN")

(defvar *foreign-code-loaders* ()
  "List of functions for loading foreign code.  Each function of this list
should be executed when starting a saved core.")

(defun reload-foreign-code ()
  (dolist (func *foreign-code-loaders*)
    (funcall func)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Foreign libraries
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; SuperLU

(defun load-superlu-library ()
  "Loads foreign code for accessing SuperLU."
  (whereas ((superlu-glue (probe-file #p"femlisp:interface;superlu.so")))
    (load-shared-object superlu-glue)
    (pushnew :superlu *features*)))

;; execute this function immediately
(load-superlu-library)
;; and register it for reinitialization
(pushnew 'load-superlu-library *foreign-code-loaders*)

;;; UMFPACK/AMD

(defun load-umfpack-library ()
  "Loads foreign code for accessing UMFPACK."
  (whereas ((umfpack-glue (probe-file #p"femlisp:interface;umfpack.so")))
    (load-shared-object umfpack-glue)
    (pushnew :umfpack *features*)))


;; load it immediately
(load-umfpack-library)
;; and register it for reinitialization
(pushnew 'load-umfpack-library *foreign-code-loaders*)

#+cmu  
(defun reload-global-table ()
  "Function mailed by Eric Marsden to cmucl-help at 18.9.2003."
  (loop :for lib-entry in sys::*global-table*
	:for (sap . lib-path) = lib-entry
        :when lib-path :do
        (let ((new-sap (sys::dlopen (namestring lib-path)
                                    (logior sys::rtld-now sys::rtld-global))))
          (when (zerop (sys:sap-int new-sap))
            (error "Couldn't open library ~S: ~S" lib-path (sys::dlerror)))
          (setf (car lib-entry) new-sap)))
  (alien:alien-funcall (alien:extern-alien "os_resolve_data_linkage"
					   (alien:function c-call:void))))

#+cmu
(pushnew 'reload-global-table ext:*after-save-initializations*)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Utility function (working around a bug in CMUCL)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmacro enumerate (name &rest items)
  (declare (ignore name))
  (cons 'progn
	(loop for i from 0 and item in items collect
	      `(defconstant ,item ,i))))

(defun direct-solver-test (solver)
  "Tests @arg{solver} on the simple example from the SuperLU User
Guide."
  (let* ((m 5) (n 5) (nnz 12) (nrhs 1)
	 (s 19.0) (u 21.0) (p 16.0) (e 5.0) (r 18.0) (l 12.0)
	 (cs (make-array (1+ n) :element-type '(signed-byte 32)
			 :initial-contents '(0 3 6 8 10 12)))
	 (ri (make-array nnz :element-type '(signed-byte 32)
			 :initial-contents '(0 1 4 1 2 4 0 2 0 3 3 4)))
	 (sa (make-array nnz :element-type 'double-float
			 :initial-contents (list s l l u l l u p u e u r)))
	 (rhs (make-array (* m nrhs) :element-type 'double-float
			  :initial-element 1.0))
	 (sol (make-array (* m nrhs) :element-type 'double-float
			  :initial-element 0.0)))
    (values (funcall solver m n nnz cs ri sa nrhs rhs sol) sol)))



