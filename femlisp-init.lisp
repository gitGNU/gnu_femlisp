;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; femlisp-init.lisp - Initialization file for Femlisp
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;; GCL specialities

(in-package :cl-user)

(defparameter *femlisp-directory*
  (concatenate 'string
	       #+cmu (cdr (assoc :FEMLISP_DIR ext:*environment-list*))
	       #+sbcl (posix-getenv "FEMLISP_DIR")
	       #+gcl "/home/neuss/CL-HOME/femlisp"
	       "/"))

(defparameter *femlisp-pathname* (pathname *femlisp-directory*))

#-gcl
(let ((directory (pathname-directory *femlisp-pathname*)))
  (setf (logical-pathname-translations "FEMLISP")
	`(("**;*.*.*" 
	   ,(make-pathname :directory `(,@directory :wild-inferiors)
			   :name :wild :type :wild :version :wild)))))

;;; Ensure the presence of other libraries

#+cmu
(progn
  (setq extensions:*gc-verbose* nil)
  #-asdf (load #p"femlisp:external;asdf"))

#+sbcl (require 'asdf)

#+gcl
(progn  ; should be dropped when gcl has logical pathnames
  (load (compile-file (concatenate 'string *femlisp-directory* "external/infix.cl")))
  (load (compile-file (concatenate 'string *femlisp-directory* "external/defsystem.lisp"))))

#-infix (load "femlisp:external;infix.cl")
#-(or asdf mk-defsystem) (load "femlisp:external;defsystem.lisp")

#+asdf
(push *femlisp-pathname* asdf::*central-registry*)
#+mk-defsystem
(push *femlisp-pathname* mk::*central-registry*)
