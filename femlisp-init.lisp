;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; femlisp-init.lisp - Initialization file for Femlisp
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(in-package :cl-user)

(defparameter *femlisp-pathname*
  (make-pathname :directory (pathname-directory *load-pathname*))
  "The pathname for the Femlisp main directory.  This should be the
location of this file when it is loaded.")

(defparameter *femlisp-directory* (namestring *femlisp-pathname)
  "The namestring for @var{*femlisp-pathname*}.")

(let ((directory (pathname-directory *femlisp-pathname*)))
  (setf (logical-pathname-translations "FEMLISP")
	`(("**;*.*.*" 
	   ,(make-pathname :directory `(,@directory :wild-inferiors)
			   :name :wild :type :wild :version :wild)))))

;;; Ensure the presence of external libraries

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
