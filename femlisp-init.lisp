;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; femlisp-init.lisp - Initialization file for Femlisp
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(in-package :cl-user)

(defparameter *cl-home*
  (concatenate 'string
	       #+cmu (cdr (assoc :CL_HOME ext:*environment-list*))
	       #+sbcl (posix-getenv "CL_HOME") "/"))

(defparameter *cl-home-pathname* (pathname *cl-home*))

(let ((directory (pathname-directory *cl-home-pathname*)))
  (setf (logical-pathname-translations "CL")
	`(("**;*.*.*" 
	   ,(make-pathname :directory `(,@directory :wild-inferiors)
			   :name :wild :type :wild :version :wild)))))

#+cmu
(progn
  (setq extensions:*gc-verbose* nil)
  (load #p"cl:lisp;asdf"))

#+sbcl
(progn
  (require 'asdf))

(setf (logical-pathname-translations "FEMLISP")
      `(("**;*.*.*"
	 ,(make-pathname :host "CL" 
			 :directory '(:absolute "FEMLISP" :wild-inferiors)
			 :name :wild :type :wild :version :wild))))


(push #p"cl:cl-ppcre" asdf::*central-registry*)
(push (translate-logical-pathname #p"femlisp:") asdf::*central-registry*)
