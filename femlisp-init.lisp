;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; femlisp-init.lisp - Initialization file for Femlisp
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(in-package :cl-user)

(defparameter *femlisp-directory*
  (concatenate 'string
	       #+cmu (cdr (assoc :FEMLISP_DIR ext:*environment-list*))
	       #+sbcl (posix-getenv "FEMLISP_DIR") "/"))

(defparameter *femlisp-pathname* (pathname *femlisp-directory*))

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

#+sbcl
(progn
  (require 'asdf))

#-infix
(load "femlisp:external;infix.cl")

#-cl-ppcre
(progn (load "femlisp:external;cl-ppcre;load" :verbose nil)
       (pushnew :cl-ppcre *features*))

(push (translate-logical-pathname #p"femlisp:") asdf::*central-registry*)
