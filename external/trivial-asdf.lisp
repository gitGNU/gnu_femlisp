;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; trivial-asdf.lisp - Trivial loader of femlisp.asd
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;
;;; Copyright (C) 2005 Nicolas Neuss, University of Heidelberg.
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

(defpackage "ASDF"
  (:use "COMMON-LISP")
  (:export "DEFSYSTEM" "OPERATE")
  (:documentation "This package provides a trivial @code{asdf}-like
loader if @code{asdf} should not be present."))

(in-package :asdf)

(defvar *central-registry* ())
(defvar *systems* ())

(defun store-system (name pathname components)
  (setq *systems* (delete name *systems* :key #'first))
  (push (list name pathname components) *systems*))

(defmacro defsystem (name &key pathname components &allow-other-keys)
  `(store-system ,name ,pathname (quote ,components)))

(defun operate (operate name)
  (declare (ignore operate))
  ;; try to load system file
  (dolist (dir *central-registry*)
    (unless (consp dir)
      (let ((file (concatenate 'string (directory-namestring dir)
			       (string-downcase (symbol-name name))
			       ".asd")))
	(when (probe-file file)
	  (load file)))))
  (let ((system (assoc name *systems*)))
    (unless system (error "System not found."))
    (load-components (directory-namestring (second system))
		     (third system))))

(defun load-components (dir components)
  (dolist (component components)
    (let ((file (getf component :file)))
      (if file
	  (let ((filename (concatenate 'string dir file ".lisp")))
	    (load (or (probe-file (compile-file-pathname filename))
		      (compile-file filename))))
	  (let ((module (getf component :module)))
	    (load-components
	     (let ((new-path (getf component :pathname)))
	       (or (and new-path (directory-namestring new-path))
		   (concatenate 'string dir module "/")))
	     (getf component :components)))))))
      
(pushnew :asdf *features*)