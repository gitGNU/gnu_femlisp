;;; -*- mode: lisp; fill-column: 64; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; setup.lisp - Femlisp setup file
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;
;;; Copyright (C) 2003-2005 Nicolas Neuss, University of Heidelberg.
;;; Copyright (C) 2006-2008 Nicolas Neuss, University of Karlsruhe.
;;; Copyright (C) 2010-     Nicolas Neuss, University Erlangen-Nuremberg.
;;; All rights reserved.
;;; 
;;; Redistribution and use in source and binary forms, with or without
;;; modification, are permitted provided that the following conditions are
;;; met:
;;; 
;;; 1. Redistributions of source code must retain the above
;;; copyright notice, this list of conditions and the following
;;; disclaimer.
;;; 
;;; 2. Redistributions in binary form must reproduce the above
;;; copyright notice, this list of conditions and the following
;;; disclaimer in the documentation and/or other materials
;;; provided with the distribution.
;;; 
;;; THIS SOFTWARE IS PROVIDED ``AS IS'' AND ANY EXPRESS OR
;;; IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
;;; IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
;;; PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
;;; AUTHOR, THE UNIVERSITIES HEIDELBERG OR KARLSRUHE, OR OTHER
;;; CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
;;; SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
;;; NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
;;; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
;;; HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
;;; CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
;;; OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
;;; SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(in-package :common-lisp-user)

(defpackage "FL.START"
  (:use "COMMON-LISP")
  (:export "*FEMLISP-VERSION*" "FEMLISP-HERALD" "FEMLISP-BANNER" "FEMLISP-PATHNAME"
           "*DX-PATH*" "*GNUPLOT-PATH*"
           "*TETGEN-PATH*" "*IMAGES-DIRECTORY*" "*MESHES-DIRECTORY*"
           "*BLAS-LIBRARY*" "*LAPACK-LIBRARY*" "*SUPERLU-LIBRARY*" "*UMFPACK-LIBRARY*")
  (:documentation "This package contains some routines called
during initialization of Femlisp."))

(in-package :fl.start)

(defparameter *femlisp-version* "2.0.1")

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Setup the logical host "FEMLISP"
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun femlisp-pathname (&optional (namestring ""))
  "Get pathname relative to the Femlisp directory."
  (asdf:system-relative-pathname
   :femlisp (concatenate 'string "../" namestring)))

(defvar *femlisp-pathname*
  (femlisp-pathname)
  "The pathname for the Femlisp main directory.  This should be the
location of this file when it is loaded.")

(defvar *femlisp-directory* (namestring *femlisp-pathname*)
  "The namestring for @var{*femlisp-pathname*}.")

;;; earlier code might depend on the Femlisp logical host being defined
(let ((directory (pathname-directory *femlisp-directory*)))
  (setf (logical-pathname-translations "FEMLISP")
        `((#+gcl "**;*.*.*" #-gcl "**;*.*" 
                 ,(make-pathname :directory `(,@directory :wild-inferiors)
                                 :name :wild :type :wild :version :wild)))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Configuration
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Implementation-dependent configurations
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; since Femlisp usually works with lots of memory, GC messages
;;; are annoying

#+allegro (setq excl:*global-gc-behavior* :auto)
#+cmu (setq extensions:*gc-verbose* nil)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; we want to work generally with double float numbers
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; We drop this here for not interfering with other systems.
;;; Instead this is replaced by loading Femlisp inside an
;;; environment where this variable is bound dynamically.  See
;;; @file{system/femlisp.asd}.

;; (setq *READ-DEFAULT-FLOAT-FORMAT* 'double-float)

