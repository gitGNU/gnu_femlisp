;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; cdrsys.lisp
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

(in-package :cl-user)

(defpackage "FL.CDRSYS"
  (:use "COMMON-LISP" "FL.MATLISP" "FL.MACROS" "FL.UTILITIES"
	"FL.MESH" "FL.PROBLEM" "FL.ELLSYS")
  (:export "<CDRSYS-PROBLEM>" "CDRSYS-MODEL-PROBLEM")
  (:documentation "This package contains some definitions for systems of
convection-diffusion-reaction equations.  These are a special case of
general elliptic systems defined in @path{ellsys.lisp}."))

(in-package :fl.cdrsys)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; <cdrsys-problem>
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass <cdrsys-problem> (<ellsys-problem>)
  ()
  (:documentation "Problem class for a system of
convection-diffusion-reaction equations."))

(defun cdrsys-model-problem (domain ncomps &key (diffusion nil diffusion-p)
			     (source nil source-p) (dirichlet nil dirichlet-p)
			     convection reaction initial evp properties)
  "Generates a system of convection-diffusion-reaction equations.  Defaults
are identity diffusion, right-hand-side equal 1, and Dirichlet zero
boundary conditions for each component.  Ordinary function are converted
into coefficient functions depending on a global coordinate.  The first
argument can be either a domain or an integer n which is interpreted as the
n-dimensional unit cube."
  (setq domain (ensure-domain domain))
  (let ((dim (dimension domain)))
    (ellsys-model-problem
     domain ncomps
     :derived-class '<cdrsys-problem>
     :a (if diffusion-p
	    diffusion
	    (isotropic-diffusion
	     dim (coerce (loop repeat ncomps collect 1.0) 'vector)))
     :b convection
     :c reaction
     :f (if source-p
	    source
	    (coerce (loop repeat ncomps collect
			  (if evp (zeros 1) (eye 1)))
		    'vector))
     :properties properties
     :dirichlet (if dirichlet-p
		    dirichlet
		    (constraint-coefficient ncomps 1))
     :initial initial :evp evp)))


;;; Testing: (test-cdrsys)

(defun test-cdrsys ()
  (cdrsys-model-problem 1 1)
  (describe (cdrsys-model-problem 1 1))
  )

(fl.tests:adjoin-test 'test-cdrsys)

