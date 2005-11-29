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
	"FL.MESH" "FL.PROBLEM")
  (:export "<CDRSYS-PROBLEM>"
	   "SYSTEM-DIFFUSION-TENSOR"
	   "CDRSYS-MODEL-PROBLEM")
  (:documentation "This package contains the problem definition of systems
of convection-diffusion-reaction equations.  All coefficients are vectors
with entries corresponding to the scalar case."))

(in-package :fl.cdrsys)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; <cdrsys-problem>
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass <cdrsys-problem> (<pde-problem>)
  ((ncomps :reader nr-of-components :initarg :nr-of-components))
  (:documentation "Systems of convection-diffusion-reaction equations.  The
coefficients should be vector-valued functions in this case."))

(defmethod interior-coefficients ((problem <cdrsys-problem>))
  "Coefficients for a CDR system."
  '(DIFFUSION CONVECTION REACTION SOURCE GAMMA CONSTRAINT))
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Generation of standard cdrsys problems
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun vector-diffusion (dim ncomps value)
  (constant-coefficient
   (make-array ncomps :initial-element (scal value (eye dim)))))

(defun cdrsys-model-problem (dim/domain ncomps &key (diffusion nil diffusion-p)
			     (source nil source-p) (dirichlet nil dirichlet-p)
			     gamma convection reaction initial evp properties)
  "Generates a system of convection-diffusion-reaction equations.  Defaults
are identity diffusion, right-hand-side equal 1, and Dirichlet zero
boundary conditions for each component.  Ordinary function are converted
into coefficient functions depending on a global coordinate.  The first
argument can be either a domain or an integer n which is interpreted as the
n-dimensional unit cube."
  (let* ((domain (if (numberp dim/domain)
		     (n-cube-domain dim/domain)
		     dim/domain))
	 (dim (dimension domain)))
    ;; set default values
    (unless diffusion-p
      (setq diffusion (vector-diffusion dim ncomps 1.0)))
    (unless source-p
      (setq source (make-array ncomps :initial-element (ensure-matlisp (if evp 0.0 1.0)))))
    (unless dirichlet-p (setq dirichlet (constraint-coefficient ncomps 1)))
    (apply #'fl.amop:make-programmatic-instance
	   (cond (initial '(<time-dependent-problem> <cdrsys-problem>))
		 (evp '(<cdrsys-problem> <evp-mixin>))
		 (t '<cdrsys-problem>))
	   :nr-of-components ncomps
	   :properties properties
	   :domain domain :patch->coefficients
	   `((:external-boundary
	      ,(when dirichlet
		     `(FL.CDRSYS::CONSTRAINT ,(ensure-coefficient dirichlet))))
	     (:d-dimensional
	      ,(append
		(when diffusion
		  `(FL.CDRSYS::DIFFUSION ,(ensure-coefficient diffusion)))
		(when source
		  `(FL.CDRSYS::SOURCE ,(ensure-coefficient source)))
		(when convection
		  `(FL.CDRSYS::CONVECTION ,(ensure-coefficient convection)))
		(when reaction
		  `(FL.CDRSYS::REACTION ,(ensure-coefficient reaction)))
		(when gamma
		  `(FL.CDRSYS::GAMMA ,(ensure-coefficient gamma)))
		(when initial
		  `(FL.CDRSYS::INITIAL ,(ensure-coefficient initial))))))
	   (append
	    (when evp (destructuring-bind (&key lambda mu) evp
			(list :lambda lambda :mu mu))))
	   )))


;;; Testing: (test-cdrsys)

(defun test-cdrsys ()
  (make-instance '<cdrsys-problem> :domain (n-cube-domain 2)
		 :patch->coefficients (constantly nil) :nr-of-components 1))

(fl.tests:adjoin-test 'test-cdrsys)

