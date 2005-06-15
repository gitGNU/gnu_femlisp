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
	   "SYSTEM-DIFFUSION-PROBLEM")
  (:documentation "This package contains the problem definition of systems
of convection-diffusion-reaction equations."))

(in-package :fl.cdrsys)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; <cdrsys-problem>
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass <cdrsys-problem> (<pde-problem>)
  ((n :reader nr-of-components :initarg :nr-of-components))
  (:documentation "Systems of convection-diffusion-reaction equations.  The
coefficients should be vector-valued functions in this case."))

(defmethod interior-coefficients ((problem <cdrsys-problem>))
  "Coefficients for a CDR system."
  '(DIFFUSION CONVECTION REACTION SOURCE GAMMA CONSTRAINT))
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Systems of diffusion equations (mainly for testing purposes)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun system-diffusion-tensor (&key dim D)
  "Returns a tensor for a system of separate diffusion equations."
  (let* ((nr-comps (length D))
	 (tensor (make-array (list nr-comps nr-comps) :initial-element nil)))
    (dotimes (k nr-comps)
      (dotimes (l nr-comps)
	(setf (aref tensor k l)
	      (if (= k l)
		  (scal! (aref D k) (eye dim))
		  (zeros dim)))))
    tensor))

(defun system-diffusion-problem (domain &key nr-comps D source)
  (let ((dim (dimension domain)))
    (make-instance
     '<cdrsys-problem>
     :domain domain
     :patch->coefficients
     #'(lambda (patch)
	 (if (member-of-skeleton? patch (domain-boundary domain))
	     (list 'FL.CDRSYS::CONSTRAINT (constraint-coefficient nr-comps 1))
	     (list 'FL.CDRSYS::DIFFUSION
		   (constant-coefficient
		    (system-diffusion-tensor :dim dim :D D))
		   'FL.CDRSYS::SOURCE
		   (ensure-coefficient source)))))))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Generation of standard cdrsys problems
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; Testing: (test-cdrsys)

(defun test-cdrsys ()
  (make-instance '<cdrsys-problem> :domain (n-cube-domain 2)
		 :patch->coefficients (constantly nil) :nr-of-components 1))

(fl.tests:adjoin-test 'test-cdrsys)

