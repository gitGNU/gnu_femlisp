;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;  cdr.lisp
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
(defpackage "CDR"
  (:use "COMMON-LISP" "MATLISP" "MACROS" "UTILITIES" "MESH" "PROBLEM")
  (:export "<CDR-PROBLEM>" "MAP-DOMAIN-TO-CDR-PROBLEM"
	   "STANDARD-CDR-PROBLEM" "SCALAR-DIFFUSION" "IDENTITY-DIFFUSION-TENSOR"
	   "LAPLACE-TEST-PROBLEM-ON-DOMAIN" "LAPLACE-TEST-PROBLEM"))
(in-package :cdr)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; <cdr-problem>
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass <cdr-problem> (<problem>)
  ()
  (:documentation "Convection-diffusion-reaction problem."))

(defmethod interior-coefficients ((problem <cdr-problem>))
  "Interior coefficients for a CDR problem."
  '(DIFFUSION CONVECTION SOURCE REACTION GAMMA FE-RHS))

(defmethod boundary-coefficients ((problem <cdr-problem>))
  "Boundary coefficients for a CDR problem."
  '(DIRICHLET))

(defmethod self-adjoint-p ((problem <cdr-problem>))
  (doskel (patch (domain problem))
    (when (getf (funcall (patch->coefficients problem) patch) 'CONVECTION)
      (return-from self-adjoint-p (values NIL T))))
  (values T T))

(defmethod dual-problem ((problem <cdr-problem>) cell->rhs)
  "Dual problem of a cdr problem.  To be improved, at the moment only for
selfadjoint problems with functional of interest being the load
functional."
  (make-instance
   '<cdr-problem> :domain (domain problem)
   :patch->coefficients
   #'(lambda (cell)
       ;; check that problem is self adjoint
       (let ((coeffs (copy-seq (funcall (patch->coefficients problem) cell))))
	 (when (getf coeffs 'DIRICHLET)
	   (setf (getf coeffs 'DIRICHLET) *cf-constantly-0.0d0*))
	 (assert (not (getf coeffs 'CONVECTION)))  ; better: change to negative
	 (unless (eq cell->rhs :load-functional)
	   (let ((dual-rhs (funcall cell->rhs cell)))
	     (setf (getf coeffs 'SOURCE) (getf dual-rhs 'SOURCE))
	     (setf (getf coeffs 'GAMMA) (getf dual-rhs 'GAMMA))
	     (setf (getf coeffs 'FE-RHS) (getf dual-rhs 'FE-RHS))))
	 coeffs))
   :multiplicity (multiplicity problem)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Generation of standard cdr problems
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun standard-cdr-problem (domain &key diffusion gamma convection reaction source
			     (dirichlet *cf-constantly-0.0d0*))
  "Generates a simple cdr problem.  Passing ordinary functions instead of
coefficient functions treats them as x-dependent coefficients."
  (flet ((ensure-coefficient (obj)
	   (if (functionp obj) (function->coefficient obj) obj)))
    (let ((dim (dimension domain))
	  (bdry (skeleton-boundary domain)))
      (make-instance
       '<cdr-problem> :domain domain
       :patch->coefficients
       #'(lambda (cell)
	   (cond ((member-of-skeleton? cell bdry)
		  (list 'DIRICHLET dirichlet))
		 ((= (dimension cell) dim)
		  (nconc
		   (and diffusion (list 'DIFFUSION (ensure-coefficient diffusion)))
		   (and convection (list 'CONVECTION (ensure-coefficient convection)))
		   (and source (list 'SOURCE (ensure-coefficient source)))
		   (and reaction (list 'REACTION (ensure-coefficient reaction)))
		   (and gamma (list 'GAMMA (ensure-coefficient gamma)))))
		 (t nil)))))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Generation of laplace test problems
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun scalar-diffusion (dim value)
  (make-instance '<coefficient>
		 :input *empty-coefficient-input*
		 :eval (constantly (scal value (eye dim)))))

(defun identity-diffusion-tensor (dim)
  (scalar-diffusion dim 1.0d0))

(defun laplace-test-problem-on-domain (domain)
  (standard-cdr-problem
   domain
   :diffusion (identity-diffusion-tensor (dimension domain))
   :source *cf-constantly-1.0d0*))

(defun laplace-test-problem (dim)
  "Generates the problem
$$ -\Delta u = 1 in [0,1]^d$$
with Dirichlet bc."
  (laplace-test-problem-on-domain (n-cube-domain dim)))

;;; Testing: (test-cdr)

(defun test-cdr ()
  (check (domain (laplace-test-problem 2)))
  (check (domain (laplace-test-problem 2)))
  (let* ((domain (n-cell-domain 1))
	 (problem (laplace-test-problem-on-domain domain)))
    (problem-info problem))
  )

(tests:adjoin-femlisp-test 'test-cdr)

