;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;  cdrsys-fe.lisp
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

(in-package "COMMON-LISP-USER")
(defpackage "FL.CDRSYS-FE"
  (:use "COMMON-LISP" "FL.MACROS" "FL.UTILITIES" "FL.MATLISP"
	"FL.MESH" "FL.PROBLEM" "FL.CDRSYS" "FL.DISCRETIZATION")
  (:export)
  (:documentation "This package specializes the finite element
discretization for systems of convection-diffusion-reaction type."))

(in-package "FL.CDRSYS-FE")

#+(or)
(defmethod discretize-locally ((problem <cdrsys-problem>) coeffs vecfe qrule fe-geometry
			       &key local-mat local-rhs local-sol local-u local-v
			       coefficient-parameters &allow-other-keys)
  "Local discretization for a convection-diffusion-reaction system."
  (let ((diffusion-function (getf coeffs 'FL.CDRSYS::DIFFUSION))
	(convection-function (getf coeffs 'FL.CDRSYS::CONVECTION))
	(gamma-function (getf coeffs 'FL.CDRSYS::GAMMA))
	(source-function (getf coeffs 'FL.CDRSYS::SOURCE))
	(reaction-function (getf coeffs 'FL.CDRSYS::REACTION)))
    (error "NYI")))


#+(or)  ; has to be converted
(defmethod essential-boundary-constraints ((problem <pde-problem>) (ansatz-space <ansatz-space>)
					    &key level sol mat rhs (where :surface))
  (declare (ignore sol rhs))
  (assert (or (not (eq where :surface)) mat))
  (let* ((problem (problem ansatz-space))
	 (h-mesh (hierarchical-mesh ansatz-space))
	 (level-skel (if level (cells-on-level h-mesh level) h-mesh))
	 (fe-class (fe-class ansatz-space))
	 (constraints-P (make-ansatz-space-automorphism ansatz-space))
	 (constraints-Q (make-ansatz-space-automorphism ansatz-space))
	 (constraints-rhs (make-ansatz-space-vector ansatz-space)))
    (doskel (cell level-skel)
      (when (ecase where
	      (:refined (refined-p cell h-mesh))
	      (:surface (or (not (refined-p cell h-mesh))
			    (matrix-row mat cell)))
	      (:all t))
	(let* ((coeffs (coefficients-of-cell cell h-mesh problem))
	       (dirichlet-function (getf coeffs 'CONSTRAINT))
	       (cell-key (cell-key cell h-mesh)))
	  (when dirichlet-function
	    (loop with fe = (get-fe fe-class cell)
		  for dof in (fe-dofs fe)
		  for j below (nr-of-inner-dofs fe)
		  for k = (dof-in-vblock-index dof)
		  for ci = (list :local (dof-coord dof)
				 :global (local->global cell (dof-gcoord dof)))
		  do
		  ;; The following is only correct for degrees of freedom of
		  ;; Lagrange type.  Perhaps one should use Hermite finite
		  ;; cells only in the interior?
		  (setf (mref (mref constraints-P cell-key cell-key) k k) 1.0)
		  ;;		       (clear-row mat cell k)
		  ;;		       (setf (mref (mref mat cell cell) k k) 1.0)
		  (setf (vref (vref constraints-rhs cell-key) k)
			(evaluate dirichlet-function ci))
		  )))))
    (values constraints-P constraints-Q constraints-rhs)))

;;; Testing
(defun cdrsys-fe-tests ()
  #+(or) ; something like the following
  (let* ((order 1) (level 2)
	 (problem (cdr-model-problem 1))
	 (h-mesh (uniformly-refined-hierarchical-mesh (domain problem) level))
	 (fedisc (lagrange-fe order)))
    (multiple-value-bind (matrix rhs)
	(discretize-globally problem h-mesh fedisc)
      (getrs (sparse-ldu matrix) rhs)))
  )

(fl.tests:adjoin-test 'cdrsys-fe-tests)




