;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; hole-domain.lisp - Definition of a domain with hole
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

(in-package :application)

(defun n-cube-with-cubic-hole (dim &key (refinements 0))
  "Generates an n-cube-domain with an n-cube hole."
  (let* ((skel1 (copy-skeleton
		 (skeleton-boundary
		  (refcell-refinement-skeleton (n-cube dim) refinements))))
	 (1->2 (nth-value 1 (linearly-transformed-skeleton
			     skel1 :A (scal 0.5 (eye dim))
			     :b (make-double-vec dim 0.25)))))
    (change-class (telescope skel1 1->2) '<domain>)))

(defun n-cell-with-cubic-hole (dim &key (refinements 0))
  "Generates an n-dimensional cell domain with an n-cube hole."
  (identify-unit-cell-faces
   (n-cube-with-cubic-hole dim :refinements refinements)))

(defun n-cube-with-ellipsoidal-hole (dim &key (refinements 0) A)
  "Generates an n-cube-domain with an ellipsoidal hole satisfying (Ax,x)=1
using n-cube patches."
  (let* ((skel1 (copy-skeleton
		 (skeleton-boundary
		  (refcell-refinement-skeleton (n-cube dim) refinements))))
	 (midpoint (make-double-vec dim 0.5))
	 (projection (project-to-ellipsoid midpoint A))
	 (1->2 (nth-value
		1 (transformed-skeleton skel1 :transformation projection))))
    (change-class (telescope skel1 1->2) '<domain>)))

(defun n-cell-with-ellipsoidal-hole (dim &key (refinements 0) A)
  "Generates an n-dimensional cell domain with an ellipsoidal hole."
  (identify-unit-cell-faces
   (n-cube-with-ellipsoidal-hole dim :refinements refinements :A A)))

(defun n-cube-with-n-ball-hole (dim &key (refinements 0) (radius 0.25))
  "Generates an n-cube-domain with an n-ball hole using n-cube patches."
  (n-cube-with-ellipsoidal-hole dim :refinements refinements
				:A (scal (/ (* radius radius)) (eye dim))))

(defun n-cell-with-n-ball-hole (dim &key (radius 0.25) (refinements 0))
  "Generates an n-dimensional cell domain with an n-ball hole."
  (identify-unit-cell-faces
   (n-cube-with-n-ball-hole dim :radius radius :refinements refinements)))

(defun patch-on-inner-boundary-p (patch)
  "Checks if the patch is part of the hole boundary."
  (patch-in-inlay-p patch))


;;; Testing: (test-hole-domain)
(defun test-hole-domain ()
  (check (n-cube-with-cubic-hole 2))
  (check (n-cube-with-n-ball-hole 2))
  (let ((skel (n-cell-with-n-ball-hole 2 :radius 0.3 :refinements 0)))
    (doskel (cell skel)
      (whereas ((id (cell-identification cell skel)))
	(format t "~A~% --> ~A~%" cell id))))
  )

(fl.tests:adjoin-test 'test-hole-domain)