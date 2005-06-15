;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; inlay-domain.lisp - Domain definitions for an inlay domain
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

(in-package :fl.domains)

(defun n-cube-with-cubic-inlay (dim)
  "Generates an n-cube-domain with an n-cube inlay."
  (multiple-value-bind (inlay)
      (linearly-transformed-skeleton (skeleton (n-cube dim))
				     :A (scal 0.5 (eye dim))
				     :b (make-double-vec dim 0.25))
    (change-class
     (skel-add! inlay (n-cube-with-cubic-hole dim))
     '<domain>)))

(defun n-cell-with-cubic-inlay (dim)
  "Generates an n-dimensional cell domain with an n-cube hole."
  (identify-unit-cell-faces (n-cube-with-cubic-inlay dim)))

(defun n-cube-with-ball-inlay (dim &key (radius 0.25))
  "Generates an n-cube-domain with an n-ball inlay using n-cube patches."
  (let* ((outer-skel (skeleton-boundary (skeleton (n-cube dim))))
	 (midpoint (make-double-vec dim 0.5))
	 (projection (project-to-sphere midpoint radius))
	 (outer->middle
	  (nth-value
	   1 (transformed-skeleton
	      outer-skel :transformation projection))))
    (multiple-value-bind (center-block unit-cell->center-block)
	(let ((factor (/ radius (sqrt dim))))
	  (linearly-transformed-skeleton (skeleton (n-cube dim))
					 :A (scal factor (eye dim))
					 :b (scal (- 1.0 factor) midpoint)))
      (let ((center-skin (skeleton-boundary center-block))
	    (center->middle (make-hash-table)))
	;; fill center->middle table
	(dohash ((outer-cell middle-cell) outer->middle)
	  (let ((center-cell (gethash outer-cell unit-cell->center-block)))
	    (assert center-cell) 
	    (setf (gethash center-cell center->middle)
		  middle-cell)))
	(skel-add! center-block (telescope center-skin center->middle))
	(skel-add! center-block (telescope outer-skel outer->middle))
	(change-class center-block '<domain>)))))

(defun n-cell-with-ball-inlay (dim &key (radius 0.25))
  "Generates an n-dimensional cell domain with an n-ball inlay."
  (identify-unit-cell-faces
   (n-cube-with-ball-inlay dim :radius radius)))

(defun patch-in-inlay-p (patch)
  "Checks if the patch is part of the inlay including its boundary."
  (let ((corners (corners patch)))
    (every #'(lambda (corner)
	       (every #'(lambda (coord) (< 0.0 coord 1.0)) corner))
	   corners)))


;;; Testing
(defun test-inlay-domain ()
  (n-cube-with-cubic-inlay 2)
  (let* ((domain (n-cell-with-ball-inlay 2 :radius 0.3))
	 (chars (domain-characteristics domain)))
    (assert (and (getf chars :exact) (getf chars :curved)))
    (doskel (cell domain)
      (whereas ((id (cell-identification cell domain)))
	(format t "~A~% --> ~A~%" cell id)))
    (fl.mesh::check-identification domain))
  )

;;; (test-inlay-domain)
(fl.tests:adjoin-test 'test-inlay-domain)

