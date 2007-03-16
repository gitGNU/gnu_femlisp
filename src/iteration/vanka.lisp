;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; vanka.lisp
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

(in-package :fl.geomg)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; <vanka>
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass <vanka> (<geometric-ssc>)
  ()
  (:documentation "Vanka-like smoother for @math{Q^{k+1}/Q^k}
discretizations of Navier-Stokes."))

(defun extended-block (asa keys)
  "Collect the next layer around keys in the matrix graph of asa."
  (let ((inner (make-hash-table))
	(outer (make-hash-table)))
    ;; set inner hash table
    (loop for key across keys do (setf (gethash key inner) t))
    ;; if it is a standard non-slave block extend it by the outer layer
    (unless (slave-or-dirichlet-dof-p (aref keys 0) asa)
      (loop for key across keys do
	    (for-each-key-in-row
	     #'(lambda (col-key)
		 (unless (gethash col-key inner)
		   (setf (gethash col-key outer) t)))
	     asa key)))
    ;; finally return the block together with the components
    (let ((extended-keys (concatenate 'vector keys (hash-table-keys outer)))
	  (as (ansatz-space asa)))
      (values
       extended-keys
       (vector-map
	#'(lambda (key)
	    (let ((fe (get-fe as (representative key))))
	      (if (gethash key inner)
		  (cons 0 (nr-of-inner-dofs fe))
		  (let ((subcell-offsets (getf (properties fe) 'SUBCELL-OFFSETS)))
		    (cons 0 (aref (vector-last subcell-offsets) 0))))))
	extended-keys)))))

(defmethod setup-blocks ((vanka <vanka>) (asa <ansatz-space-automorphism>))
  "Adds to the usual vertex centered blocks all surrounding velocity
degrees of freedom."
  (let ((inner-blocks (call-next-method)) ; blocks from vertex-centered BGS
	(result-blocks ())
	(result-components ()))
    (dolist (inner inner-blocks)
      (multiple-value-bind (keys components)
	  (extended-block asa inner)
	(push keys result-blocks)
	(push components result-components)))
    (dbg :iter "Vanka blocks: ~A" result-blocks)
    (values result-blocks result-components)))


