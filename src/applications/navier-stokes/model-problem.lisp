;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; model-problem.lisp - Model problems for Navier-Stokes
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

(defvar *stokes-demo*
  (make-demo :name "Stokes"
	     :short "Solve the Stokes equation on different domains"))
(adjoin-demo *stokes-demo* *equation-demo*)

(defun stokes-driven-cavity-demo (dim order levels &key plot output)
  "Solves the driven cavity problem for Stokes' equation."
  (setf (getf *strategy-output* :observe) nil)
  (defparameter *result*
    (solve-with
     (make-instance
      '<fe-strategy> :fe-class (navier-stokes-lagrange-fe order dim 1)
      :estimator (make-instance '<projection-error-estimator>)
      :indicator (make-instance '<largest-eta-indicator> :fraction 0.5 :from-level levels)
      :appraise (stop-if :nr-levels>= levels)
      :solver
      #+(or)
      (make-instance '<special-solver> :solver-function
		     #'(lambda (&key matrix rhs &allow-other-keys)
			 (m* (sparse-ldu matrix :ordering (boundary-last-order matrix))
			     rhs)))
      #-(or)  ; no local mg here
      (make-instance
       '<linear-solver> :iteration
       (let ((smoother (make-instance '<vanka>)))
	 (geometric-cs
	  :coarse-grid-iteration
	  (make-instance '<multi-iteration> :nr-steps 1 :base smoother)
	  :pre-steps 1 :pre-smooth smoother
	  :post-steps 1 :post-smooth smoother
	  :gamma 2))
       :success-if `(or (> :step 10) (< :defnorm 1.0e-10))
       :failure-if `(and (> :step 2) (> :step-reduction 0.9))
       :output (eq output :all))
      :output t)
     (driven-cavity dim)
     :base-level 0))
  (when plot
    (let ((solution (getf *result* :solution)))
      ;; plot components of cell solution tensor
      (dotimes (i (1+ dim))
	(plot solution :component i :depth 2)
	(sleep 1.0)))))

;;; Testing: (stokes-driven-cavity-demo 2 4 4 :output :all)

(let* ((dim 2) (order 3) (levels 4)
       (demo (make-demo
	      :name "driven-cavity"
	      :short "Solves the driven cavity problem for Stokes' equation"
	      :long
	      (format nil "~A~%~%Parameters: dim=~D, order=~D, levels=~D~%~%"
		      (documentation 'stokes-driven-cavity-demo 'function)
		      dim order levels)
	      :execute (lambda () (stokes-driven-cavity-demo dim order levels :plot t)))))
    (adjoin-demo demo *stokes-demo*))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Testing
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; Driven cavity

;;; Multigrid works!
#+(or)
(defparameter *result*
  (time
   (let* ((dim 2) (level 3) (order 3) (delta 1)
	  (problem (driven-cavity dim :smooth-p nil))
	  (h-mesh (uniformly-refined-hierarchical-mesh (domain problem) level))
	  (fedisc (navier-stokes-lagrange-fe order dim delta)))
     (multiple-value-bind (matrix rhs)
	(discretize-globally problem h-mesh fedisc)
       (let* ((smoother (make-instance '<ns-vanka>))
	      (cs (geometric-cs
		   :gamma 1 :base-level 1 :coarse-grid-iteration
		   (make-instance '<multi-iteration> :base smoother :nr-steps 4)
		   :pre-steps 1 :pre-smooth smoother
		   :post-steps 1 :post-smooth smoother)))
	 (nth-value 1 (linsolve matrix rhs :output t :iteration cs :threshold 1.0e-10)))
	 ))))

;;(plot (getf *result* :res) :component 0)
;;(plot (getf *result* :solution) :component 0 :depth 3)
;;(plot (mesh (getf *result* :res)))

;;; Periodic cavity
#+(or)
(defparameter *result*
  (time
   (let* ((dim 2) (level 3) (order 1) (delta 1)
	  (problem (periodic-cavity dim))
	  (h-mesh (uniformly-refined-hierarchical-mesh (domain problem) level))
	  (fedisc (navier-stokes-lagrange-fe order dim delta)))
     (multiple-value-bind (matrix rhs)
	 (discretize-globally problem h-mesh fedisc)
       (let ((order (boundary-last-order matrix)))
	 #+(or)(display matrix :order order)
	 #+(or)(show rhs)
	 #+(or)(show constraints-P)
	 #+(or)(show constraints-r)
	 #+(or)(m* (sparse-ldu matrix	; :diagonal-inverter #'pressure-diagonal-inverter
			       :ordering order) rhs)
	 #-(or)
	 (let* ((smoother (make-instance '<ns-vanka>))
		(cs (geometric-cs
		     :base-level 2 :coarse-grid-iteration
		     (make-instance '<multi-iteration> :base smoother :nr-steps 2)
		     :gamma 1
		     :pre-steps 1 :pre-smooth smoother
		     :post-steps 1 :post-smooth smoother)))
	   #+(or) ; sieht richtig aus fuer Ordnung 1 und 2
	   (multiple-value-bind (blocks blocks-ranges)
	       (setup-blocks smoother matrix)
	     (loop for block in blocks and ranges in blocks-ranges do
		   (format t "*********~%")
		   (loop for key across block and range across ranges do
			 (format t "~A : ~A~%" key range))))
	   #-(or)
	   (linsolve matrix rhs :output t :iteration cs :maxsteps 20)
	   ))))))


