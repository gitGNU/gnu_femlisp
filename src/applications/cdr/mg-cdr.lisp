;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; mg-cdr.lisp - Solving CDR problems with multigrid
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;
;;; Copyright (C) 2003- Nicolas Neuss, University of Heidelberg.
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

(in-package :fl.application)

(defun test-v-cycle-convergence (dim order level &key smoother (output t) cr-max)
  "Solves with a V-cycle and prints the average convergence rate.  If
cr-max is provided, it is checked if the convergence rate is smaller than
this value."
  (let* ((problem (cdr-model-problem dim))
	 (cs (geometric-cs
	      :gamma 1 :pre-steps 1 :pre-smooth smoother :post-steps 0
	      :base-level (1- level)))
	 (solver (make-instance '<linear-solver> :iteration cs
				:success-if '(or (< :defnorm 1.0e-12) (> :step 20))))
	 (mesh (uniformly-refined-hierarchical-mesh (domain problem) level))
	 (fedisc (lagrange-fe order)))
    (multiple-value-bind (A b)
	(discretize-globally problem mesh fedisc)
      (setq *result*
	    (solve (blackboard :matrix A :rhs b :solver solver :output output)))))
  (let ((cr (getbb (getbb *result* :report) :convergence-rate)))
    (format t "CR(~D,~D,~D)=~F" dim order level cr)
    (when cr-max (assert (< cr cr-max)))))

(defun mg-cdr-tests ()
  (time (test-v-cycle-convergence
	 3 1 2 :smoother (make-instance '<gauss-seidel>) :cr-max 0.15))
  (time (test-v-cycle-convergence
	 2 4 3 :smoother (make-instance '<gauss-seidel>) :cr-max 0.25))
  (time (test-v-cycle-convergence
	 2 4 3 :smoother (geometric-ssc) :cr-max 0.06))
  (time (test-v-cycle-convergence
	 2 1 3 :smoother (make-instance '<gauss-seidel>) :cr-max 0.15))
  )

(fl.tests:adjoin-test 'mg-cdr-tests)