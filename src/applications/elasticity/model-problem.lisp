;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; model-problem.lisp - Model problems for linear elasticity
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

(in-package :fl.application)

(defun test-elasticity-model-problem ()

;;; Linear elasticity problem
(defparameter *result*
  (time
   (let ((dim 2))
     (solve
      (blackboard
       :problem
       (standard-elasticity-problem
	(n-cube-domain dim) :lambda 1.0 :mu 1.0
	:force (constant-coefficient (make-array dim :initial-element (ones 1))))
       :output t :success-if '(> :nr-levels 2))))))
(plot (getbb *result* :solution) :component 0)

;;; Test for diffusion problem
(defparameter *result*
  (time
   (let ((dim 2) (nr-comps 1))
     (solve
      (blackboard
       :problem
       (system-diffusion-problem
	(n-cube-domain dim) :D 1.0 :nr-comps nr-comps
	:force (constant-coefficient
		(make-array nr-comps :initial-element (eye 1))))
       :fe-class (lagrange-fe 3 :nr-comps 1)
       :output t :success-if '(> :nr-levels 2))))))
(plot (getbb *result* :solution) :component 0)

;;; comparison with scalar case
(defparameter *result*
  (time
   (let* ((dim 2))
     (solve
      (blackboard
       :problem (cdr-model-problem dim)
       :output t :success-if '(> :nr-levels 2))))))
(plot (getbb *result* :solution) :component 0)

)

;;; (fl.application::test-elasticity-model-problem)
(fl.tests:adjoin-test 'test-elasticity-model-problem)
