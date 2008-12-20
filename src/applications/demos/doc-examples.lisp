;;; -*- mode: lisp; fill-column: 64; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; doc-examples.lisp
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;
;;; Copyright (C) 2007 Nicolas Neuss, University of Karlsruhe.
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

;;; Elasticity
(let* ((problem
        (elasticity-model-problem
         (n-cube-domain 2) :lambda 1.0 :mu 1.0
         :force (vector #m(1.0) #m(0.0))))
       (blackboard
        (blackboard :problem problem :output t
                    :success-if '(> :time 20.0))))
  (storing (solve blackboard)))

(plot (getbb *result* :solution) :component 'fl.elasticity::u :rank 1 :shape 2)

;;; Navier-Stokes
(let ((dim 2))
  (storing
    (solve 
     (blackboard
      :problem (driven-cavity dim :reynolds 10.0) :base-level 0
      :success-if '(> :time 20.0) :output t :observe
      (append *stationary-fe-strategy-observe*
              (list
               (list (format nil "~{                 u~1D~}" '(1 2))
                     "~{~19,10,2E~}"
                     #'(lambda (blackboard)
                         (let ((val (fe-value (getbb blackboard :solution)
                                              #d(0.5 0.5))))
                           (loop for i below dim collect
                                (vref (aref val i) 0)))))))))))

;;; time-dependent
(let* ((dim 1) (levels 4) (order 2)
       (problem (cdr-model-problem
                 dim :initial #'(lambda (x) #I(sin(2*pi*x[0]^^2)))
                 :reaction 0.0 :source #m(0.0)))
       (rothe (make-instance
               '<rothe> :model-time 0.0 :time-step 0.01
               :stationary-success-if `(> :nr-levels ,levels)
               :success-if '(>= :step 20)
               :output t :plot t)))
  (storing
    (iterate rothe (blackboard
                    :problem problem :fe-class (lagrange-fe order)
                    :plot-mesh nil :output t))))

;;; eigenvalues
(let ((problem (cdr-model-problem 1 :evp (list :lambda (box 10.0)
                                               :mu (box 1.0)))))
  (storing
    (solve (blackboard :problem problem
                       :success-if '(or (>= :time 5) (>= :nr-levels 5))
                       :output :all))))
(slot-value (getbb *result* :problem) 'lambda)
(plot (getbb *result* :solution))