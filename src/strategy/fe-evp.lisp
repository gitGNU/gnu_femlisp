;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; fe-evp.lisp - Solution strategy for PDE eigenvalue problems
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;
;;; Copyright (C) 2004 Nicolas Neuss, University of Heidelberg.
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

(in-package :fl.strategy)

(file-documentation
 "This file contains solution strategies for pde eigenvalue problems.")

(defclass <fe-evp-strategy> (<stationary-fe-strategy>)
  ()
  (:documentation "This class is a modification of
@class{<stationary-fe-strategy>} for solving PDE eigenvalue problems by a
Wielandt iteration."))

#+nil
(defmethod mass ((evp <evp>) (solution <ansatz-space-vector>))
  (let ((bb (blackboard :ansatz-space (ansatz-space solution) :solution solution
			:mass-factor 1.0 :stiffness-factor 0.0)))
    (fe-discretize bb)
    (let ((lse (linearize (getbb bb :discretized-problem) solution)))
      (dot solution (m* (matrix lse) solution)))))

#+nil
(defmethod energy ((evp <evp>) (solution <ansatz-space-vector>))
  (let ((bb (blackboard :ansatz-space (ansatz-space solution) :solution solution
			:mass-factor 0.0 :stiffness-factor 1.0)))
    (fe-discretize bb)
    (let ((lse (linearize (getbb bb :discretized-problem) solution)))
      (dot solution (m* (matrix lse) solution)))))

(defmethod approximate :after ((fe-strategy <fe-evp-strategy>) blackboard)
  "Ensures accuracy of the solution and the error estimate."
  (with-items (&key discretized-problem solver-blackboard
		    solution residual output mass-factor stiffness-factor)
      blackboard
    ;; assemble (better would be a local assembly)
    (fe-discretize blackboard)
    ;; improve approximation by solving
    (dbg :iter "Approximate: lambda=~A mu=~A"
	 (unbox (slot-value discretized-problem 'lambda))
	 (unbox (slot-value discretized-problem 'mu)))
    (setf solver-blackboard (blackboard :problem discretized-problem
					:solution solution :residual residual
					:output output))
    (solve (slot-value fe-strategy 'solver) solver-blackboard)
    (setf solution (getbb solver-blackboard :solution))))

