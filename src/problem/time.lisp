;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; time.lisp - Time-dependent problems
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

(in-package :problem)

(defclass <time-dependent-problem> (<problem>)
  ((stationary-problem :accessor stationary-problem
    :initarg :stationary-problem :documentation
    "Stationary problem F(u)=A(u)-f=0 for which this problem is the
time-dependent version."))
  (:documentation "The problem d/dt(alpha u) + F(u) = 0, where F(u) is
given by stationary-problem."))

(defmethod interior-coefficients ((problem <time-dependent-problem>))
  "Derived interior coefficients for a time-dependent problem."
  (append '(INITIAL ALPHA)
	  (interior-coefficients (slot-value problem 'stationary-problem))))

(defmethod boundary-coefficients ((problem <time-dependent-problem>))
  "Boundary coefficients for a time-dependent problem."
  (boundary-coefficients (slot-value problem 'stationary-problem)))

(defmethod self-adjoint-p ((problem <time-dependent-problem>))
  nil)

(defmethod initialize-instance :after ((problem <time-dependent-problem>) &key &allow-other-keys)
  (setf (slot-value problem 'multiplicity)
	(multiplicity (slot-value problem 'stationary-problem))))

