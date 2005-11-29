;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; strategy-defp.lisp
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

(defpackage "FL.STRATEGY"
  (:use "COMMON-LISP" "FL.MACROS" "FL.UTILITIES" "FL.DEBUG"
	"FL.MATLISP" "FL.ALGEBRA" "FL.FUNCTION" "FL.MESH"
	"FL.PROBLEM" "FL.DISCRETIZATION"
	"FL.CDR" "FL.CDRSYS"
	"FL.CDR-FE" "FL.ITERATION" "FL.GEOMG"
	"FL.PLOT")
  (:export				; stationary.lisp
   "<STRATEGY>" "<FE-APPROXIMATION>"
   "<FE-INTERPOLATION>" "<STATIONARY-FE-STRATEGY>"
   "*STATIONARY-FE-STRATEGY-OBSERVE*" "*ETA-OBSERVE*")
  (:export				; error-estimator.lisp
   "ESTIMATE" "<PROJECTION-ERROR-ESTIMATOR>" "<DUALITY-ERROR-ESTIMATOR>")
  (:export				; error-indicator.lisp
   "<REFINEMENT-INDICATOR>" "INDICATE" "<UNIFORM-REFINEMENT-INDICATOR>"
   "<REGION-INDICATOR>"
   "<LARGEST-ETA-INDICATOR>" "<ABOVE-THRESHOLD-INDICATOR>")
  (:export   ; rothe.lisp
   "<ROTHE>")
  (:documentation "This package provides methods for solving problems by
adaptive FEM."))
