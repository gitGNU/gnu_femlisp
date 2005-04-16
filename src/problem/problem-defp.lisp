;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; problem-defp.lisp
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

(defpackage "FL.PROBLEM"
  (:use "COMMON-LISP" "FL.MACROS" "FL.UTILITIES" "FL.MATLISP"
	"FL.MESH" "FL.ALGEBRA" "FL.FUNCTION")
  (:export				; problem.lisp
   "<PROBLEM>" "GET-PROPERTY" "LINEAR-P" "<NONLINEAR-PROBLEM>"
   "<SOLVER>" "SOLVE" "SELECT-SOLVER" "SELECT-LINEAR-SOLVER"
   "LINEAR-P" "COERCIVE"
   "<LSE>" "LSE" "MATRIX" "RHS" "ENSURE-SOLUTION" "ENSURE-RESIDUAL"
   "<NLSE>" "NLSE" "LINEARIZE")
  (:export				; evp.lisp
   "<EVP-MIXIN>" "LAMBDA" "MU" "<EVP>" "<LS-EVP>" "MASS" "ENERGY")
  (:export				; pde-problem.lisp
   "<DOMAIN-PROBLEM>" "<INTERPOLATION-PROBLEM>" "<PDE-PROBLEM>"
   "DOMAIN-DIMENSION" "NR-OF-COMPONENTS"
   "PATCH->COEFFICIENTS" "MULTIPLICITY"
   "SELF-ADJOINT-P" "DUAL-PROBLEM"
   "COEFFICIENTS-OF-PATCH" "COEFFICIENTS-OF-CELL"
   "COEFFICIENTS" "INTERIOR-COEFFICIENTS" "BOUNDARY-COEFFICIENTS"
   "INTERIOR-COEFFICIENT-P" "BOUNDARY-COEFFICIENT-P" "COEFFICIENT-P"
   "CONSTRAINT-IDENTIFIER"
   "<COEFFICIENT>" "DEMANDS" "FUNCTION->COEFFICIENT" "ENSURE-COEFFICIENT"
   "CONSTANT-COEFFICIENT" "CONSTRAINT-COEFFICIENT" "ZERO-CONSTRAINTS"
   "*EMPTY-COEFFICIENT-INPUT*")
  (:export				; time.lisp
   "<TIME-DEPENDENT-PROBLEM>" "STATIONARY-PROBLEM" "STATIONARY-PROBLEM-CLASS"
   "INITIAL" "ALPHA")
  (:documentation " The @package{FL.PROBLEM} package introduces the
general class @class{<problem>} and some subclasses.  The most interesting
subclass is @class{<pde-problem>}.  A @class{<pde-problem>} is defined on a
domain and provides a table mapping the domain patches to property lists
containing the coefficient functions.

Several subclasses of @class{<pde-problem>} are defined in own packages,
e.g. @class{<cdr-problem>} in @package{FL.CDR},
@class{<elasticity-problem>} in @package{FL.ELASTICITY} and
@class{<navier-stokes-problem} in @package{FL.NAVIER-STOKES}."))


