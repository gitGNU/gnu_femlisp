
;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; geomg-defp.lisp
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

(defpackage "FL.GEOMG"
  (:use "COMMON-LISP" "COMMON-LISP-USER" "FL.UTILITIES" "FL.MACROS" "FL.ALGEBRA" "FL.FUNCTION"
	"FL.MATLISP" "FL.ITERATION" "FL.MULTIGRID" "FL.MESH" "FL.DISCRETIZATION" "FL.DEBUG")
  (:export
   ;; geomg.lisp
   "<GEOMETRIC-CS>" "GEOMETRIC-CS" "<GEOMETRIC-FAS>" "FAS"
   "<S1-REDUCTION>" "<S1-COARSE-GRID-ITERATOR>" "S1-REDUCTION-AMG-SOLVER"
   ;; geoblock.lisp
   "GEOMETRIC-SSC" "GEOMETRIC-PSC"
   ;; vanka.lisp
   "<VANKA>")

  (:documentation "The @code{GEOMG} package contains iterations which
depend on geometric information, obtained for example from the
discretization.  At the moment, these are the geometric multigrid
iteration, an AMG-like scheme for preconditioning high-order
discretizations with low-order ones, and some block smoothers with
overlapping blocks."))
