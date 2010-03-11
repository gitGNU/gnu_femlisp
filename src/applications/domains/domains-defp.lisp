;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; porous-domains.lisp - Generating periodic porous domains
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;
;;; Copyright (C) 2005 Nicolas Neuss, University of Heidelberg.
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

(defpackage "FL.DOMAINS"
  (:use "COMMON-LISP"
	"FL.MACROS" "FL.UTILITIES" "FL.MATLISP"
	"FL.DEBUG" "FL.DEMO"
	"FL.FUNCTION" "FL.MESH")
  (:export ; inlay-domain.lisp
   "N-CUBE-WITH-CUBIC-INLAY" "N-CELL-WITH-CUBIC-INLAY"
   "N-CUBE-WITH-BALL-INLAY" "N-CELL-WITH-BALL-INLAY"
   "PATCH-IN-INLAY-P" "PATCH-ON-N-CUBE-BOUNDARY-P")
  (:export ; hole-domain.lisp
   "N-CUBE-WITH-CUBIC-HOLE" "N-CELL-WITH-CUBIC-HOLE"
   "N-CUBE-WITH-ELLIPSOIDAL-HOLE" "N-CELL-WITH-ELLIPSOIDAL-HOLE"
   "N-CUBE-WITH-BALL-HOLE" "N-CELL-WITH-BALL-HOLE"
   "PATCH-ON-INNER-BOUNDARY-P")
  (:export ; bl-cell.lisp
   "OSCILLATING-BOUNDARY-DOMAIN" "SINUSOIDAL-BL-CELL"
   "SPLINE-INTERPOLATED-BL-CELL"
   "BL-PATCH-ON-LOWER-BOUNDARY" "BL-PATCH-ON-PELLET-BOUNDARY"
   "BL-PATCH-ON-UPPER-BOUNDARY" "BL-PATCH-ON-ARTIFICIAL-BOUNDARY")
  (:documentation "Femlisp package for domain definitions."))
