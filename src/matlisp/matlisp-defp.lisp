;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; matlisp-defp.lisp - Defines the Femlisp interface to Matlisp 
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

(in-package :cl-user)

(eval-when (:read-toplevel :compile-toplevel :execute)
  (defparameter *matlisp-symbols*
    '("STANDARD-MATRIX" "REAL-MATRIX"
      "NROWS" "NCOLS"
      "MAKE-FLOAT-MATRIX" "MAKE-REAL-MATRIX"
      "*PRINT-MATRIX*"
      "EYE" "ZEROS" "ONES"
      "SQUARE-MATRIX-P"
      "MREF" "MATRIX-REF"
      "JOIN"
      "COPY" "COPY!" "TRANSPOSE"
      "SCAL" "SCAL!" "M+" "M+!" "M-" "AXPY" "AXPY!"
      "GEMM" "GEMM!" "M*" "M*!"
      "GETRF!" "GETRS!" "GETRS" "GESV" "GESV!" "M/" "M/!"
      "M./" "M./!" "MAP-MATRIX"
      "DOT" "NORM"
      "HELP")))

(defpackage "FEMLISP.MATLISP"
  (:nicknames "FL.MATLISP")
  (:use "COMMON-LISP" "MACROS" "UTILITIES" "FEMLISP-DEBUG")
  #+matlisp
  (:import-from "MATLISP" . #.*matlisp-symbols*)
  (:export . #.*matlisp-symbols*))

