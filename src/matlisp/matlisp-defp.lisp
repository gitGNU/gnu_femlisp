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

(defpackage "FL.MATLISP"
  (:use "COMMON-LISP" "FL.MACROS" "FL.UTILITIES" "FL.DEBUG")
  (:export
   ;; Matlisp symbols
   "VLENGTH" "NROWS" "NCOLS"
   "STANDARD-MATRIX" "MAKE-REAL-MATRIX" "MAKE-REAL-VECTOR"
   "STORE"
   "*PRINT-MATRIX*"
   "EYE" "ZEROS" "ONES" "DIAG" "LAPLACE-FULL-MATRIX"
   "MREF"
   "FILL!" "FILL-RANDOM!"
   "JOIN" "TRANSPOSE" "TRANSPOSE!"
   "COPY" "COPY!"
   "SCAL" "SCAL!" "M+" "M+!" "M-" "M-!" "AXPY" "AXPY!"
   "GEMM!" "GEMM-NN!" "GEMM-NT!" "GEMM-TN!" "GEMM-TT!"
   "GEMM" "M*" "M*-TN" "M*-PRODUCT-INSTANCE"
   "GETRF!" "GETRF" "GETRS!" "GETRS"
   "GESV!" "GESV" "M/" "M/!"
   "M./" "M./!" "MAP-MATRIX"
   "DOT" "NORM" "L2-NORM" "LINF-NORM" "LP-NORM"
   
   ;; new Matlisp symbols

   ;; vector.lisp
   "<VECTOR>" "<STORE-VECTOR>"
   "VREF" "TOTAL-ENTRIES"

   ;; matrix.lisp
   "<MATRIX>"
   "FOR-EACH-ROW-KEY" "FOR-EACH-COL-KEY"
   "FOR-EACH-KEY-IN-ROW" "FOR-EACH-KEY-IN-COL"
   "FOR-EACH-ENTRY-IN-ROW" "FOR-EACH-ENTRY-IN-COL"
   "FOR-EACH-KEY-AND-ENTRY-IN-ROW" "FOR-EACH-KEY-AND-ENTRY-IN-COL"
   "NR-OF-ENTRIES" "DOROWS" "DOROW"
   "<SUBMATRIX>"
   "X<-AY" "X+=AY" "X-=AY" "X-ON-RANGE-OF-A<-0" "X-ON-DOMAIN-OF-A<-0"
   "MATRIX-SLICE" "VECTOR-SLICE" "CLEAR-ROW" "CLEAR-COLUMN"
   "MAKE-IMAGE-VECTOR-FOR" "MAKE-DOMAIN-VECTOR-FOR"
   "ELEMENT-TYPE" "SCALAR-TYPE"
   "MATRIX-TRANSPOSE-INSTANCE" "M*-PRODUCT-INSTANCE"
   "X<-0" "M-INCF"
   "MZEROP" "MEQUALP" "MIDENTITY-P"
   "FOR-EACH-KEY" "FOR-EACH-ENTRY" "FOR-EACH-KEY-AND-ENTRY"
   "FOR-EACH-ENTRY-OF-VEC1" "FOR-EACH-ENTRY-OF-VEC2"
   "DOVEC" "MEXTRACT" "MINJECT"
   
   ;; array-blas.lisp
   "DOUBLE-VEC" "LIST->DOUBLE-VEC" "MAKE-DOUBLE-VEC" "UNIT-VECTOR"
   
   ;; compat.lisp
   "AREA-OF-SPAN" "DET" "DET-FROM-LR"
   "DOT-ABS"   "SUBMATRIX"
   "ENSURE-MATLISP"
   "MULTIPLICITY"

   ))

