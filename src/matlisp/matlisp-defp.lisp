;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; matlisp-defp.lisp - Defines the Femlisp interface to Matlisp 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;
;;; Copyright (C) 2003-2005 Nicolas Neuss, University of Heidelberg.
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
  (:use "COMMON-LISP" "FL.MACROS" "FL.UTILITIES"  "FL.DEBUG"
        "FL.PARALLEL" "FL.DICTIONARY"
        "FL.LAPACK")
  (:export
   ;; Matlisp symbols
   "VLENGTH" "NROWS" "NCOLS"
   "STANDARD-MATRIX" "STANDARD-MATRIX-P"
   "MAKE-REAL-MATRIX" "MAKE-REAL-VECTOR"
   "STORE"
   "*PRINT-MATRIX*" "*PRINT-MATRIX-ELEMENT-FORMAT*"
   "EYE" "ZEROS" "ONES" "DIAG" "DIAGONAL" "COLUMN" "ROW" "LAPLACE-FULL-MATRIX"
   "MREF"
   "FILL!" "FILL-RANDOM!"
   "JOIN" "JOIN-INSTANCE" "JOIN-HORIZONTAL!" "JOIN-VERTICAL!"
   "TRANSPOSE" "TRANSPOSE!"
   "COPY" "COPY!"
   "SCAL" "SCAL!" "M+" "M+!" "M-" "M-!" "M.*" "AXPY" "AXPY!"
   "GEMM!" "GEMM-NN!" "GEMM-NT!" "GEMM-TN!" "GEMM-TT!"
   "GEMM" "M*" "M*-TN" "M*-NT" "M*-PRODUCT-INSTANCE" "M*-TN-PRODUCT-INSTANCE"
   "GETRF!" "GETRF" "GETRS!" "GETRS"
   "GESV!" "GESV" "M/"
   "M./" "M./!" "MAP-MATRIX"
   "DOT" "NORM" "L2-NORM" "LINF-NORM" "LP-NORM"
   "NORMALIZE!" "NORMALIZE"
   
   ;; new Matlisp symbols

   ;; vector.lisp
   "<VECTOR>" "VREF" "NR-OF-ENTRIES" "TOTAL-ENTRIES"

   ;; store-vector.lisp
   "STORE-VECTOR"

   ;; matrix.lisp
   "<MATRIX>" "SHOW" "DISPLAY"
   "ACCESS-TYPE"
   "FOR-EACH-ROW-KEY" "FOR-EACH-COL-KEY"
   "FOR-EACH-KEY-IN-ROW" "FOR-EACH-KEY-IN-COL"
   "FOR-EACH-ENTRY-IN-ROW" "FOR-EACH-ENTRY-IN-COL"
   "FOR-EACH-KEY-AND-ENTRY-IN-ROW" "FOR-EACH-KEY-AND-ENTRY-IN-COL"
   "DOROWS" "DOROW"
   "<SUBMATRIX>"
   "X-ON-RANGE-OF-A<-0" "X-ON-DOMAIN-OF-A<-0"
   "MATRIX-SLICE" "VECTOR-SLICE"
   "MAKE-IMAGE-VECTOR-FOR" "MAKE-DOMAIN-VECTOR-FOR"
   "ELEMENT-TYPE" "SCALAR-TYPE"
   "MATRIX-TRANSPOSE-INSTANCE" "M*-PRODUCT-INSTANCE"
   "X<-0" "M-INCF"
   "MZEROP" "*MZEROP-THRESHOLD*" "MEQUALP" "MAT-DIFF" "MIDENTITY-P"
   "MSQUARE-P" "MSYMMETRIC-P"
   "FOR-EACH-KEY" "FOR-EACH-ENTRY" "FOR-EACH-ENTRY-AND-KEY"
   "FOR-EACH-ENTRY-AND-VECTOR-INDEX"
   "DOVEC" "MEXTRACT!" "MINJECT!"

   ;; standard-matrix.lisp
   "MRANDOM" "EXTEND-MATLISP-FUNCTION"

   ;; standard-matrix-lr.lisp
   "CHOLESKY"
   
   ;; ctypes.lisp
   "INT-VEC" "MAKE-INT-VEC"
   "UINT" "UINT-VEC" "MAKE-UINT-VEC"
   "DOUBLE-VEC" "MAKE-DOUBLE-VEC" "UNIT-VECTOR"
   
   ;; compat.lisp
   "AREA-OF-SPAN" "DET" "DET-FROM-LR"
   "DOT-ABS"   "SUBMATRIX"
   "ENSURE-MATLISP"
   "MULTIPLICITY"

   ;; compressed.lisp
   "COMPRESSED-PATTERN" "SIZES" "NUMBER-NONZERO-ENTRIES" "COMPRESSED-MATRIX"
   "FULL-COMPRESSED-PATTERN" "FULL-CCS-PATTERN" "FULL-CRS-PATTERN" "PATTERN"
   "MAKE-FULL-COMPRESED-MATRIX" "MAKE-FULL-CRS-MATRIX" "COMPRESSED->MATLISP"
   "TRANSPOSED-PATTERN"
   
   ;; tensor.lisp
   "TENSOR"  "IN-PATTERN-P" "ENTRY-ALLOWED-P" "FULL-TENSOR" "MAKE-REAL-TENSOR"
   "TENSOR-REF" "DIMENSIONS" "RANK" "*PRINT-TENSOR*"
   "SLICE" "T+" "REARRANGE-TENSOR" "T*"
   "DOTENSOR" "TENSOR-FOR-EACH" "TENSOR-MAP"

   ;; sparse-vector.lisp
   "SHOW" "DISPLAY" "MAT-DIFF"
   "<SPARSE-VECTOR>" "VECTOR-BLOCK" "BLOCKS" "KEY->SIZE" "<MULTIPLICITY-MIXIN>"
   "MAKE-SPARSE-VECTOR" "KEYS"
   
   ;; sparse-matrix.lisp
   "<SPARSE-MATRIX>" "MATRIX-BLOCK" "ROW-TABLE" "COLUMN-TABLE" "KEYS-OF-ROW" "KEYS-OF-COLUMN"
   "ROW-KEY->SIZE" "COL-KEY->SIZE"
   "KEYS->PATTERN" "MAKE-SPARSE-MATRIX" "MAKE-SPARSE-AUTOMORPHISM"
   "MATRIX-ROW-P" "MATRIX-COL-P" "MATRIX-ROW" "MATRIX-COLUMN"
   "MAKE-FULL-BLOCK-ANALOG"
   "SPARSE-VECTOR->MATLISP"
   "SET-SVEC-TO-LOCAL-BLOCK" "ADD-SVEC-TO-LOCAL-BLOCK"
   "VALUE-BLOCKS-IN-REGION"
   "MATRIX-BLOCK" "MAKE-MATRIX-BLOCK" "TOTAL-NROWS" "TOTAL-ENTRIES"
   "INDEX-RANGE-DISJOINT-P" "RANGE-AND-DOMAIN-DISJOINT-P"
   "REMOVE-ENTRY" "REMOVE-ROW" "REMOVE-COLUMN" "REMOVE-KEY" "REMOVE-KEYS" "ROW<-ID" "COLUMN<-ID"
   "KEYS->MBLOCKS"
   "PRINT-SMAT"
   "SYMMETRIC-P"
   "SPARSE-MATRIX->MATLISP"
   "EXTENDED-EXTRACT" "EXTRACT-IF" "EXTRACT-MATRIX-BLOCK" "EXTRACT-VALUE-BLOCKS"
   "COL-KEYS" "ROW-KEYS"
   "COMBINED-PROJECTION" "REMOVE-PROJECTION-RANGE" "EXTEND-BY-IDENTITY"
   "LAPLACE-SPARSE-MATRIX"
   
   ;; sparselu.lisp
   "SPARSE-LDU" "SPARSE-M*" "SHIFT-DIAGONAL-INVERTER"
   "SPARSE-MATRIX->CCS"
   
   ;; sparse-tensor.lisp
   "<SPARSE-TENSOR>" "SPARSE-TENSOR" "DIAGONAL-SPARSE-TENSOR" "INDICES" "ENTRIES"
   
   ;; LAPACK routines
   "HEGV" "GGEV"
   )
  (:documentation "This package provides a Common Lisp version of full
  matrices with elements being numbers of a given type.  Those classes are
automatically generated when needed.  It provides also part of the BLAS and
LAPACK operations for those matrices.  The corresponding methods are
automatically compiled for the given matrix classes.  The interface is very
similar to the library Matlisp @cite{(Matlisp)}, which provides a CLOS
interface to the Fortran BLAS and LAPACK routines."))
