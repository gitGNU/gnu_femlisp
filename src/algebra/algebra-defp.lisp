;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; algebra-defp.lisp
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

(defpackage "ALGEBRA"
  (:use "COMMON-LISP" "FL.MACROS" "FL.UTILITIES" "FL.MATLISP")
  (:export

   ;; tensor.lisp
   "<TENSOR>" "<REAL-TENSOR>" "<COMPLEX-TENSOR>" "MAKE-TENSOR"
   "MAKE-GENERAL-TENSOR" "MAKE-REAL-TENSOR" "LIST->REAL-TENSOR"
   "TENSOR-REF" "RANK" "*PRINT-TENSOR*"
   "SLICE" "COPY" "T+" "REARRANGE-TENSOR" "T*"
   "DOTENSOR" "TENSOR-FOR-EACH" "TENSOR-MAP"
   "K-JET" "EVALUATE-K-JET"
   
   ;; sparse-tensor.lisp
   "<SPARSE-TENSOR>" "IN-PATTERN-P" 
   
   ;; crs.lisp
   "CRS-PATTERN" "STORE-SIZE" "ROW-STARTS" "COL-INDS" "OFFSETS"
   "MAKE-CRS-PATTERN" "FULL-CRS-PATTERN" "PATTERN->FULL-PATTERN" "SHIFT-PATTERN" "CRS-MATRIX"
   "MAKE-FULL-CRS-MATRIX"

   ;; sparse.lisp
   "SHOW" "DISPLAY" "MAT-DIFF"
   "<SPARSE-VECTOR>" "BLOCKS" "KEY->SIZE" "PRINT-KEY" "MULTIPLICITY"
   "MAKE-SPARSE-VECTOR" "VECTOR-BLOCK" "KEYS"
   "<SPARSE-MATRIX>" "ROW-TABLE" "COLUMN-TABLE" "KEYS-OF-ROW" "KEYS-OF-COLUMN"
   "PRINT-ROW-KEY" "PRINT-COL-KEY" "ROW-KEY->SIZE" "COL-KEY->SIZE"
   "KEYS->PATTERN" "MAKE-SPARSE-MATRIX" "MAKE-SPARSE-AUTOMORPHISM" "MATRIX-ROW"
   "MATRIX-COLUMN" "MAKE-SPARSE-ANALOG" "MAKE-FULL-BLOCK-ANALOG"
   "SPARSE-VECTOR->MATLISP"
   "SET-SVEC-TO-LOCAL-BLOCK" "ADD-SVEC-TO-LOCAL-BLOCK"
   "MATRIX-BLOCK" "TOTAL-NROWS" "TOTAL-ENTRIES"
   "INDEX-RANGE-DISJOINT-P" "RANGE-AND-DOMAIN-DISJOINT-P"
   "REMOVE-ENTRY" "REMOVE-ROW" "REMOVE-COLUMN" "REMOVE-KEY" "REMOVE-KEYS" "ROW<-ID" "COLUMN<-ID"
   "KEYS->MBLOCKS"
   "PRINT-SMAT"
   "SYMMETRIC-P"
   "SPARSE-MATRIX->MATLISP"
   "EXTENDED-EXTRACT" "EXTRACT-IF" "EXTRACT-MATRIX-BLOCK"
   "COL-KEYS" "ROW-KEYS"
   "NR-NONEMPTY-ROWS" "NR-NONEMPTY-COLUMNS"
   "COMBINED-PROJECTION" "REMOVE-PROJECTION-RANGE" "EXTEND-BY-IDENTITY"

   ;; old sparse matrix representation
   ;;"<CONTAINER>"
   ;;"<VBLOCK>" "VBLOCK-TYPE" "CREATE-VBLOCK"
   ;;"<EXTENSIBLE-CONTAINER>" "TYPE->KEY" "KEY->TYPE"
   ;;"<CONTAINER-OBJECT>" "VBLOCK-SLICE" "DISPOSE" "KEYS->VBLOCKS"
   ;;"WITH-COBJS"
   ;;"<VSPACE>" "MAKE-EMPTY-VSPACE"
   ;;"<SVEC>" "VSPACE" "ALLOCATE-SVEC" "ALLOCATE-COPY" "VEC-SCALAR" "VECTOR-BLOCK"
   ;;"<GRAPH>" "MAKE-EMPTY-GRAPH" "ENSURE-MBLOCK"
   ;;"<SMAT>" "GRAPH" "ALLOCATE-SMAT" "MATRIX-BLOCK"

   ;; sparselu.lisp
   "SPARSE-LDU" "SPARSE-M*" "SHIFT-DIAGONAL-INVERTER"
   ))
