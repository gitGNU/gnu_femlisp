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

(defpackage "FL.ALGEBRA"
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
   "MAKE-SPARSE-VECTOR" "KEYS"
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
   "LAPLACE-SPARSE-MATRIX"

   ;; sparselu.lisp
   "SPARSE-LDU" "SPARSE-M*" "SHIFT-DIAGONAL-INVERTER"
   "SPARSE-MATRIX->CCS"
   )
  (:documentation "This package defines classes for sparse matrices and
methods operating on them.  The interface is mostly the one used in the
package @package{FL.MATLISP} extended suitably.

This module contains definitions for doing linear algebra in @femlisp{} and
consists of several files.  Besides others, these are @path{vector.lisp}
and @path{matrix.lisp}, where an abstract interface for linear algebra on
vectors and matrices is defined.  In @path{matlisp.lisp}, this interface is
realised for Matlisp matrices, and in addition the Matlisp operations are
partially extended to arrays.  In @path{crs.lisp}, the well-known compact
row-ordered storage is defined for storing sparse matrices, and
@path{tensor.lisp} and @path{sparse-tensor.lisp} contain class and method
definitions for full and sparse tensors of arbitrary rank.

The file @path{sparse.lisp} then introduces a sparse storage scheme for
block vectors and block matrices over arbitrary index sets.  This is very
convenient for functions and linear operators defined on unstructured grids
because the geometric grid objects themselves can index their degrees of
freedom.  It also allows for local updates when advanced adaptive schemes
like the one proposed in @cite{URuede_1993b} are used.

At the moment, the classes @class{<sparse-matrix>} and
@class{<sparse-matrix>} which are defined in @file{sparse.lisp} are used
almost exclusively.  Those classes are hash-table based, which implies that
almost every index set can be used.  Basic methods for those classes are
also defined in @path{sparse.lisp}.  An LU decomposition for
@class{<sparse-matrix>} is implemented in @file{sparselu.lisp}.

Unfortunately, the use of hash-tables instead of arrays is much slower than
working with, for example, compact row-ordered storage.  Thus, good
performance can only be expected if the inner blocks are relatively large.
This is the case for systems of equations and/or approximations of higher
order.  Future versions of @femlisp{} will probably improve on that by
using an array-based scheme similar to the one in
@path{sparse-tensor.lisp}."))
