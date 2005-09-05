;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; superlu.lisp - interface to SuperLU
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;
;;; Copyright (C) 2004 Nicolas Neuss, University of Heidelberg.
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

(in-package "FL.ALIEN")

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; From superlu/util.h
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; SuperLU options
;; (enumerate yes_no_t NO YES)
;; (enumerate fact_t DOFACT SamePattern SamePattern_SameRowPerm FACTORED)
;; (enumerate rowperm_t NOROWPERM LargeDiag MY_PERMR)
;; (enumerate colperm_t NATURAL MMD_ATA MMD_AT_PLUS_A COLAMD MY_PERMC)
;; (enumerate trans_t NOTRANS TRANS CONJ)
;; (enumerate DiagScale_t NOEQUIL ROW COL BOTH)
;; (enumerate IterRefine_t NOREFINE SINGLE=1 DOUBLE EXTRA)
;; (enumerate MemType LUSUP UCOL LSUB USUB)
;; (enumerate stack_end_t HEAD TAIL)
;; (enumerate LU_space_t SYSTEM USER)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; SuperLU interface function
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

#+superlu
(fl.port:def-function ("c_superlu" c-superlu)
    ((m :int) (n :int) (nnz :int)
     (colptr (* :int))
     (rowind (* :int))
     (nzval (* :double))
     (nrhs :int)
     (rhs (* :double))
     (sol (* :double)))
  :returning :int)

#+superlu
(defun superlu (m n nnz cs ri store nrhs rhs sol)
  "Calls SuperLU."
  (fl.port:foreign-call-wrapper
   (c-superlu m n nnz (vector-sap cs) (vector-sap ri)
	      (vector-sap store) nrhs (vector-sap rhs)
	      (vector-sap sol))))

;;; Testing

(defun test-superlu ()
  #+superlu (direct-solver-test 'superlu)
  )

(fl.tests:adjoin-test 'test-superlu)

