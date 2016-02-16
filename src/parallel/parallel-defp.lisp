;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; parallel-defp.lisp
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;
;;; Copyright (C) 2015 Nicolas Neuss, FAU Erlangen-Nuremberg.
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

(defpackage "FL.PARALLEL"
  (:use "COMMON-LISP" "FL.UTILITIES" "FL.MACROS" "FL.DEBUG"
        "LPARALLEL"  "LPARALLEL.KERNEL" "LPARALLEL.QUEUE"
        "BORDEAUX-THREADS")
  (:export
   "MUTEX-MIXIN" "WITH-MUTUAL-EXCLUSION"
   "WAITQUEUE-MIXIN" "MAKE-WAITQUEUE"
   "LOCKED-REGION-MIXIN" "WITH-REGION"
   
   "MP-DBG"
   "*KERNEL*" "WITH-WORKERS" "WORK-ON" "WITH-ACCUMULATORS"

   "PWORK"

   ;; multiprocessing.lisp
   "WITH-FEMLISP-WORKERS"
   "PARQUEUE"
;;   "WITH-PARALLEL-NETWORK" "PARCELL" "NAMED-PARCELL" "VALUE"
   )
  (:documentation
   "This package provides some Femlisp extensions for parallel execution
building on the BORDEAUX-THREADS and LPARALLEL libraries")
  )


