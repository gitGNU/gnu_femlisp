;;; -*- mode: lisp; fill-column: 70; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; femlisp-matlisp.asd - System definition file
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;
;;; Copyright (C) 2017-
;;; Nicolas Neuss, Friedrich-Alexander-Universitaet Erlangen-Nuernberg
;;; All rights reserved.
;;; 
;;; Redistribution and use in source and binary forms, with or without
;;; modification, are permitted provided that the following conditions
;;; are met:
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
;;; MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
;;; IN NO EVENT SHALL THE AUTHOR, THE FAU ERLANGEN-NUERNBERG, OR OTHER
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

(asdf:defsystem :femlisp-matlisp
  :author "Nicolas Neuss"
  :license "Modified BSD"
  :depends-on (:femlisp-basic :femlisp-parallel :femlisp-dictionary)
  :pathname "../src"
  :around-compile call-with-read-double-float-environment
  :components
  ((:module
    "alien"
    :components
    (;; alien should be recompiled if superlu.so or umfpack.so has changed
     (:file "alien" :depends-on ())
     (:file "alienc" :depends-on ("alien"))
     (:file "lapack" :depends-on ("alien"))
     (:file "superlu" :depends-on ("alien"))
     (:file "umfpack" :depends-on ("alien"))))
   (:module
    "matlisp"
    :depends-on ("alien")
    :components
    ((:file "matlisp-defp")
     (:file "ctypes" :depends-on ("matlisp-defp"))
     (:file "vector" :depends-on ("matlisp-defp"))
     (:file "blas-basic" :depends-on ("vector"))
     (:file "matrix" :depends-on ("vector"))
     (:file "number-blas" :depends-on ("matrix" "blas-basic"))
     (:file "array-blas" :depends-on ("matrix" "ctypes" "number-blas"))
     (:file "store-vector" :depends-on ("matrix" "ctypes" "blas-basic"))
     (:file "standard-matrix" :depends-on ("matrix" "store-vector"))
     (:file "standard-matrix-blas" :depends-on ("standard-matrix"))
     (:file "standard-matrix-lr" :depends-on ("standard-matrix-blas"))
     (:file "compat" :depends-on ("standard-matrix"))
     (:file "call-matlisp" :depends-on ("standard-matrix"))
     (:file "tensor" :depends-on ("store-vector" "standard-matrix"))
     (:file "sparse-tensor" :depends-on ("tensor"))
     (:file "compressed" :depends-on ("store-vector" "standard-matrix"))
     (:file "ggev" :depends-on ("standard-matrix"))
     (:file "hegv" :depends-on ("standard-matrix"))
     (:file "sparse-vector" :depends-on ("compat" "sparse-tensor"))
     (:file "sparse-matrix" :depends-on ("sparse-vector" "standard-matrix-blas" "compressed"))
     (:file "sparselu" :depends-on ("sparse-matrix"))))))


