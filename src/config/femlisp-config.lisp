;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; femlisp-config.lisp - configuration file for femlisp
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

(in-package :fl.start)

;;; Set the following variables only if DX or Gnuplot should not be found
;;; automatically or if you are not satisfied with the standard directory
;;; where images or meshes are written to.

(defparameter *dx-path*
  nil ; default: search in system path
  "Path to the @program{DX} executable.")

(defparameter *gnuplot-path*
  nil ; default: search in system path
  "Path to the @program{Gnuplot} executable.")

(defparameter *triangle-path*
  nil ; default: search in system path
  "Path to the @program{triangle} executable.")
  
(defparameter *tetgen-path*
  nil ; default: search in system path
  "Path to the @program{tetgen} executable.")
  
(defparameter *images-directory*
  nil ; default: (femlisp-pathname "images/")
  "Directory where images are put by default.")

(defparameter *meshes-directory*
  nil ; default: (femlisp-pathname "meshes/")
  "Directory where meshes are put by default.")

(defparameter *blas-library*
  (or
   ;; :none
   (probe-file #p"femlisp:interface;mkl.so")
   #+linux "libblas.so"
   #+os-macosx "libblas.dylib")
  "BLAS library path.
If the value is NIL, the BLAS library is searched for in standard library paths.
If the value is :none, no external BLAS routines are used.")

(defparameter *lapack-library*
  ;; default: search in system libraries, e.g. #p"/usr/lib/liblapack.so"
  (or
   ;;:none
   ;;(probe-file #p"femlisp:interface;mkl.so")
   #+linux "liblapack.so"
   #+os-macosx "liblapack.dylib"
   )
  "LAPACK library path.
If the value is NIL, the LAPACK library is searched for in standard library paths.
If the value is :none, no external LAPACK routines are used.")

(defparameter *superlu-library*
  (probe-file (femlisp-pathname "interface/superlu.so"))
  "Wrapper for SuperLU, if available.")

(defparameter *umfpack-library*
  nil ; default: (femlisp-pathname "interface/umfpack.so")
  "Wrapper for UMFPACK, if available.")

