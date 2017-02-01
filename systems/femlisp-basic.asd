;;; -*- mode: lisp; fill-column: 70; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; femlisp-basic.asd - System definition file
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

(asdf:defsystem :femlisp-basic
  :author "Nicolas Neuss"
  :license "Modified BSD"
  :depends-on (#+(or clisp ccl) :cffi
               #+sbcl :sb-posix #+sbcl :sb-introspect
               #+allegro (:require "osi")
               :closer-mop :fiveam
               )
  :pathname "../src"
  :components
  ((:module
    "config"
    :depends-on ()
    :components
    ((:file "setup")
     (:file "femlisp-config" :depends-on ("setup"))))
   (:module
    "basic"
    :depends-on ("config")
    :components
    ((:file "debug" :depends-on ())
     (:file "tests" :depends-on ())
     (:file "patches" :depends-on ())
     (:file "macros" :depends-on ())
     (:file "port" :depends-on ("debug"))
     (:file "utilities-defp" :depends-on ("patches" "macros" "debug"))
     (:file "utilities" :depends-on ("utilities-defp" "macros" "tests"))
     (:file "amop" :depends-on ("debug" "port" "utilities"))
     (:file "mflop" :depends-on ("debug" "utilities"))
     (:file "general" :depends-on ("amop" "utilities-defp"))
     (:file "demo" :depends-on ("tests" "mflop" "macros" "utilities"))))))
