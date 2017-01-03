;;; -*- Mode: LISP; Syntax: COMMON-LISP; Package: CL-USER; Base: 10 -*-

;;; Copyright (C) 2015, Dr. Nicolas Neuss.  All rights reserved.

;;; Redistribution and use in source and binary forms, with or without
;;; modification, are permitted provided that the following conditions
;;; are met:

;;;   * Redistributions of source code must retain the above copyright
;;;     notice, this list of conditions and the following disclaimer.

;;;   * Redistributions in binary form must reproduce the above
;;;     copyright notice, this list of conditions and the following
;;;     disclaimer in the documentation and/or other materials
;;;     provided with the distribution.

;;; THIS SOFTWARE IS PROVIDED BY THE AUTHOR 'AS IS' AND ANY EXPRESSED
;;; OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
;;; WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
;;; ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY
;;; DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
;;; DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
;;; GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
;;; INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
;;; WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
;;; NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
;;; SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

(in-package :cl-user)

(asdf:defsystem :ddo-femlisp
  :serial t
  :version "0.2.0"
  :components (
	       (:file "packages" :depends-on ())
	       (:file "ddo-test" :depends-on ("packages"))
	       (:file "mesh-ddo" :depends-on ("packages"))
	       (:file "sparseas-ddo" :depends-on ("packages"))
	       (:file "discretize-ddo" :depends-on ("ddo-sparseas"))
	       (:file "solve-ddo" :depends-on ("ddo-sparseas"))
	       (:file "strategy-ddo" :depends-on ("packages"))
	       (:file "hom-ddo" :depends-on ("packages"))
	       (:file "elahom-testing" :depends-on ("ddo-hom" "ddo-solve" "ddo-strategy"))
               )
  :depends-on (:lfarm-server
               :lfarm-admin :lfarm-client
               :cl-mpi :cl-mpi-extensions
               :uiop :trees :alexandria
               :ddo :femlisp))

(asdf:defsystem :femlisp-initialize
  :serial t
  :version "0.2.0"
  :components (
	       (:file "elahom-initialize")
               )
  :depends-on (:femlisp-niko))

(asdf:defsystem :femlisp-save-core
  :serial t
  :version "0.2.0"
  :components (
	       (:file "save-core")
               )
  :depends-on (:femlisp-initialize))

(asdf:defsystem :ddo-femlisp-initialize
  :serial t
  :version "0.2.0"
  :components (
	       (:file "elahom-initialize")
               )
  :depends-on (:ddo :femlisp-initialize))

(asdf:defsystem :ddo-femlisp-save-core
  :serial t
  :version "0.2.0"
  :components (
	       (:file "save-core")
               )
  :depends-on (:ddo-femlisp-initialize))
