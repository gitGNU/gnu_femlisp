;;; -*- Mode: LISP; Syntax: COMMON-LISP; Package: CL-USER; Base: 10 -*-

;;; Copyright (C) 2016, Dr. Nicolas Neuss.  All rights reserved.

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

(asdf:defsystem :femlisp-ddo
  :author "Nicolas Neuss"
  :license "Modified BSD"
  :components (
	       (:file "packages" :depends-on ())
	       (:file "mesh-ddo" :depends-on ("packages"))
	       (:file "sparseas-ddo" :depends-on ("packages"))
	       (:file "discretize-ddo" :depends-on ("sparseas-ddo"))
	       (:file "solve-ddo" :depends-on ("sparseas-ddo"))
	       (:file "strategy-ddo" :depends-on ("packages"))
	       (:file "hom-ddo" :depends-on ("packages"))
	       (:file "elahom-testing" :depends-on ("hom-ddo" "solve-ddo" "strategy-ddo"))
	       (:file "heisig-neuss-2017-2" :depends-on ("elahom-testing"))
               )
  :depends-on (:femlisp
               :lfarm-server :lfarm-admin :lfarm-client
               :cl-mpi :cl-mpi-extensions
               :uiop :trees :alexandria
               :net.scipolis.graphs
               :ddo))

(asdf:defsystem :femlisp-mpi-worker
  :author "Nicolas Neuss"
  :license "Modified BSD"
  :components ()  ; no components here, because :mpi-worker saves a core
  :depends-on (:femlisp :ddo :femlisp-ddo :mpi-worker))
