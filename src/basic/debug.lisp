;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; debug.lisp - debug facility (from Peter Norvig's PAIP book)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;
;;; Copyright (C) 1998-2002 Peter Norvig.
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

;;; This utility is a modified version from code contained in
;;; http://www.norvig.com/paip/auxfns.lisp.

(in-package "COMMON-LISP-USER")

(defpackage "FL.DEBUG"
  (:use "COMMON-LISP")
  (:export "DBG-ON" "DBG-OFF" "DBG-WHEN" "DBG" "DBG-INDENT"))

(in-package :fl.debug)

(defvar *dbg-ids* () "Identifiers used by dbg.")

(defun dbg-on (&rest ids)
  "Register ids for dbg."
  (setf *dbg-ids* (union ids *dbg-ids*)))

(defun dbg-off (&rest ids)
  "Stop dbg on ids.  With no ids, stop dbg altogether."
  (setf *dbg-ids* (and ids (set-difference *dbg-ids* ids))))

(defmacro dbg-when (id &body body)
  "Perform a check only if debugging."
  `(when (member ,id *dbg-ids*)
    ,@body))

(defun dbg (id format-string &rest args)
  "Output of status information."
  (dbg-when id (format *debug-io* "~&~?" format-string args)
	    (force-output *debug-io*)))

(defun dbg-indent (id indent format-string &rest args)
  "Indented output of status information."
  (dbg-when id (format *debug-io* "~&~VT~?" indent format-string args)
	    (force-output *debug-io*)))


