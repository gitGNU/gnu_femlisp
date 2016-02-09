;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; dictionary-defp.lisp
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;
;;; Copyright (C) 2003-2006 Nicolas Neuss, University of Heidelberg.
;;; Copyright (C) 2006-2011 Nicolas Neuss, KIT Karlsruhe.
;;; Copyright (C) 2011-
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
;;; IN NO EVENT SHALL THE AUTHOR, THE KARLSRUHE INSTITUTE OF TECHNOLOGY,
;;; OR OTHER CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
;;; SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
;;; LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
;;; DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
;;; THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
;;; (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
;;; OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defpackage "FL.DICTIONARY"
  (:use "COMMON-LISP" "FL.MACROS" "FL.DEBUG" "FL.UTILITIES" "FL.PARALLEL")
  (:export  ; dictionary.lisp
   "MAKE-ANALOG"
   "DIC-REF" "DIC-FOR-EACH" "DIC-FOR-EACH-KEY" "DIC-FOR-EACH-VALUE"
   "DIC-REMOVE" "DIC-POP" "DIC-EMPTY-P" "KEYS" "DODIC"
   "SMALL-CACHE-DICTIONARY" "SORTED-HASH-TABLE" "CACHE-DICTIONARY"
   "COMPUTED-VALUE-DICTIONARY"
   "WITH-MEMOIZATION" "MEMOIZING" "MEMOIZING-LET"  "MEMOIZING-LET*"
   ;; parallel-heap.lisp
   "EXTRACT-COMPLETE-SUBGRAPH" "INVERT-GRAPH"
   "MAKE-PARALLEL-HEAP" "TAKE-OBJECT" "DROP-OBJECT"
   "PROCESS-IN-PARALLEL"
   )
  (:documentation
   "This package contains functions for dictionaries as a general concept for
association lists, hash tables, etc."))
