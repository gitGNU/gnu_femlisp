;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; tests.lisp - Femlisp test suite definitions
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

;;;; The idea of the Femlisp test suite is that every package or file is
;;;; able to add locally consistency checks.  Those tests can then be
;;;; performed before a CVS commit or a release.

(in-package "COMMON-LISP-USER")

(defpackage "TESTS"
  (:use "COMMON-LISP")
  (:export "ADJOIN-FEMLISP-TEST" "CLEAR-FEMLISP-TESTS" "TEST-FEMLISP"
	   "*FEMLISP-TEST-INTERNAL*"))

(in-package :tests)

(defparameter *femlisp-tests* ()
  "List of Femlisp tests.")

(defparameter *femlisp-bugs* ()
  "List of Femlisp bugs, i.e. failing tests that have not yet been fixed.")

(defun adjoin-femlisp-test (fsym)
  "Adjoins a test to the Femlisp test suite."
  (pushnew fsym *femlisp-tests*))
  
(defun clear-femlisp-tests ()
  "Clears the Femlisp test suite."
  (setq *femlisp-tests* ()))

(defun adjoin-femlisp-bug-test (fsym)
  "Register a bug test."
  (pushnew fsym *femlisp-bugs*))
  
(defun clear-femlisp-bug-tests ()
  "Clears the Femlisp bug register."
  (setq *femlisp-tests* ()))

(defun test-femlisp ()
  "Runs the Femlisp test suite.  The result is printed to
*standard-output*."
  (loop for fsym in (reverse *femlisp-tests*)
	for result =
	(catch 'trap
	  (handler-bind ((error #'(lambda (condition) (throw 'trap condition))))
	    (format t "~&~%***** Testing ~A *****~%~%" fsym)
	    (funcall fsym)
	    nil))
	collect (cons fsym result) into report
	do
	(format t "~A~%" fsym)
	(when result (format t "~A~%" result))
	finally
	(let ((failed (remove nil report :key #'cdr)))
	  (if failed
	      (loop initially (format t "~&~%The following tests failed:~%")
		    for (sym . error) in failed do
		    (format t "~%~%Test: ~A.~%Result:~A~%" sym error))
	      (format t "~%All tests passed.~%"))))
  
  ;; report also remaining bugs
  (when *femlisp-bugs*
    (format t "The following bugs are still open:~%")
    (loop for sym in *femlisp-bugs* do (format t "~A~%" sym)))
  )

;;; (tests::test-femlisp)
;;;
;;; Test without plotting:
;;; (let ((plot::*plot* nil)) (tests::test-femlisp))