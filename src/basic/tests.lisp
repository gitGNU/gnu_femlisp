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

(defpackage "FL.TESTS"
  (:use "COMMON-LISP")
  (:export "ADJOIN-TEST" "REMOVE-TEST" "TEST-FEMLISP")
  (:documentation "This package provides routines for building a simple
regression test suite.  Most files in @femlisp{} contain a test function at
the end which checks several critical features which the file or module
provides.  By calling the function @function{adjoin-test} at load time,
this function is added to a list of functions to be checked.  After loading
@femlisp{}, all functions in this list can be executed one after the other
by calling the function @function{test-femlisp}.  Errors and exceptions are
registered and finally reported.  It is very much recommended to run this
test suite before a release."))

(in-package :fl.tests)

(defvar *tests* ()
  "List of Femlisp tests.")

(defvar *bugs* ()
  "List of Femlisp bugs, i.e. failing tests that have not yet been fixed.")

(defvar *testing* ()
  "List of routines which still have to be tested.")

(defvar *failed* ()
  "List of routines which failed for last run of the test suite.")

(defun adjoin-test (fsym)
  "Adjoins a test to the Femlisp test suite."
  (pushnew fsym *tests*))
  
(defun remove-test (fsym)
  "Adjoins a test to the Femlisp test suite."
  (setq *tests* (remove fsym *tests*)))
  
(defun clear-tests ()
  "Clears the Femlisp test suite."
  (setq *tests* ()))

(defun adjoin-bug-test (fsym)
  "Register a bug test."
  (pushnew fsym *bugs*))
  
(defun clear-bug-tests ()
  "Clears the Femlisp bug register."
  (setq *tests* ()))

(defun test-femlisp (&key continue package (logfile #p"femlisp:fltest.log"))
  "Runs the Femlisp test suite.  The result is printed to
*standard-output*."
  (unless continue
    (setq *testing* (if package
			(remove (find-package package) *tests*
				:key #'symbol-package :test-not #'eq)
			*tests*))
    (setq *failed* nil))
  (flet ((run-tests ()
	   (loop for fsym = (pop *testing*) while fsym
	      for result =
		#+gcl (funcall fsym)
		#-gcl
		(catch 'trap
		  (handler-bind ((error #'(lambda (condition) (throw 'trap condition))))
		    (format t "~&~%***** Testing ~A *****~%~%" fsym)
		    (funcall fsym)
		    nil))
	      when result do (push (cons fsym result) *failed*)
	      finally (return *failed*))))
    (let ((output-string
	   (with-output-to-string (stream)
	     (multiple-value-bind (secs mins hrs day month year)
		 (get-decoded-time)
	       (format stream "~&Started tests at ~D.~D.~D at ~D:~2,'0D:~2,'0D~%"
		       day month year hrs mins secs))
	     (run-tests)
	     (multiple-value-bind (secs mins hrs day month year)
		 (get-decoded-time)
	       (format stream "~&Finished tests at ~D.~D.~D at ~D:~2,'0D:~2,'0D~%"
		       day month year hrs mins secs))
	     (if *failed*
		 (loop initially (format stream "~&~%The following tests failed:~%")
		    for (sym . error) in *failed* do
		      (format stream "~%~%Test: ~A.~%Result:~A~%" sym error))
		 (format stream "~%All tests passed.~%"))
	     (when *bugs*
	       (format stream "The following bugs are still open:~%")
	       (loop for sym in *bugs* do (format stream "~A~%" sym))))))
      (write-string output-string)
      (when logfile
	(with-open-file (stream logfile :direction :output :if-exists :append)
	  (write-string output-string stream))))))


;;; (time (fl.tests:test-femlisp :package :fl.mesh))
;;;
;;; Test without plotting:
;;; (time (let ((fl.plot::*plot* nil)) (fl.tests:test-femlisp)))