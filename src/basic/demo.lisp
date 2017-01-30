;;; -*- mode: lisp; fill-column: 64; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; demo.lisp - Femlisp demo suite
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

(in-package "COMMON-LISP-USER")

(defpackage "FL.DEMO"
  (:use "COMMON-LISP" "FL.MACROS" "FL.UTILITIES")
  (:export "FEMLISP-DEMO" "LEAVES" "MAKE-DEMO" "ADJOIN-DEMO"
	   "REMOVE-DEMO" "FIND-DEMO"
	   "*DEMO-ROOT*" "*DEMO-TIME*"
	   "USER-INPUT" "USER-INPUT-TEXTFIELD"
	   "EXTRACT-DEMO-STRINGS"
	   "TEST-ALL-DEMOS")
  (:documentation "This package provides routines for building a
demo suite.  Wherever something interesting can be demonstrated,
a small demo node should be generated with the
@function{make-demo} and added to the tree of all demos with
@function{adjoin-demo}.  After loading @femlisp{}, the whole
demo suite is available and can be run with the command
@function{femlisp-demo}."))

(in-package :fl.demo)

(defclass <demo> ()
  ((name :accessor name :initarg :name :type string
	 :documentation "The name of the demo.")
   (short :initarg :short :initform nil
	  :documentation "A short description.")
   (long :initarg :long :initform nil
	 :documentation "A long description.")
   (leaves :accessor leaves :initform () :type list
	   :documentation "The child nodes of the demo.")
   (execute :accessor execute :initarg :execute :initform nil
	    :documentation "A function performing the demo.")
   (input :initarg :test-input :initform nil
	  :documentation "Sample string of user input for testing."))
  (:documentation "Femlisp demo node."))

(defun make-demo (&rest initargs)
  (apply #'make-instance '<demo> initargs))

(defvar *demo-help*
  "?, help    Displays this help.
up, back   Goes up in the demo hierarchy.
quit       Leave the demo program.

Other words you enter are matched with leave names.
Abbreviations are possible, e.g. 'hom-2' for
'homogenization-2d'.  Note that case is ignored.")

(defvar *demo-root*
  (make-demo
   :name "Demo"
   :short "Femlisp demo suite"
   :long
   "This is the root of the standard Femlisp demos.  Choose the
demo you want to see by typing its name or an abbreviation.
Such abbreviations can be of the form `h-2' for
`homogenization-2D'.

Type `up' or `back' to go up in the demo hierarchy, `quit' to
quit, and '?' or 'help' for help."))

(defvar *visited-demos* (make-hash-table :test 'equal)
  "Table of demos already visited during this Lisp session.")

(defvar *demo-time* 10.0
  "Suggested time for a demo which should be used as a
termination criterion in demo iterations.")

(defun update-demo-status (demo)
  "Checks if all leaves have been visited.  If yes, the demo is
added to the visited demos."
  (when (every (lambda (leaf) (gethash (name leaf) *visited-demos*)) (leaves demo))
    (setf (gethash (name demo) *visited-demos*) t)))

(defun fuzzy-match-strings (string1 string2)
  "Match two strings allowing abbreviations of the form `h-2'
for `homogenization-2d'."
  (declare (optimize debug))
  (loop for pos1 below (length string1)
	and pos2 below (length string2) do
	;;until (or (= pos1 (length string1)) (= pos2 (length string2)))
	(when (eql (aref string1 pos1) #\-)
	  (setq pos2 (or (position #\- string2 :start pos2)
			 (return nil))))
	(when (eql (aref string2 pos2) #\-)
	  (setq pos1 (or (position #\- string1 :start pos1)
			 (return nil))))
	(unless (char-equal (aref string1 pos1) (aref string2 pos2))
	  (return nil))
	finally (return t)))

(defun match-input (input leaves)
  ;; precise match
  (dolist (leaf leaves)
    (when (string-equal (name leaf) input)
      (return-from match-input leaf)))
  ;; fuzzy (nonempty) match
  (dolist (leaf leaves)
    (when (and (fuzzy-match-strings (name leaf) input)
	       (plusp (length input)))
      (return-from match-input leaf))))


(defun show-demo (demo)
  "Shows the given demo."
  (with-slots (long short execute leaves)
    demo
    (format t "~&~64,,,'*<~>~%~%")
    (when short (format t "~A~%~%"
                        (typecase short
                         (function (funcall short))
                         (t short))))
    (when long (format t "~A~%~%"
                       (typecase long
                         (function (funcall long))
                         (t long))))
    (when execute (funcall execute))
    (when (or leaves (eq demo *demo-root*))
      (loop
       (loop for leaf in leaves
	     for name = (name leaf) do
	     (format t "~&~C ~A - ~A~%"
		     (if (gethash name *visited-demos*) #\+ #\*)
		     name (slot-value leaf 'short)))
       (format t "~&~%Your choice (? for help): ")
       (finish-output)
       (let ((input (read-line *standard-input* nil)))
	 (cond
	   ((or (null input) (string-equal input "quit"))
	    (update-demo-status demo)
	    (throw 'quit nil))
	   ((member input '("back" "up") :test #'equalp)
	    (update-demo-status demo)
	    (return))
	   ((member input '("help" "?") :test #'equalp)
	    (format t "~&~A~%~%" *demo-help*)
	    (format t "Type ENTER to continue.~%")
	    (finish-output)
	    (read-line))
	   ((match-input input leaves)
	    (show-demo (match-input input leaves)))
	   (t (format t "~&There is no demo with this name here.  Try again.~%~%"))))
       (format t "~&~64,,,'*<~>~%~%")
       (when short (format t "~A~%~%" short))))
    (update-demo-status demo)))

(defun femlisp-demo (&optional (demo *demo-root*))
  "Shows all demos below the given demo root."
  (catch 'quit (show-demo demo))
  (values))

(defun adjoin-demo (demo parent)
  "Adjoins the demo @arg{demo} to @arg{parent}."
  (let ((tail (member (name demo) (leaves parent) :test #'string-equal
		      :key #'name)))
    (if tail
	(setf (car tail) demo)
	(setf (leaves parent) (nconc (leaves parent) (list demo))))
    nil))

(defun remove-demo (demo parent)
  "Remove @arg{demo} from @arg{parent}."
  (setf (leaves parent)
	(unless (eq demo :all)
	  (if (stringp demo)
	      (remove demo (leaves parent) :test #'string-equal :key #'name)
	      (remove demo (leaves parent))))))

(defun find-demo (name parent)
  (find name (leaves parent) :test #'string-equal :key #'name))

(defparameter *user-input-stream* *standard-input*
  "This is an input stream for user input which is normally
bound to *standard-input*.  During testing this can be bound to
a string stream to provide sample input.")

(defun user-input (prompt &optional (converter #'identity) (test-p (constantly t)))
  "User input for demo functions.  Reads lines until
@arg{test-p} returns t on the item read."
  (loop for line = (progn (format t "~&~A" prompt)
			  (finish-output)
			  (read-line *user-input-stream*)) do
       (cond
	 ((string-equal line "quit") (throw 'quit nil))
	 ((string-equal line "up") (return :up))
	 (t (ignore-errors
	      (let ((item (funcall converter line)))
		(when (funcall test-p item)
		  (return item))))))))

(defun user-input-textfield (&optional (test-p (constantly t)))
  "User textfield input for demo functions.  Reads and concatenates lines until
an empty line is read."
  (loop 
    (format t "~&Textfield input (finished by empty line):~%")
    (finish-output)
    (let ((result
            (format nil "~{~A~^~%~}" 
                    (loop for line = (read-line *user-input-stream*)
                          until (or (null line) (zerop (length line)))
                          collect line))))
      (when (funcall test-p result)
        (return result)))))

(defun extract-demo-strings (string &optional translations)
  "Extract demo information from the documentation string of the
generating function."
  (let ((results
	 ;; #+cl-ppcre (nth-value 1 (cl-ppcre:scan-to-strings "(.*) - (.*)\\s*((?s).*)" string))
	 (let ((pos1 (search " - " string))
	       (pos2 (search (format nil "~%~%") string)))
	   (list (subseq string 0 pos1)
		 (and pos1 (subseq string (+ pos1 3) pos2))
		 (and pos2 (subseq string (+ pos2 2)))))))
    (values-list (map 'list (rcurry #'translate translations) results))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Usage
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; The following usage has turned out to yield good results in
;;; more complicated cases:

;;; First a function is written to do the demo, e.g.

#+(or)
(defun do-the-demo (parameters)
  "Documentation"
  ...)

;;; Then a demo performing this function and using its
;;; documentation is generated and put into the demo graph.
#+(or)
(let ((demo
	(make-demo
	 :name "..."
	   :short "..."
	   :long
	   (format nil "~A~%~%Parameters: ... ~%~%"
		   (documentation 'do-the-demo 'function)
		   parameters ...)
	  :execute (lambda () (do-the-demo parameters)))))
    (adjoin-demo demo *where-you-want-it-1*)
    (adjoin-demo demo *where-you-want-it-2*) ...)

;;; For multiple use, this may be wrapped into a function or
;;; macro, e.g.

#+(or)
(defun make-cdr-bl-demo (dim)
  (let ((demo
	  (make-demo
	   :name (format nil "name-~DD" dim)
	   :short "Diffusion problem on a domain with oscillating boundary"
	   :long
	   (format nil "~A~%~%Parameters: dim=~D: ~A~%"
		   (documentation 'cdr-bl-computation 'function)
		   dim)
	   :execute (lambda () (cdr-bl-computation dim)))))
    (adjoin-demo demo *laplace-demo*)
    (adjoin-demo demo *adaptivity-demo*)
    (adjoin-demo demo *boundary-coeffs-demo*)))


;;;; Testing: 

(defun test-all-demos (&optional (demo *demo-root*))
  "Performs all demos reachable from @arg{demo}."
  (let ((*visited-demos* *visited-demos* #+(or) (make-hash-table :test 'equal))
	(buggy-demos ()))
    (with-slots (name long short execute leaves input)
	demo
      (format t "~&~64,,,'*<~>~%~%")
      (when short (format t "~A~%~%" short))
      (when long (format t "~A~%~%" long))
      (update-demo-status demo)
      (append
       (when execute
	 (let ((*user-input-stream*
		(and input (make-string-input-stream input))))
	   (when (catch 'trap
		   (handler-bind ((error #'(lambda (condition) (throw 'trap condition)))
                                  (warning #'(lambda (condition) (throw 'trap condition))))
		     (format t "~&~%***** Testing ~A *****~%~%" name)
		     (funcall execute)
		     nil))
	     (push demo buggy-demos))))
       (loop for leaf in leaves
	     for name = (name leaf)
	     unless (gethash name *visited-demos*)
	     nconcing (test-all-demos leaf))))
    buggy-demos))
  
(defun test-demo ()
  (assert (not (fuzzy-match-strings "unit-interval" "--2")))
  ;; generate and delete a demo
  (let ((demo (make-demo :name "hello" :short "Hello, world!" :long "Simple test"
			 :execute (lambda () (princ "Hello, world!~%")))))
    (adjoin-demo demo *demo-root*)
    (remove-demo "hello" *demo-root*))
  )

;;; (fl.demo::test-demo)
(fl.tests:adjoin-test 'test-demo)
