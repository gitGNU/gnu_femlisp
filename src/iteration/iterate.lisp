;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; iterate.lisp - Iteration framework
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

(in-package :iteration)

(defclass <iteration> ()
  ((initially :initform () :initarg :initially
	    :documentation "List of actions to be performed initially.")
   (intermediate :initform () :initarg :intermediate
		 :documentation "List of actions to be performed after
initialization and after each step.")
   (finally :initform () :initarg :finally
	  :documentation "List of actions to be performed finally.")
   (observe :initform () :initarg :observe
	    :documentation "The initform depends on the subclass.")
   (success-if :initform nil :initarg :success-if
	       :documentation "A form specifying a success criterion.")
   (failure-if :initform nil :initarg :failure-if
	       :documentation "A form specifying a failure criterion.")
   (output :initform nil :initarg :output :documentation "A boolean
indicating if output is to be done."))
  (:documentation "The iteration base class."))

(defgeneric name (iter)
  (:documentation "Name of the iteration."))

(defgeneric initially (iter blackboard)
  (:documentation "Performs initial operations."))

(defgeneric intermediate (iter blackboard)
  (:documentation "Is called after initialization and after each step."))

(defgeneric terminate-p (iter blackboard)
  (:documentation "Tests terminating conditions.  Returns either NIL or
:success or :failure."))

(defgeneric finally (iter blackboard)
  (:documentation "Performs final operations."))

(defgeneric next-step (iter blackboard)
  (:documentation "Does a step of the iteration."))

(defgeneric iterate (iter blackboard)
  (:documentation "Iterates on the data in the blackboard according to the
iteration iter."))

(defmethod name (iter)
  "Default name is the class name of the iteration."
  (class-name (class-of iter)))

(defmethod initially ((iter <iteration>) blackboard)
  "Default method.  Prints the header line for observed quantities and
performs all actions in initial."
  (setf (getbb blackboard :step) 0)
  (with-slots (output observe) iter
    (when output
      (indented-format t "Iteration ~A" (name iter))
      (let ((fstr (make-array '(0) :element-type 'base-char
			      :fill-pointer 0 :adjustable t)))
	(with-output-to-string (s fstr)
	  (dolist (item observe)
	    (let ((title (first item)))
	      (format s (if (functionp title) (funcall title) title))
	      (format s "  "))))
	(indented-format t "~A" fstr)
	(indented-format t "~V,,,'-<~>" (length fstr))))))

(defmethod initially :after ((iter <iteration>) blackboard)
  "Perform user-defined actions."
  (dolist (action (slot-value iter 'initially))
    (funcall action blackboard)))

(defmethod intermediate ((iter <iteration>) blackboard)
  "Default method.  Prints observed quantities and
performs all actions in intermediate."
  (with-slots (output observe) iter
    (when output
      (let ((fstr (make-array '(0) :element-type 'base-char
			      :fill-pointer 0 :adjustable t)))
	(with-output-to-string (s fstr)
	  (dolist (item observe)
	    (let ((formatter (second item))
		  (evaluator (third item)))
	      (format s formatter (funcall evaluator blackboard))
	      (format s "  "))))
	(indented-format t fstr)
	(force-output)))))

(defmethod intermediate :after ((iter <iteration>) blackboard)
  "Performs user-defined actions."
  (dolist (action (slot-value iter 'intermediate))
    (funcall action blackboard)))

(defmethod terminate-p ((iter <iteration>) blackboard)
  "Default method evaluating success-if and failure-if expressions."
  (flet ((replace-symbols (tree)
	   "Replaces all keyword symbols by blackboard values."
	   (map-tree
	    #'(lambda (atom)
		(if (and (symbolp atom)
			 (eq (symbol-package atom)
			     (find-package :keyword)))
		    (getbb blackboard atom)
		    atom))
	    tree)))
    (with-slots (success-if failure-if) iter
      (cond ((eval (replace-symbols success-if))
	     (setf (getbb blackboard :status) :success))
	    ((eval (replace-symbols failure-if))
	     (setf (getbb blackboard :status) :failure))
	    (t (setf (getbb blackboard :status) nil))))))

(defmethod next-step :after ((iter <iteration>) blackboard)
  "Increment step counter."
  (incf (getbb blackboard :step)))

(defmethod finally ((iter <iteration>) blackboard)
  "Setup a report."
  (with-items (&key report status step) blackboard
    (setq report (blackboard :status status :steps step))))

(defmethod finally :after ((iter <iteration>) blackboard)
  "Executes all actions in finally."
  (dolist (action (slot-value iter 'finally))
    (funcall action blackboard))
  (when (slot-value iter 'output)
    ;; and ensure a fresh line
    (format t "~&")))

(defmethod iterate ((iter <iteration>) blackboard)
  "Default method for performing an iteration.  Increases indentation level
and calls initialization, termination check, stepping, and finalization in
the correct order."
  (let ((*indentation-level* (1+ *indentation-level*)))
    (dbg :iter (format nil "Iteration ~A: initializing" (name iter)))
    (initially iter blackboard)
    (dbg :iter (format nil "Iteration ~A: looping" (name iter)))
    (loop (intermediate iter blackboard)
	  (when (terminate-p iter blackboard)
	    (finally iter blackboard)
	    (return blackboard))
	  (next-step iter blackboard))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Testing success and failure
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;; Testing: (test-iterate)
(defun test-iterate ()
  (let ((iter (make-instance
	       '<iteration>
	       :success-if '(and (> :step 5) (< :defnorm 1.0e-8)))))
    (describe iter)
    (terminate-p iter (blackboard :step 10 :defnorm 1.0e-9)))
  )
(tests::adjoin-femlisp-test 'test-iterate)