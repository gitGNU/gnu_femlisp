;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; macros.lisp - Useful macro definitions
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

(in-package :macros)

(defmacro with-gensyms (syms &body body)
  "From Graham's book."
  `(let ,(mapcar #'(lambda (s) `(,s (gensym ,(symbol-name s))))
		 syms)
     ,@body))

(defmacro aif (test-form then-form &optional else-form)
  `(let ((it ,test-form))
     (if it ,then-form ,else-form)))

(defmacro whereas (clauses &rest body)
  "Own implementation of Erik Naggum's whereas-macro (c.l.l., 4.12.2002)."
  (if (null clauses)
      `(progn ,@body)
      (destructuring-bind ((var expr &optional type) . rest-clauses)
	  clauses
	`(let ((,var ,expr))
	  ,(when type `(declare (type ,type ,var)))
	  (when ,var
	    (whereas ,rest-clauses ,@body))))))

(defmacro awhen (test-form &body body)
  `(let ((it ,test-form))
     (when it ,@body)))

(defmacro awhile (expr &body body)
  `(do ((it ,expr ,expr))
    ((not it))
    ,@body))

(defmacro aand (&rest args)
  (cond ((null args) t)
	((null (cdr args)) (car args))
	(t `(aif ,(car args) (aand ,@(cdr args))))))

(defmacro _f (op place &rest args)
  (multiple-value-bind (vars forms var set access) 
                       (get-setf-expansion place)
    `(let* (,@(mapcar #'list vars forms)
            (,(car var) (,op ,access ,@args)))
       ,set)))

;;; Others

(defmacro short-remove-method (gf-name qualifiers specializers)
  "Syntax: (short-remove-method m* (:before) (<sparse-matrix> <sparse-vector>))"
  `(remove-method
    (function ,gf-name)
    (find-method (function ,gf-name) (quote ,qualifiers)
     (mapcar #'find-class (quote ,specializers)))))

(defun with-items-expander-1 (prop-vars prop-names prop-defaults plist body)
  "Expander for with-items.  Does not work because assembly-line is not
modified immediately."
  `(let ,(mapcar #'(lambda (var prop default)
		     `(,var (getf ,plist ',prop ,default)))
		 prop-vars prop-names prop-defaults)
    ,@body
    ,@(mapcar #'(lambda (var prop)
		  `(setf (getf ,plist ',prop) ,var))
	      prop-vars prop-names)))

(defun with-items-expander-2 (prop-vars prop-names prop-defaults plist body)
  "Expander for with-items."
  `(progn
    ,@(mapcar #'(lambda (prop default)
		  `(setf (getf ,plist ',prop) (getf ,plist ',prop ,default)))
	      prop-names prop-defaults)
    (symbol-macrolet
	,(mapcar #'(lambda (var prop) `(,var (getf ,plist ',prop)))
		 prop-vars prop-names)
      ,@body)))

(defmacro with-items (props assembly-line &body body)
  "Introduce property list members as variables.  Use:
   (with-items (&key sol (rhs nil rhs-high))
       assembly-line
     (setq sol rhs-high))"
  (let (key-p prop-vars prop-defaults prop-names)
    (dolist (prop props)
      (cond
	((eq prop '&key) (assert (not key-p))
	 (setq key-p t))
	(t (let ((prop-var
		  (if (atom prop)
		      prop
		      (or (nth 2 prop) (nth 0 prop))))
		 (prop-default (and (listp prop) (nth 1 prop)))
		 (prop-name (if (atom prop) prop (nth 0 prop))))
	     ;;(format t "~A ~A ~A" prop-var prop-default prop-name)
	     (push prop-var prop-vars)
	     (push prop-default prop-defaults)
	     (push (if key-p
		       (intern (symbol-name prop-name) :keyword)
		       prop-name)
		   prop-names)))))
  `(progn
    (assert (eq (car ,assembly-line) :assembly-line))
    ,(with-items-expander-2 prop-vars prop-names prop-defaults
			    `(cddr ,assembly-line) body))))

(defmacro with-properties (props plist &body body)
  "Deprecated."
  `(symbol-macrolet
    ,(if (eq (car props) '&key)
	 (mapcar #'(lambda (prop) `(,prop (getf ,plist ,(intern (symbol-name prop) :keyword)))) (cdr props))
	 (mapcar #'(lambda (prop) `(,prop (getf ,plist ',prop))) props))
    ,@body))

(defmacro for ((var start end) &body body)
  "Syntax: (for (i 1 10) (princ i)).
Loops for i from 1 to 10."
  (let ((limit (gensym)))
    `(let ((,limit ,end))
       (do ((,var ,start (+ ,var 1)))
	   ((> ,var ,limit))
	 ,@body))))

(defmacro for< ((var start end) &body body)
  "Syntax: (for< (i 1 10) (princ i)).
Loops for i from 1 to 9."
  (let ((limit (gensym)))
    `(let ((,limit ,end))
       (do ((,var ,start (+ ,var 1)))
	   ((>= ,var ,limit))
	 ,@body))))

(defmacro multi-for ((var start stop) &body body)
  "multi-for: This macro loops through vectors of (integer) values
between the (integer) vectors start and stop.
Example: (multi-for (x #(1 1) #(3 3)) (princ x) (terpri))"
  (with-gensyms
      (inc! begin end inside)
    `(let ((,begin ,start)
	   (,end ,stop))
      (flet ((,inc! (x)
	       (do ((i 0 (1+ i)))
		   ((= i (length ,begin)) nil)
		 (cond ((< (aref x i) (aref ,end i))
			(setf (aref x i) (1+ (aref x i)))
			(return t))
		       (t (setf (aref x i) (aref ,begin i)))))))
	(loop with ,var = (copy-seq ,begin)
	      for ,inside = (every #'<= ,begin ,end) then (setf ,inside (,inc! ,var))
	      until (not ,inside)
	      do ,@body)))))

;;; delay and force
(defmacro delay (form)
  (with-gensyms (value computed)
    `#'(lambda ()
	 (let ((,computed nil)
	       (,value nil))
	   (if ,computed
	       ,value
	       (prog1
		   (setq ,value ,form)
		 (setq ,computed t)))))))

(defmacro force (delayed)
  `(funcall ,delayed))

(defmacro definline (name &rest rest)
  `(progn
    (declaim (inline ,name))
    (defun ,name ,@rest)
    ))

;;; some macros for choosing between possibilities
(defmacro ?1 (&rest args) (first args))
(defmacro ?2 (&rest args) (second args))
(defmacro ?3 (&rest args) (third args))

;;;; Testing:

(defun test-macros ()
  (let ((assembly-line (list :assembly-line t)))
    (with-items (&key (a 2)) assembly-line
      (princ a)))
  )



