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

(in-package :fl.macros)

(defmacro with-gensyms (syms &body body)
  "From Graham's book."
  `(let ,(mapcar #'(lambda (s) `(,s (gensym ,(symbol-name s))))
		 syms)
     ,@body))

(defun symconc (&rest args)
  "This function builds a symbol from its arguments and interns it.  This
is for use in some macros."
  (intern
   (apply #'concatenate 'string
	  (mapcar #'(lambda (arg) (if (symbolp arg)
				      (symbol-name arg)
				      arg))
		  args))
   (let ((sym (find-if #'symbolp args)))
     (if sym
	 (symbol-package sym)
	 *package*))))

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

(defmacro acond (&rest clauses)
  (if (null clauses)
      nil
      (let ((cl1 (car clauses))
            (sym (gensym)))
        `(let ((,sym ,(car cl1)))
           (if ,sym
               (let ((it ,sym))
		 (declare (ignorable it))
		 ,@(cdr cl1))
               (acond ,@(cdr clauses)))))))

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

(defmacro remove-this-method (gf-name &rest rest)
  "It should be possible to use this directly on a copied first line of a
DEFMETHOD definition, e.g.:
> (remove-this-method m* :before ((mat <matrix>) (x <vector>)))"
  (let ((next (first rest)))
    (multiple-value-bind (qualifiers args)
	(if (member next '(:before :after :around))
	    (values (list next) (second rest))
	    (values () next))
      (let ((specializers (mapcar #'(lambda (arg)
				      (if (consp arg)
					  (second arg)
					  t))
				  args)))
	`(remove-method
	  (function ,gf-name)
	  (find-method (function ,gf-name) (quote ,qualifiers)
	   (mapcar #'find-class (quote ,specializers))))))))

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
  (let ((fixnum-vec '(simple-array fixnum (*))))
    (with-gensyms
	(inc! begin end inside)
      `(let ((,begin (coerce ,start ',fixnum-vec))
	     (,end (coerce ,stop ',fixnum-vec)))
	;;(declare (type ,fixnum-vec ,begin ,start))
	(flet ((,inc! (x)
		 (declare (type ,fixnum-vec x))
		 (dotimes (i (length ,begin))
		   (cond ((< (aref x i) (aref ,end i))
			  (incf (aref x i))
			  (return t))
			 (t (setf (aref x i) (aref ,begin i)))))))
	  (do ((,var (copy-seq ,begin))
	       (,inside (every #'<= ,begin ,end) (,inc! ,var)))
	      ((not ,inside)) ,@body))))))

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
  )



