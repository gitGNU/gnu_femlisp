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

(defpackage "FL.MACROS"
  (:use "COMMON-LISP")
  (:export
   "WITH-GENSYMS" "SYMCONC" "AWHEN" "WHEREAS" "AIF"
   "AAND" "ACOND" "_F" "IT" "ENSURE" "REMOVE-THIS-METHOD"
   "FOR" "FOR<" "MULTI-FOR" "DEFINLINE"
   "?1" "?2" "?3"
   "DELAY" "FORCE"
   "FLUID-LET")
  (:documentation
   "This package contains some basic macro definitions used in Femlisp."))

(in-package :fl.macros)

(defmacro with-gensyms (syms &body body)
  "Standard macro providing the freshly generated symbols @arg{syms} to the
code in @arg{body}."
  `(let ,(mapcar #'(lambda (s) `(,s (gensym ,(symbol-name s))))
		 syms)
    ,@body))

(defun symconc (&rest args)
  "This function builds a symbol from its arguments and interns it.  This
is used in some macros."
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

(defmacro whereas (clauses &rest body)
  "Own implementation of the macro @function{whereas} suggested by Erik
Naggum (c.l.l., 4.12.2002)."
  (if (null clauses)
      `(progn ,@body)
      (destructuring-bind ((var expr &optional type) . rest-clauses)
	  clauses
	`(let ((,var ,expr))
	  ,(when type `(declare (type ,type ,var)))
	  (when ,var
	    (whereas ,rest-clauses ,@body))))))

(defmacro aif (test-form then-form &optional else-form)
  `(let ((it ,test-form))
     (if it ,then-form ,else-form)))

(defmacro awhen (test-form &body body)
  "Anaphoric macro from @cite{(Graham XXX)}."
  `(let ((it ,test-form))
     (when it ,@body)))

(defmacro awhile (expr &body body)
  "Anaphoric macro from @cite{(Graham XXX)}."
  `(do ((it ,expr ,expr))
    ((not it))
    ,@body))

(defmacro aand (&rest args)
  "Anaphoric macro from @cite{(Graham XXX)}."
  (cond ((null args) t)
	((null (cdr args)) (car args))
	(t `(aif ,(car args) (aand ,@(cdr args))))))

(defmacro acond (&rest clauses)
  "Anaphoric macro from @cite{(Graham XXX)}."
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
  "Macro from @cite{(Graham XXX)}.  Turns the operator @arg{op} into a
modifying form, e.g. @code{(_f + a b) @equiv{} (incf a b)}."
  (multiple-value-bind (vars forms var set access) 
                       (get-setf-expansion place)
    `(let* (,@(mapcar #'list vars forms)
            (,(car var) (,op ,access ,@args)))
       ,set)))

(define-modify-macro ensure (&rest args) or
   "Ensures that some place is set.  It expands as follows:
@lisp
  (ensure place value) @expansion{} (or place (setf place value))
@end lisp
It is not clear if the implementation used here works everywhere.  If not,
the workaround below the macro definition should be used.")

#+(or)
(defmacro ensure (place newval &environment env)
  "Essentially (or place (setf place newval)).  Posted by Erling Alf to
c.l.l. on 11.8.2004, implementing an idea of myself posted on c.l.l. on 30
Jul 2004 in a probably more ANSI conforming way."
  (multiple-value-bind (vars vals putvars putform getform) 
      (get-setf-expansion place env)
    `(let* ,(mapcar #'list vars vals)
       (or ,getform
	   (multiple-value-bind ,putvars
	       ,newval
	     ,putform)))))

;;; Others

(defmacro short-remove-method (gf-name qualifiers specializers)
  "Removes the method for the generic function @arg{gf-name} which is
specified by @arg{qualifiers} and @arg{specializers}.  Example:
@lisp
  (short-remove-method m* (:before) (<sparse-matrix> <sparse-vector>))
@end lisp"
  `(remove-method
    (function ,gf-name)
    (find-method (function ,gf-name) (quote ,qualifiers)
     (mapcar #'find-class (quote ,specializers)))))

(defmacro remove-this-method (gf-name &rest rest)
  "Removes the method for the generic function @arg{gf-name} which is
specified by @arg{qualifiers} and @arg{specializers}.  Example:
@lisp
  (remove-this-method m* :before ((mat <matrix>) (x <vector>)))
@end lisp
It should be possible to use this directly on a copied first line of a
DEFMETHOD definition."
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
  "Loops for @arg{var} from @arg{start} upto @arg{end}."
  (let ((limit (gensym)))
    `(let ((,limit ,end))
       (do ((,var ,start (+ ,var 1)))
	   ((> ,var ,limit))
	 ,@body))))

(defmacro for< ((var start end) &body body)
  "Loops for @arg{var} from @arg{start} below @arg{end}."
  (let ((limit (gensym)))
    `(let ((,limit ,end))
       (do ((,var ,start (+ ,var 1)))
	   ((>= ,var ,limit))
	 ,@body))))

(defmacro multi-for ((var start stop) &body body)
  "Loops for @arg{var} being an integer vector starting from @arg{start}
below @arg{end}.  Example:
@lisp
  (multi-for (x #(1 1) #(3 3)) (princ x) (terpri))
@end lisp"
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
  "Delays the evaluation of @arg{form}."
  (with-gensyms (value computed)
    `#'(lambda ()
	 (let ((,computed nil)
	       (,value nil))
	   (if ,computed
	       ,value
	       (prog1
		   (setq ,value ,form)
		 (setq ,computed t)))))))

(defmacro force (delayed-form)
  "Forces the value of a @arg{delayed-form}."
  (with-gensyms (form)
    `(let ((,form ,delayed-form))
      (if (functionp ,form) (funcall ,form) ,form))))

(defmacro definline (name &rest rest)
  "Short form for defining an inlined function.  It should probably be
deprecated, because it won't be recognized by default by editors."
  `(progn
    (declaim (inline ,name))
    (defun ,name ,@rest)
    ))

;;; some macros for choosing between possibilities
(defmacro ?1 (&rest args)
  "A macro returning the first of its arguments."
  (first args))
(defmacro ?2 (&rest args)
  "A macro returning the second of its arguments."
  (second args))
(defmacro ?3 (&rest args)
  "A macro returning the third of its arguments."
  (third args))

(defmacro fluid-let (bindings &body body)
  "Sets temporary bindings."
  (let ((identifiers (mapcar #'car bindings)))
    (with-gensyms (saved)
      `(let ((,saved (list ,@identifiers)))
	(unwind-protect
	     (progn ,@(loop for (place value) in bindings collect
			    `(setf ,place ,value))
		    ,@body)
	  ,@(loop for (place value) in bindings and i from 0 collect
		  `(setf ,place (nth ,i ,saved))))))))

;;;; Testing:
(defun test-macros ()
  (let ((x 5))
    (ensure x 1))
  (let ((a 1) (b 2))
    (fluid-let ((a 3) (b 4))
      (list a b))))




