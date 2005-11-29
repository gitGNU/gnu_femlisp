;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; blas-macros.lisp - Macros for defining BLAS routines 
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

(in-package :fl.matlisp)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; BLAS macros
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defvar *blas-macro-table* (make-hash-table)
  "Table of all BLAS macros.")


(defgeneric blas-macro (object name)
  (:documentation "Interns a BLAS macro for instances of @arg{object} with
name being the symbol @arg{name}.")
  (:method (object name)
    "The default is no macro for this class."
    (declare (ignore object name))
    nil))

(defun define-blas-macro (class macro-definition)
  (let ((symbol (first macro-definition)))
    (eval `(defmethod blas-macro ((z ,class) (sym (eql ',symbol)))
	    ',macro-definition))
    (setf (gethash symbol *blas-macro-table*) t)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; define-blas-template
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun new-blas-method-code (name class-name template-args actual-args blas-macros body)
  (flet ((matrix-p (template-arg)
	   (and (consp template-arg)
		(eq class-name (second template-arg))))
	 (number-p (template-arg)
	   (and (consp template-arg)
		(eq 'number (second template-arg)))))
    (let* ((pos (or (position-if #'(lambda (x)
				     (member x '(&optional &key &allow-other-keys)))
				 template-args)
		    (length template-args)))
	   (primary-args (subseq template-args 0 pos))
	   (rest-args (subseq template-args pos)))
      (let (matrices number-args numbers specialized-class element-type)
	(loop for template-arg in primary-args
	      and actual-arg in actual-args
	      for arg = (if (consp template-arg)
			    (car template-arg)
			    template-arg)
	      do
	      (cond
		((matrix-p template-arg)
		 (cond (specialized-class
			(unless (eq specialized-class (class-of actual-arg))
			  (error "Template depends on different classes.")))
		       (t
			(setq specialized-class (class-of actual-arg))
			(setq element-type (element-type actual-arg))))
		 (push arg matrices))
		((number-p template-arg)
		 (push arg number-args)
		 (push actual-arg numbers))
		;; otherwise: do nothing
		))
	(dolist (number numbers)
	  (unless (subtypep (type-of number) element-type)
	    (error "Type of number does not fit with element-type.")))
	`(defmethod
	  ,name
	  (,@(loop for arg in primary-args collect
		   (if (matrix-p arg)
		       `(,(car arg) ,(class-name specialized-class))
		       arg))
	   ,@rest-args)
	  (declare (type ,element-type ,@number-args))
	  (declare (optimize (speed 3) (safety 0) (debug 0)))
	  ,(subst element-type 'element-type
		  `(macrolet ,blas-macros
		    ,@body)))))))

(eval-when (:compile-toplevel :load-toplevel :execute)

  (defun new-dispatcher-code (name template-args body)
    "Generates a method which generates code for a type-specialized method."
    (let* ((class-pos (position-if
		       #'(lambda (x) (and (consp x)
					  (not (eq (second x) 'number))))
		       template-args))
	   (class-name (second (elt template-args class-pos)))
	   (instance-name (first (elt template-args class-pos))))
      (whereas ((gf (and (fboundp name) (symbol-function name))))
	(fl.amop:remove-subclass-methods gf template-args))
      (let* ((actual-args (gensym "ACTUAL-ARGS"))
	     (blas-macros (gensym "BLAS-MACROS")))
	`(defmethod ,name ,template-args
	  ,@(when (stringp (car body)) (list (car body)))
	  (let ((,actual-args
		 (list ,@(loop for arg in template-args
			       unless (member arg '(&optional))
			       collect (if (consp arg) (car arg) arg)
			       do (assert (not (member arg '(&key &allow-other-keys)))))))
		(,blas-macros
		 (loop for sym being each hash-key of *blas-macro-table*
		       when (blas-macro ,instance-name sym)
		       collect it)))
	    ;; define specialized method
	    (fl.amop:compile-and-eval
	     (new-blas-method-code
	      ',name ',class-name ',template-args ,actual-args ,blas-macros
	      ',(if (stringp (car body)) (cdr body) body)))
	    ;; retry call
	    (apply #',name ,actual-args)))))))

(defmacro define-blas-template (name args &body body)
  (new-dispatcher-code name args body))

;;;; Performance test

(defvar *test-blas-vector-generator* nil
  "Function of an argument @arg{N} which returns a vector of size @arg{N}
of a type that can be handled by the BLAS routines which are tested.")

(defun test-blas (fsym size &key (generator *test-blas-vector-generator*)
		  (flop-calculator (curry #'* 2)) (nr-of-arguments 2))
  "Test the performance of BLAS routines."
  (let ((x (funcall generator size))
	(y (funcall generator size))
	(z (when (= nr-of-arguments 3)
	     (funcall generator size)))
	(fn (symbol-function fsym)))
    (format
     t "~A-~D: ~$ MFLOPS~%" fsym size
     (loop with after = 0
	   for before = (get-internal-run-time) then after
	   and count of-type fixnum = 1 then (* count 2)
	   do
	   (if z
	       (loop repeat count do (funcall fn x y z))
	       (loop repeat count do (funcall fn x y)))
	   (setq after (get-internal-run-time))
	   (when (> (/ (- after before) internal-time-units-per-second)
		    fl.utilities::*mflop-delta*)
	       (return (/ (* (funcall flop-calculator size)
			     count
			     internal-time-units-per-second)
			  (* 1e6 (- after before)))))))))

