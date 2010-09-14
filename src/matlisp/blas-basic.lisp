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
;;; float scalars
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun float-accuracy (type)
  (or (if (symbolp type)
          (case type
            (single-float 'single-float)
            (double-float 'double-float))
          (when (complex-type type)
            (float-accuracy (second type))))
      'number))

(defun complex-type (type)
  (and (listp type) (eq (first type) 'complex)))

(defun float-type-union (type1 type2)
  (let ((accuracy-1 (float-accuracy type1))
        (accuracy-2 (float-accuracy type2)))
    (cond ((eql accuracy-1 'number) 'number)
          ((eql accuracy-2 'number) 'number)
          (t (let* ((accuracy
                     (if (or (eql accuracy-1 'double-float)
                             (eql accuracy-2 'double-float))
                         'double-float
                         'single-float)))
               (if (or (complex-type type1) (complex-type type2))
                   (list 'complex accuracy)
                   accuracy))))))

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
      (let (matrices number-args numbers element-types)
	(loop for template-arg in primary-args
	   and actual-arg in actual-args
	   for arg = (if (consp template-arg)
			 (car template-arg)
			 template-arg)
	   do
	   (cond
	     ((matrix-p template-arg)
	      (push (element-type actual-arg) element-types)
	      (push arg matrices))
	     ((number-p template-arg)
	      (push arg number-args)
	      (push actual-arg numbers))
	     ;; otherwise: do nothing
	     ))
        (let ((element-type (reduce #'float-type-union element-types)))
          (dolist (number numbers)
            (unless (subtypep (type-of number) element-type)
              (error "Type of float number has to fit with the most general matrix
	    element-type.")))
          `(defmethod
               ,name
               (,@(loop for template-arg in primary-args
                     and actual-arg in actual-args
                     collect
		     (if (matrix-p template-arg)
			 `(,(car template-arg) ,(class-name (class-of actual-arg)))
			 template-arg))
                ,@rest-args)
             (declare (type ,element-type ,@number-args))
             #+lispworks (declare (optimize (float 0)))
             (very-quickly
               ,(subst element-type 'element-type
                       `(macrolet ,blas-macros
                          ,@body)))))))))

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
	    (fl.port:compile-and-eval
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
    (format t "~A-~D: ~$ MFLOPS~%" fsym size
	    (multiple-value-bind (time count)
		(fl.utilities::measure-time-repeated
		 (lambda ()
		   (if z
		       (funcall fn x y z)
		       (funcall fn x y))))
	      (/ (* (funcall flop-calculator size)
		    count)
		 (* 1e6 time))))))

;;; Testing

(defun test-blas-basic ()
  (float-type-union 'single-float 'number)
  )