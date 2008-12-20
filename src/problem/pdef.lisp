;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; pdef.lisp
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;
;;; Copyright (C) 2007- Nicolas Neuss, University of Karlsruhe.
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

(in-package :fl.problem)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; defines an embedded language for simplifying problems definitions
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defgeneric make-coefficients-for (problem name patch demands evaluator)
  (:documentation "Generates a coefficient while dispatching on problem and
coefficient name.  May return a single coefficient or a list of several
coefficients."))

(defun extract-from (solution from ncomps &optional scalar-p index)
  "Extracts numbers or subvectors from the solution vector."
  (let ((multiplicity (ncols (vref solution 0))))
    ;; no scalar return value is possible for multiple solutions
    (cond
      (scalar-p (assert (or index (= multiplicity 1)))
		(tensor-ref solution from (or index 0)))
      (t (lret ((result (zeros ncomps (if index 1 multiplicity))))
	   (dotimes (i ncomps)
	     (if index
		 (setf (mref result i 0)
		       (tensor-ref solution (+ from i) index))
		 (dotimes (j multiplicity)
		   (setf (mref result i j)
			 (tensor-ref solution (+ from i) j))))))))))

(defun extract-from-gradient (gradient from ncomps)
  "Extract a numbers or subvectors from the solution vector."
  (let ((x (vref gradient from)))
    (assert (= 1 (ncols x)) ()
            "Multiplicity >1 is not allowed here - one would have to use
            general tensors here and at several other places.")
    (lret* ((dim (nrows x))
            (result (zeros ncomps dim)))
      (loop repeat ncomps
         for i from from
         for y = (vref gradient i) do
           (dotimes (j dim)
             (setf (mref result i j) (mref y j 0)))))))

(defun prepare-coefficient-arguments (components args)
  "Prepares arguments for the given coefficient function."
  (let ((source
	 `(lambda (&key global solution time &allow-other-keys)
	   (declare (ignorable global solution time))
	   (values
	    ,@(loop
	       for sym in args collecting
	       (let ((name (symbol-name sym)))
		 (cond
		   ((string-equal name "X") 'global)
		   ((string-equal name "TIME") 'time)
		   (t (let* ((grad-p (> (mismatch name "GRAD_" :test #'string-equal) 4))
			     (sym (if grad-p
				      (find-symbol (subseq name 5) (symbol-package sym))
				      sym)))
			(multiple-value-bind (from ncomps flag)
			    (extraction-information components sym)
                          (if grad-p
                              `(extract-from-gradient (second solution) ,from ,ncomps)
                              `(extract-from (first solution) ,from ,ncomps ,flag))))))))))))
    (fl.debug:dbg :compile "Compiling:~%~A" source)
    (compile nil source)))
	     
(defmethod make-coefficients-for ((problem <pde-problem>) coeff-name patch args eval)
  (let ((demands
	 (remove-duplicates
	  (mapcar (lambda (sym)
		    (let ((name (symbol-name sym)))
		      (cond
			((string-equal name "X") :global)
			((string-equal name "TIME") :time)
			(t '(:fe-parameters :solution)))))
		  args)
	  :test 'equalp))
	(prepare-args (prepare-coefficient-arguments
		       (components-of-patch patch problem) args)))
    (list
     (make-instance
      '<coefficient> :name coeff-name :dimension (dimension patch)
      :demands demands
      :eval (lambda (&rest rest)
	      (multiple-value-call
		  eval (apply prepare-args rest)))))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; problem definition macro
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(with-gensyms (internal-problem internal-patch)

  (defmacro coeff (name args &body body)
    "A local @macro{coeff} defines a coefficient function inside
@macro{setup-coefficients}.  It is defined here at the toplevel such that
the Lisp editor indents the definitions correctly."
    `(make-coefficients-for
      ,internal-problem ',name ,internal-patch ',args
      (lambda ,args ,@body)))

  (defmacro select-on-patch ((patch) &body clauses)
    `(cond
       ,@(loop for (pattern . coeffs) in clauses collect
	      `((patch-if ,patch ',pattern) ,@coeffs))))
  
  (defmacro setup-components ((patch) &body patch-definitions)
    "Defines components dispatching on @arg{patch}."
    `(setf (slot-value ,internal-problem 'components)
	 (memoize-1 (lambda (,internal-patch)
		      (let ((,patch ,internal-patch))
			,@patch-definitions)))))

  (defmacro setup-coefficients ((patch) &body patch-definitions)
    "Defines coefficients dispatching on @arg{patch}."
    `(setf (slot-value ,internal-problem 'coefficients)
           (memoize-1
            (lambda (,internal-patch)
              (flatten (let ((,patch ,internal-patch))
                         ,@patch-definitions))))))

  (defmacro create-problem (type (&key domain components (multiplicity 1) properties)
			    &body body)
    "Creates a PDE problem.  @arg{type} is the type of the problem which
can be the name of a problem class or a list of class names.  @arg{domain}
is the domain for this problem, @arg{multiplicity} is the multiplicity of
the solution, e.g. the number of eigenvectors we search for.  In
@arg{body}, patch-dependent coefficients should be defined with
@macro{setup-coefficients}.  It is also possible to define patch-dependent
components with @macro{setup-components}."
  `(lret ((,internal-problem
	   (fl.amop:make-programmatic-instance
	    ,type :domain ,domain :multiplicity ,multiplicity
	    :classify nil :properties ,properties
	    :components ,components)))
     (flet ((patch-if (patch condition)
              (test-condition
               condition (patch-classification
                          patch (domain ,internal-problem)))))
       ;; the missing of a call to this local function should not be
       ;; reported, even if it is unusual for coefficients not to depend on
       ;; a patch
       (declare (ignorable #'patch-if))
       ,@body)
     (classify-problem ,internal-problem)
     ))
  )

;;;; Testing
(defun test-pdef ()

  (extract-from (vector (ones 1 2) (zeros 1 2)) 0 1 nil 1)

  ;; general problem
  (flet ((central-patch (domain)
	   (find-cell (lambda (cell) (= (dimension cell) (dimension domain)))
		      domain)))
    (let* ((components '((u 3) p))
	   (domain (n-cube-domain 3))
	   (problem (make-instance '<pde-problem> :domain domain :components components))
	   (coeffs (make-coefficients-for 
		    problem 'X (central-patch domain) '(x time u p)
		    (lambda (x time u p)
		      (format t "x=~A, t=~A, u=~A, p=~A~%" x time u p)))))
      (evaluate
       (first coeffs)
       (list :global #(1.0) :time 10.0
	     :solution (list (vector #m(0.0) #m(0.0) #m(0.0) #m(0.0))
			     (vector #m(0.0 0.0) #m(0.0 0.0) #m(0.0 0.0) #m(0.0 0.0)))
	     ))))
    
  )

;;; (test-pdef)
(fl.tests:adjoin-test 'test-pdef)