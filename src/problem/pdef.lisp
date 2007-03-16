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
			  `(extract-from (,(if grad-p 'second 'first) solution) ,from ,ncomps ,flag)))))))))))
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
;;;; nonlinear right-hand side for <ellsys-problem>
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defparameter *linearization-factor* 1.0
  "Damping factor used for linearizing nonlinear right-hand sides.  0.0
yields the fixed point iteration, 1.0 is full Newton approximation.")

(defmethod make-coefficients-for ((problem fl.ellsys::<ellsys-problem>)
				  (coeff (eql 'FL.ELLSYS::NONLINEAR-F))
				  patch args eval)
  (let ((grad (sparse-real-derivative eval)))
    (append
     (make-coefficients-for
      problem 'FL.ELLSYS::R patch args
      (lambda (&rest args)
	;; Transform grad from args to components!
	(scal (- *linearization-factor*) (apply grad args))))
     (make-coefficients-for
      problem 'FL.ELLSYS::F patch args
      (lambda (&rest args)
	(let* ((u (coerce args 'vector))
	       (f (apply eval args))
	       (Df (apply grad args)))
	  (axpy (-  *linearization-factor*) (m* Df u) f)))))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; cdr problems
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmethod make-coefficients-for ((problem fl.cdr::<cdr-problem>)
				  (coeff (eql 'FL.CDR::DIFFUSION))
				  patch args eval)
  (call-next-method problem 'FL.ELLSYS::A patch args
		    (lambda (&rest args)
		      (diagonal-sparse-tensor (apply eval args) 1))))
		    
(defmethod make-coefficients-for ((problem fl.cdr::<cdr-problem>)
				 (coeff (eql 'FL.CDR::ISOTROPIC-DIFFUSION))
				  patch args eval)
  (call-next-method problem 'FL.ELLSYS::A patch args
		    (let ((dim (dimension (domain problem))))
		      (lambda (&rest args)
			(diagonal-sparse-tensor (scal (apply eval args) (eye dim)) 1)))))

(defmethod make-coefficients-for ((problem fl.cdr::<cdr-problem>)
				 (coeff (eql 'FL.CDR::SCALAR-SOURCE)) patch args eval)
  (call-next-method problem 'FL.ELLSYS::F patch args
		    (lambda (&rest args)
		      (vector (ensure-matlisp (apply eval args) :row-vector)))))
		    
(defmethod make-coefficients-for ((problem fl.cdr::<cdr-problem>)
				 (coeff (eql 'FL.CDR::SCALAR-NONLINEAR-SOURCE))
				  patch args eval)
  (call-next-method problem 'FL.ELLSYS::NONLINEAR-F patch args
		    (lambda (&rest args)
		      (vector (ensure-matlisp (apply eval args) :row-vector)))))
		    
(defmethod make-coefficients-for ((problem fl.cdr::<cdr-problem>)
				 (coeff (eql 'FL.CDR::SCALAR-CONSTRAINT))
				  patch args eval)
  (call-next-method problem 'FL.PROBLEM::CONSTRAINT patch args
		    (lambda (&rest args)
		      (values #(t) (vector (ensure-matlisp (apply eval args) :row-vector))))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Navier-Stokes problems
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(in-package :fl.navier-stokes-ellsys)

(defmethod make-coefficients-for
    ((problem fl.navier-stokes-ellsys::<navier-stokes-problem>)
     (coeff (eql 'VISCOSITY)) patch args eval)
  (let ((dim (dimension (domain problem))))
    (call-next-method problem 'FL.ELLSYS::A patch args
		      (lambda (&rest args)
			(diagonal-sparse-tensor
			 (scal (apply eval args) (eye dim)) dim)))))

(defmethod make-coefficients-for
    ((problem fl.navier-stokes-ellsys::<navier-stokes-problem>)
     (coeff (eql 'REYNOLDS)) patch args eval)
  (assert (null args))
  (let ((reynolds (funcall eval))
	(dim (dimension (domain problem))))
    (list (make-coefficients-for
	   problem 'FL.ELLSYS::C patch '(u)
	   (lambda (u)
	     (diagonal-sparse-tensor (scal (* *alpha* reynolds) u) dim)))
	  (make-coefficients-for
	   problem 'FL.ELLSYS::R patch '(du)
	   (lambda (du)
	     (scal (* *beta* reynolds) du)))
	  (make-coefficients-for
	   problem 'FL.ELLSYS::F patch '(u du)
	   (lambda (u du)
	     (scal (* (+ *alpha* *beta* -1) reynolds) (m* du u)))))))

(defmethod make-coefficients-for ((problem fl.navier-stokes-ellsys::<navier-stokes-problem>)
				 (coeff (eql 'FL.NAVIER-STOKES-ELLSYS::FORCE))
				  patch args eval)
  (make-coefficients-for problem 'FL.ELLSYS::F patch args
			(lambda (&rest args)
			  (vector (apply eval args)))))

(defmethod make-coefficients-for ((problem fl.navier-stokes-ellsys::<navier-stokes-problem>)
				 (coeff (eql 'FL.NAVIER-STOKES-ELLSYS::PRESCRIBED-VELOCITY))
				 patch args eval)
  (let ((dim (dimension (domain problem)))
	(multiplicity (multiplicity problem)))
    (make-coefficients-for
     problem 'FL.PROBLEM::CONSTRAINT patch args
     (lambda (&rest args)
       (let ((velocity (apply eval args))
	     (values (make-array (1+ dim) :initial-element (zeros 1 multiplicity)))
	     (flags (make-array (1+ dim) :initial-element nil)))
	 (dotimes (i dim)
	   (mextract (aref values i) velocity i 0)
	   (setf (aref flags i) t))
	 (values flags values))))))
		    
(defmethod make-coefficients-for ((problem fl.navier-stokes-ellsys::<navier-stokes-problem>)
				 (coeff (eql 'FL.NAVIER-STOKES-ELLSYS::NO-SLIP))
				 patch args eval)
  (declare (ignore eval))
  (make-coefficients-for
   problem 'FL.NAVIER-STOKES-ELLSYS::PRESCRIBED-VELOCITY patch args
   (constantly (zeros (dimension (domain problem)) (multiplicity problem)))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; problem definition macro
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(in-package :fl.problem)

(defmacro coeff (name args &body body)
  "This macro defines a coefficient function inside a problem definition.
It is mainly defined such that the Lisp editor indents the function
definition correctly."
  (declare (ignore name args body))
  (error "Should not be called."))

(defmacro create-problem (&key type domain components (multiplicity 1) coefficients)
  (with-gensyms (problem patch classification)
    `(lret ((,problem (fl.amop:make-programmatic-instance
		       ,type :domain ,domain :multiplicity ,multiplicity
		       :classify nil)))
      ,(if (eq (first components) 'select-on-patch)
	     `(progn
	       (setf (slot-value ,problem 'components) (make-hash-table))
	       (doskel (,patch (domain ,problem))
		 (let ((,classification (patch-classification ,patch (domain ,problem))))
		   (setf (components-of-patch ,patch ,problem)
			 (cond
			   ,@(loop for (condition comps) in (rest components) collect
				   `((test-condition ',condition ,classification) ,comps)))))))
	     `(setf (slot-value ,problem 'components) ,components))
      ,(if (eq (first coefficients) 'select-on-patch)
	   `(progn
	     (setf (slot-value ,problem 'coefficients) (make-hash-table))
	     (doskel (,patch (domain ,problem))
	       (let ((,classification (patch-classification ,patch (domain ,problem))))
		 (setf (coefficients-of-patch ,patch ,problem)
		       (cond
			 ,@(loop for (condition . coeffs) in (rest coefficients) collect
				 `((test-condition ',condition ,classification)
				   (append
				    ,@(loop for (coeff name args . body) in coeffs 
					    do (assert (string-equal (symbol-name coeff) "COEFF"))
					    collect
					    `(make-coefficients-for
					      ,problem ',name ,patch
					      ',args (lambda ,args ,@body)))))))))))
	   `(setf (slot-value ,problem 'coefficients) ,coefficients))
      ;;
      (classify-problem ,problem)
      )))

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
	     )))
    (let ((domain (n-cube-domain 1)))
      (make-coefficients-for
       (make-instance 'fl.cdr::<cdr-problem> :domain domain)
       'fl.cdr::isotropic-diffusion
       (central-patch domain) '() (constantly 1.0))))
    
  #+(or)
  (let ((problem
	 (create-problem
	  :type 'FL.CDR::<CDR-PROBLEM>
	  :domain (n-cube-domain 2) :components '(u) :multiplicity 1
	  :coefficients
	  (select-on-patch
	   (:d-dimensional
	    (coeff FL.CDR::ISOTROPIC-DIFFUSION () 1.0)
	    (coeff FL.CDR::SCALAR-SOURCE () 1.0))
	   (:external-boundary
	    (coeff FL.CDR::SCALAR-CONSTRAINT () 0.0))))))
    (defparameter *result*
      (solve (blackboard :problem problem :solver (fl.iteration::lu-solver)
			 :success-if '(> :step 3) :output :all))))

  #+(or)
  (fl.plot:plot (getbb *result* :solution))

  ;; Navier-Stokes
  #+(or)
  (let ((problem
	 (create-problem
	  :type 'FL.NAVIER-STOKES-ELLSYS::<NAVIER-STOKES-PROBLEM>
	  :domain (n-cube-domain 2) :components '((u 2) p) :multiplicity 1
	  :coefficients
	  (select-on-patch
	   (:d-dimensional
	    (coeff FL.NAVIER-STOKES-ELLSYS::VISCOSITY () 1.0)
	    (coeff FL.NAVIER-STOKES-ELLSYS::REYNOLDS () 1.0))
	   ((and :d-1-dimensional :upper)
	    (coeff FL.NAVIER-STOKES-ELLSYS::PRESCRIBED-VELOCITY ()
	      #m((1.0) (0.0))))
	   (:external-boundary
	    (coeff FL.NAVIER-STOKES-ELLSYS::NO-SLIP ()))))))
    (describe problem))

  (print(macroexpand-1
   '(create-problem
     :type 'FL.NAVIER-STOKES-ELLSYS::<NAVIER-STOKES-PROBLEM>
     :domain (n-cube-domain 2) :components '((u 2) p) :multiplicity 1
     :coefficients
     (select-on-patch
      (:d-dimensional
       (coeff FL.NAVIER-STOKES-ELLSYS::VISCOSITY () 1.0)
       (coeff FL.NAVIER-STOKES-ELLSYS::REYNOLDS () 1.0))
      ((and :d-1-dimensional :upper)
       (coeff FL.NAVIER-STOKES-ELLSYS::PRESCRIBED-VELOCITY ()
	 #m((1.0) (0.0))))
      (:external-boundary
       (coeff FL.NAVIER-STOKES-ELLSYS::NO-SLIP ()))))))

  )

;;; (test-pdef)
(fl.tests:adjoin-test 'test-pdef)