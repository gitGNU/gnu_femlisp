;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; pde-problem.lisp
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;
;;; Copyright (C) 2003, 2004 Nicolas Neuss, University of Heidelberg.
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

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; coefficient
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass <coefficient> ()
  ((dimension :reader dimension :initform nil :initarg :dimension :documentation
	      "The dimension of the cell on which this coefficient is
active.  The value T means that it is active on all cells lying on the
patch.")
   (demands :reader demands :initform () :initarg :demands
    :documentation "A list of keywords indicating which information the
evaluation function needs.  Possible choices depend on problem and
discretization, e.g. :local, :global, :solution are choices which will
probably be understood.")
   (eval :accessor coeff-eval :initarg :eval :type function
    :documentation "The evaluation funtion.  It accepts a list of
keyword parameters which should correspond to the list in DEMANDS.")
   (residual :initform t :initarg :residual
    :documentation "T means evaluation for computing the residual.")
   (jacobian :initform t :initarg :jacobian
    :documentation "T means evaluation for computing the Jacobian."))
  (:documentation "The coefficient class."))

(defmethod evaluate ((coeff <coefficient>) (input list))
  "The pairing between coefficient and input."
  (apply (slot-value coeff 'eval) input))

(defmethod demands ((coeffs list))
  "Returns unified demands for all coefficients in the list."
  (let ((demands nil))
    (loop for coeff in coeffs do
	  (loop for (demand value) on (demands coeff) by #'cddr do
		(aif (getf demands demand)
		     (unless (eql it value)
		       (error "coefficient demands contradict"))
		     (setf (getf demands demand) value))))
    demands))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; some trivial coefficient functions
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(declaim (inline constant-coefficient))
(defun constant-coefficient (value &rest other-values)
  "Returns a coefficient which takes the given value.  Several values can
be passed which is needed, for example, for returning also the type of a
boundary condition."
  (make-instance
   '<coefficient> :demands () :eval
   (if other-values
       #'(lambda (&key &allow-other-keys)
	   (apply #'values value other-values))
       (constantly value))))

(defun function->coefficient (func)
  "The function argument @arg{func} is transformed into a coefficient
depending on global coordinates."
  (make-instance '<coefficient> :demands '(:global)
		 :eval #'(lambda (&key global &allow-other-keys)
			   (funcall func global))))

(defun ensure-coefficient (obj)
  "Returns @arg{obj} if it is a coefficient, converts @arg{obj} into a
coefficient depending on the space variable if @arg{obj} is a function;
otherwise, @arg{obj} is made into a constant coefficient."
  (cond ((typep obj '<coefficient>) obj)
	((functionp obj) (function->coefficient obj))
	((typep obj '<function>)
	 (make-instance '<coefficient> :demands '(:global)
			:eval #'(lambda (&key global &allow-other-keys)
				  (evaluate obj global))))
	(t (constant-coefficient obj))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; <domain-problem>
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass <domain-problem> (<problem>)
  ((domain :reader domain :initform (required-argument)
	   :initarg :domain :type <domain>)
   (coefficients :accessor coefficients :initform (make-hash-table)
		 :initarg :coefficients :documentation
		 "Hash table which maps domain patches to coefficients.")
   (multiplicity :reader multiplicity :initform 1 :initarg :multiplicity))
  ;;
  (:documentation "Base-class for problems posed on a domain.  The slot
DOMAIN contains the domain on which the problem lives.  The slot
COEFFICIENTS contains a table from domain patches to coefficients on this
patch which are property lists of the form (SYM1 coefficient1 SYM2
coefficient2 ...).  When the problem instance is initialized this table is
set up by calling the function PATCH->COEFFICIENTS which has to be provided
as a key argument.  The multiplicity slot can be chosen as n>1 if the
problem is posed with n different right hand sides simultaneously."))

(defmethod initialize-instance ((problem <domain-problem>)
				&key patch->coefficients &allow-other-keys)
  "Setup the coefficient table, if the coefficients are given as a function
mapping domain patches to coefficient property lists."
  (call-next-method)
  (with-slots (domain coefficients) problem
    (when patch->coefficients
      (doskel (patch domain)
	(setf (gethash patch coefficients)
	      (funcall patch->coefficients patch))))))

(defmethod describe-object :after ((problem <domain-problem>) stream)
  (doskel ((patch properties) (domain problem))
    (format t "~&Cell ~A  [Mapping ~A]~%Properties: ~A~%Coeffs: ~A~2%"
	    patch (mapping patch) properties
	    (coefficients-of-patch patch problem))))

(declaim (inline coefficients-of-patch coefficients-of-cell))
(defun coefficients-of-patch (patch problem)
  "An accessor for the coefficients."
  (the list (gethash patch (slot-value problem 'coefficients))))

(defun coefficients-of-cell (cell mesh problem)
  "An accessor for the coefficients."
  (coefficients-of-patch (patch-of-cell cell mesh) problem))

(defgeneric interior-coefficients (problem)
  (:documentation "Yields a list of possible interior coefficients for PROBLEM.")
  (:method (problem) ()))

(defgeneric boundary-coefficients (problem)
  (:documentation "Yields a list of possible boundary coefficients for PROBLEM.")
  (:method (problem) "The following coefficients make sense for many pde
problems.  PERIODIC: periodic boundary conditions for non-identified
boundaries.  This is not yet implemented (not needed?).  CONSTRAINT:
essential boundary conditions."  '(PERIODIC CONSTRAINT)))

(defgeneric all-coefficients (problem)
  (:documentation "Yields a list of possible coefficients for @arg{problem}.")
  (:method (problem)
	   (append (interior-coefficients problem)
		   (boundary-coefficients problem))))

(defgeneric interior-coefficient-p (problem coeff)
  (:documentation "Tests, if @arg{coeff} is an interior coefficient of @arg{problem}.")
  (:method (problem coeff) (member coeff (interior-coefficients problem))))

(defgeneric boundary-coefficient-p (problem coeff)
  (:documentation "Tests, if @arg{coeff} is a boundary coefficient of @arg{problem}.")
  (:method (problem coeff) (member coeff (boundary-coefficients problem))))

(defgeneric coefficient-p (problem coeff)
  (:documentation "Test if @arg{coeff} is a coefficient of @arg{problem}.")
  (:method (problem coeff) (member coeff (coefficients problem))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; <interpolation-problem>
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass <interpolation-problem> (<domain-problem>)
  ()
  (:documentation "Interpolation problem on a domain.  The function which
is to be interpolated is given as a coefficient with key FUNCTION in the
coefficient list."))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; <pde-problem>
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass <pde-problem> (<domain-problem>)
  ()
  (:documentation "Base-class for a pde-problem."))

(defmethod initialize-instance
    ((problem <pde-problem>) &key (linear-p t linear-p-supplied)
     (coercive nil coercive-supplied) &allow-other-keys)
  "Set properties which are explicitly initialized."
  (call-next-method)
  (when linear-p-supplied
    (setf (get-property problem 'linear-p) linear-p))
  (when coercive-supplied
    (setf (get-property problem 'coercive) coercive)))

(defmethod initialize-instance :after ((problem <pde-problem>) &key &allow-other-keys)
  "Flag linearity and coercivity of the problem.  Linearity is determined
by looking at dependencies of the coefficient functions on the solution.
Coercivity is determined by looking for essential constraints on the ansatz
space."
  (unless (property-set-p problem 'linear-p)
    (setf (get-property problem 'linear-p) t)
    (dohash ((coeffs) (coefficients problem))
      (loop for (nil coeff) on coeffs by #'cddr
	    when (member :solution (demands coeff)) do
	    (setf (get-property problem 'linear-p) nil))))
  (unless (property-set-p problem 'coercive)
    (dohash ((coeffs) (coefficients problem))
      (when (member 'CONSTRAINT coeffs)
	(setf (get-property problem 'coercive) t)))))

(defgeneric dual-problem (problem functional)
  (:documentation "Returns the dual problem for @arg{problem} with the
right-hand side given by @arg{functional}.  The solution of this problem
measures the sensitivity of @arg{functional} applied to the solution of
problem with respect to errors in the solution."))

(defgeneric self-adjoint-p (problem)
  (:documentation "Returns two values.  The first says if @arg{problem} is
self-adjoint, the second says if that value has really been checked."))

(defmethod self-adjoint-p (problem)
  "The default method says that @arg{problem} is not self-adjoint and that no
check has been performed."
  (values nil nil))

;;; Testing: (test-pde-problem)

(defun test-pde-problem ()
  (function->coefficient
   #'(lambda (&key global &allow-other-keys) (exp (aref global 0))))
  (make-instance '<pde-problem> :domain *unit-interval-domain*)
  )

(fl.tests:adjoin-test 'test-pde-problem)