;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; problem.lisp
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

(in-package :problem)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; <problem>
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass <problem> ()
  ((domain :reader domain :initform (ext:required-argument)
	   :initarg :domain :type <domain>)
   (p->c :reader patch->coefficients :initform (constantly nil)
	 :initarg :patch->coefficients)
   (memoize :initform t :initarg :memoize)
   (multiplicity :reader multiplicity :initform 1 :initarg :multiplicity))
  ;;
  (:documentation "Base-class for a pde-problem.  The domain slot contains
the domain on which the problem lives.  The p->c slot contains a map from
the domain patches to problem coefficients.  Those are property lists of
the form (SYM1 coefficient1 SYM2 coefficient2 ...).  The multiplicity slot
can be chosen as n>1 if the problem is posed with n different right hand
sides simultaneously."))

(defmethod initialize-instance :after ((problem <problem>) &key &allow-other-keys)
  "Memoize the coefficient definition.  This might create problems for
strange applications with dynamically changing problems."
  (when (slot-value problem 'memoize)
    (setf (slot-value problem 'p->c)
	  (memoize-1 (slot-value problem 'p->c)))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Coefficient functions
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(definline coefficients-of-patch (patch problem)
  "An accessor for the coefficients."
  (the list (funcall (patch->coefficients problem) patch)))

(definline coefficients-of-cell (cell mesh problem)
  "An accessor for the coefficients."
  (coefficients-of-patch (patch-of-cell cell mesh) problem))

(defparameter *general-problem-keywords* '(PERIODIC CONSTRAINT)
  "This list contains keywords which make sense for several problems.
PERIODIC: periodic boundary conditions for non-identified boundaries.  This
is not yet implemented (not needed?).
CONSTRAINT: essential boundary conditions.")

(defgeneric interior-coefficients (problem)
  (:documentation "Yields a list of possible interior coefficients for problem."))

(defgeneric boundary-coefficients (problem)
  (:documentation "Yields a list of possible boundary coefficients for problem."))

(defgeneric coefficients (problem)
  (:documentation "Yields a list of possible coefficients for problem."))

(defgeneric interior-coefficient-p (problem coeff)
  (:documentation "Tests, if coeff is an interior coefficient of problem."))

(defgeneric boundary-coefficient-p (problem coeff)
  (:documentation "Tests, if coeff is a boundary coefficient of problem."))

(defgeneric coefficient-p (problem coeff)
  (:documentation "Tests, if coeff is a coefficient of problem."))

(defgeneric dual-problem (problem functional)
  (:documentation "Returns the dual problem for problem with the right-hand
side given by functional.  The solution of this problem measures the
sensitivity of functional applied to the solution of problem with respect
to errors in the solution."))

(defgeneric self-adjoint-p (problem)
  (:documentation "Returns two values.  The first says if the problem is
self-adjoint, the second says if that value has really been checked."))

(defmethod self-adjoint-p (problem)
  "Default method says that problem is not self-adjoint and that no check has been performed."
  (values nil nil))

(defmethod coefficients (problem)
  "Standard method: tests if coeff is interior or boundary coefficient."
  (append (interior-coefficients problem)
	  (boundary-coefficients problem)))

(defmethod interior-coefficients (problem) ())

(defmethod boundary-coefficients (problem) ())

(defmethod interior-coefficient-p (problem coeff)
  (member coeff (interior-coefficients problem)))

(defmethod boundary-coefficient-p (problem coeff)
  (member coeff (boundary-coefficients problem)))

(defmethod coefficient-p (problem coeff)
  "Standard method: tests if coeff is interior or boundary coefficient."
  (member coeff (coefficients problem)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; coefficient
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass <coefficient> ()
  ((demands :reader demands :initform () :initarg :demands
	    :documentation "A list of keywords indicating which information
the evaluation function needs.  Possible choices depend on problem and
discretization, e.g. :local, :global, :solution are choices which will
probably be understood.")
   (eval :accessor coeff-eval :initarg :eval :type function
	 :documentation "The evaluation funtion.  It accepts a list of
keyword parameters which should correspond to the list in @code{needed}.")
   (residual :initform t :initarg :residual)
   (jacobian :initform t :initarg :jacobian))
  (:documentation "The coefficient class."))

(defmethod evaluate ((coeff <coefficient>) (input list))
  "The pairing between coefficient and input."
  (apply (slot-value coeff 'eval) input))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; some trivial coefficient functions
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(declaim (inline constant-coefficient))
(defun constant-coefficient (value &rest other-values)
  (make-instance
   '<coefficient> :demands () :eval
   (if other-values
       #'(lambda (&rest args)
	   (declare (ignore args))
	   (apply #'values value other-values))
       (constantly value))))

(defparameter *cf-constantly-0.0d0* (constant-coefficient 0.0d0))
(defparameter *cf-constantly-1.0d0* (constant-coefficient 1.0d0))

(defun function->coefficient (func)
  "Returns a coefficient for the given function depending on global
coordinates."
  (make-instance '<coefficient> :demands '(:global)
		 :eval #'(lambda (&key global &allow-other-keys)
			   (funcall func global))))

(defun ensure-coefficient (obj)
  "Filter ensuring that the result is a coefficient function."
  (if (functionp obj) (function->coefficient obj) obj))

(defun problem-info (problem &optional mesh &key (where :all))
  "Displays information about the given problem."
  (doskel ((cell properties) (or mesh (domain problem)))
    (let ((patch (if mesh (patch-of-cell cell mesh) cell)))
      (when (or (eq where :all)
		(null mesh)
		(eq where (if (refined-p cell mesh) :refined :surface)))
	(format t "~&Cell ~A  [Mapping ~A]~%Properties: ~A~%Coeffs: ~A~2%"
		cell (mapping cell) properties
		(funcall (patch->coefficients problem) patch))))))

;;; Testing: (test-problem)

(defun test-problem ()
  (function->coefficient
   #'(lambda (&key global &allow-other-keys) (exp (aref global 0))))
  (make-instance '<problem> :domain *unit-interval-domain*
		 :patch->coefficients nil)
  )

(tests::adjoin-femlisp-test 'test-problem)