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
    :documentation "A list indicating which information the evaluation
function needs.  Possible choices depend on problem and discretization,
e.g. @code{:local}, @code{:global}, @code{:fe}, @code{:cell} are possible
choices.  One element can also be a list starting with the keyword
@code{:fe-parameters} and followed by symbols indicating names of finite
element functions on the discretization blackboard.")
   (eval :accessor coeff-eval :initarg :eval :type function
    :documentation "The evaluation funtion.  It accepts a list of
keyword parameters which should correspond to the list in DEMANDS.")
   (residual :initform t :initarg :residual
    :documentation "T means evaluation for computing the residual.")
   (jacobian :initform t :initarg :jacobian
    :documentation "T means evaluation for computing the Jacobian."))
  (:documentation "The coefficient class."))

(defun solution-dependent (coeff)
  "Tests if the coefficient is solution-dependent, i.e. if the
corresponding problem is nonlinear."
  (member :solution
	  (cdr (find-if #'(lambda (demand)
			    (and (consp demand)
				 (eq :fe-parameters (car demand))))
			(demands coeff)))))

(defmethod evaluate ((coeff <coefficient>) (input list))
  "The pairing between coefficient and input."
  (apply (slot-value coeff 'eval) input))

(defun get-coefficient (coeffs name)
  "Get coefficient @arg{name} from the list @arg{coeffs}."
  (getf coeffs name))

(defmethod demands ((coeffs list))
  "Returns unified demands for all coefficients in the list."
  (let ((demands nil))
    (dolist (coeff coeffs)
      (dolist (demand (demands coeff))
	(pushnew demand demands :test #'equalp)))
    demands))

(defun filter-applicable-coefficients (coeffs cell patch &key (constraints t))
  "Filters out the applicable coefficients for the respective cell with the
given patch."
  (loop
     for (symbol coeff) on coeffs by #'cddr
     when (and coeff
	       (or (eq (dimension coeff) t)
		   (= (dimension cell)
		      (or (dimension coeff) (dimension patch))))
	       (or constraints
		   (not (equal (symbol-name symbol) "CONSTRAINT"))))
     collect symbol and collect coeff))

(defun fe-parameter-p (demand)
  (and (consp demand)
       (eq :fe-parameters (car demand))))

(defun required-fe-functions (coeffs)
  "Returns a list of finite element functions required by the coefficients
in the property list @arg{coeffs}."
  (let ((fe-functions ()))
    (loop for (symbol coeff) on coeffs by #'cddr do
	  (dolist (demand (demands coeff))
	    (when (fe-parameter-p demand)
	      (_f union fe-functions (cdr demand)))))
    fe-functions))

(defun add-fe-parameters-demand (demands new-paras)
  "Adds a list of fe-functions to the demands."
  (let ((old-paras (cdr (find-if #'fe-parameter-p demands))))
    (cons (cons :fe-parameters (union old-paras new-paras))
	  (remove-if #'fe-parameter-p demands))))
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; some coefficient functions
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

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

(defun f[x]->coefficient (func)
  "The function argument @arg{func} is transformed into a coefficient
depending on global coordinates."
  (make-instance '<coefficient> :demands '(:global)
		 :eval #'(lambda (&key global &allow-other-keys)
			   (funcall func global))))

(defun f[u]->coefficient (func)
  "The function argument @arg{func} is transformed into a coefficient
depending on the solution."
  (make-instance '<coefficient> :demands '((:fe-parameters :solution))
		 :eval #'(lambda (&key solution &allow-other-keys)
			   (funcall func solution))))

(defun f[xu]->coefficient (func)
  "The function argument @arg{func} is transformed into a coefficient
depending on position and solution."
  (make-instance '<coefficient> :demands '(:global (:fe-parameters :solution))
		 :eval #'(lambda (&key global solution &allow-other-keys)
			   (funcall func global solution))))

(defun function->coefficient (func)
  (f[x]->coefficient func))

(defun ensure-coefficient (obj)
  "Returns @arg{obj} if it is a coefficient, converts @arg{obj} into a
coefficient depending on the space variable if @arg{obj} is a function;
otherwise, @arg{obj} is made into a constant coefficient."
  (cond ((typep obj '<coefficient>) obj)
	((functionp obj) (f[x]->coefficient obj))
	((typep obj '<function>)
	 (make-instance '<coefficient> :demands '(:global)
			:eval #'(lambda (&key global &allow-other-keys)
				  (evaluate obj global))))
	(t (constant-coefficient obj))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; <domain-problem>
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass <domain-problem> (<problem>)
  ((domain :reader domain :initarg :domain :type <domain>)
   (coefficients :accessor coefficients :initform (make-hash-table)
		 :initarg :coefficients :documentation
		 "Hash table which maps domain patches to coefficients.")
   (multiplicity :reader multiplicity :initform 1 :initarg :multiplicity))
  ;;
  (:documentation "An instance of this class describes a problem posed on
the domain @slot{domain}.  The slot @slot{coefficients} contains a table
mapping domain patches to property lists of the form (@symbol{identifier1}
@code{coefficient1} @symbol{identifier2} @code{coefficient2} ...).  Here
identifiers are special symbols and the coefficients are objects of type
@class{<coefficient>}.  When the problem instance is initialized this table
is usually set up by calling the function @function{patch->coefficients}
which has to be provided as a key argument.  The slot @slot{multiplicity}
can be chosen as a positive integer @math{n} if the problem is posed with
@math{n} different right hand sides simultaneously."))

(defmethod initialize-instance ((problem <domain-problem>)
				&key patch->coefficients &allow-other-keys)
  "Setup the coefficient table, if the coefficients are given as a function
mapping domain patches to coefficient property lists.  Instead of a
function, this mapping can also be given as a list describing the
association of patch classifications to coefficient functions."
  (call-next-method)
  (with-slots (domain coefficients) problem
    (when patch->coefficients
      (doskel (patch domain)
	(setf (gethash patch coefficients)
	      (typecase patch->coefficients
		(function (funcall patch->coefficients patch))
		(list (loop with classifications = (patch-classification patch domain)
			    for (id coeffs) in patch->coefficients
			    when (subsetp (if (consp id) id (list id)) classifications)
			    do (return coeffs)))
		(t (error "Unknown mapping."))))))))

(defmethod domain-dimension ((problem <domain-problem>))
  (dimension (domain problem)))

(defmethod describe-object :after ((problem <domain-problem>) stream)
  (doskel ((patch properties) (domain problem))
    (format stream "~&Cell ~A  [Mapping ~A]~%Properties: ~A~%Coeffs: ~A~2%"
	    patch (mapping patch) properties
	    (coefficients-of-patch patch problem))))

(defgeneric nr-of-components (problem)
  (:documentation "Returns the number of components for @arg{problem}."))

(defun coefficients-of-patch (patch problem)
  "An accessor for the coefficients of @arg{patch} for @arg{problem}."
  (the list (gethash patch (slot-value problem 'coefficients))))

(defgeneric coefficients-of-cell (cell mesh problem)
  (:documentation "An accessor for the coefficients of @arg{problem} valid
for @arg{cell}.")
  (:method (cell mesh problem)
    "This default method returns the coefficients of the associated patch."
    (coefficients-of-patch (patch-of-cell cell mesh) problem)))

(defgeneric interior-coefficients (problem)
  (:documentation "Returns a list of possible interior coefficients for
@arg{problem}.")
  (:method (problem)
    "Default method returns no coefficients."
    (declare (ignore problem))
    ()))

(defgeneric constraint-identifier (problem &optional cell)
  (:documentation
   "Returns the symbol which identifies constraints for this problem.")
  (:method (problem &optional cell)
    "Returns the CONSTRAINT symbol in the package of the class name."
    (declare (ignore cell))
    (find-symbol "CONSTRAINT"
		 (symbol-package (class-name (class-of problem))))))

(defgeneric boundary-coefficients (problem)
  (:documentation "Returns a list of possible boundary coefficients for
@arg{problem}.")
  (:method (problem)
    "This default method returns the constraint identifier for
@arg{problem}."
    (list (constraint-identifier problem))))

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

(defun constraint-coefficient (components multiplicity)
  "Returns a coefficient function which sets Dirichlet zero boundary
conditions for all components of a PDE system."
  (constant-coefficient
   (make-array components :initial-element t)
   (make-array components :initial-element (zeros 1 multiplicity))))

(defgeneric zero-constraints (problem)
  (:documentation "Returns a coefficient function which constrains all
system components to zero.")
  (:method (problem)
    (constraint-coefficient
     (nr-of-components problem) (multiplicity problem))))

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
      (when (member :solution (required-fe-functions coeffs))
	(setf (get-property problem 'linear-p) nil))))
  ;; very coarse test which will be often false for systems
  (let ((constraint-id (constraint-identifier problem)))
    (unless (property-set-p problem 'coercive)
      (dohash ((coeffs) (coefficients problem))
	(when (member constraint-id coeffs)
	  (setf (get-property problem 'coercive) t))))))

(defgeneric dual-problem (problem functional)
  (:documentation "Returns the dual problem for @arg{problem} with the
right-hand side given by @arg{functional}.  The solution of this problem
measures the sensitivity of @arg{functional} applied to the solution of
problem with respect to errors in the solution."))

(defgeneric self-adjoint-p (problem)
  (:documentation "Returns two values.  The first says if @arg{problem} is
self-adjoint, the second says if that value has really been checked.")
  (:method (problem)
    "The default method says that @arg{problem} is not self-adjoint and
that no check has been performed."
    (declare (ignore problem))
    (values nil nil)))

;;; Testing: (test-pde-problem)

(defun test-pde-problem ()
  (function->coefficient
   #'(lambda (&key global &allow-other-keys) (exp (aref global 0))))
  (make-instance '<pde-problem> :domain (n-cube-domain 1))
  )

(fl.tests:adjoin-test 'test-pde-problem)