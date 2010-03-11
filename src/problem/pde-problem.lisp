;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; pde-problem.lisp
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;
;;; Copyright (C) 2003-2006 Nicolas Neuss, University of Heidelberg.
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

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; coefficient
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; Q: would it be reasonable to include the problem?
(defclass <coefficient> (property-mixin)
  ((name :accessor coefficient-name :initform nil :initarg :name)
   (dimension :reader dimension :initform nil :initarg :dimension :documentation
	      "The dimension of the cell on which this coefficient is
active.  The value T means that it is active on all cells lying on the
patch.  The default value NIL means that it is active on cells with the
same dimension as the patch.")
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

#+(or)
(defmethod initialize-instance :after ((coeff <coefficient>) &key &allow-other-keys)
  "For compatibility with older coefficient definitions."
  (mapf (slot-value coeff 'demands)
	(lambda (x) (if (atom x) (cons x 1) x))))

(defmethod print-object :after ((coeff <coefficient>) stream)
  "Prints the coefficient name."
  (format stream "{~A}" (coefficient-name coeff)))

(defmethod evaluate ((coeff <coefficient>) (input list))
  "The pairing between coefficient and input."
  (apply (slot-value coeff 'eval) input))

(defun get-coefficient (coeffs name)
  "Get coefficient @arg{name} from the list @arg{coeffs}."
  (find name coeffs :key #'coefficient-name))

(defmethod demands ((coeffs list))
  "Returns unified demands for all coefficients in the list."
  (let ((demands nil))
    (dolist (coeff coeffs)
      (dolist (demand (demands coeff))
	(pushnew demand demands :test #'equalp)))
    demands))

(defun copy-coefficient (coeff &key demands name eval)
  (with-slots (dimension residual jacobian) coeff
    (make-instance
     '<coefficient>
     :name (or name (coefficient-name coeff))
     :dimension dimension
     :demands (or demands (demands coeff))
     :eval (or eval (slot-value coeff 'eval))
     :residual residual :jacobian jacobian)))
    
(defun constant-coefficient (name value &rest other-values)
  "Returns a coefficient which takes the given value.  Several values can
be passed which is needed, for example, for returning also the type of a
boundary condition."
  (make-instance
   '<coefficient> :name name :demands () :eval
   (if other-values
       #'(lambda (&key &allow-other-keys)
	   (apply #'values value other-values))
       (constantly value))))

(defun constraint-coefficient (components multiplicity)
  "Returns a coefficient function which sets Dirichlet zero boundary
conditions for all components of a PDE system."
  (make-instance
   '<coefficient> :name 'CONSTRAINT :dimension t :demands ()
   :eval (lambda (&key &allow-other-keys)
	   (values
	    (make-array components :initial-element t)
	    (make-array components :initial-element (zeros 1 multiplicity))))))

(defgeneric zero-constraints (problem)
  (:documentation "Returns a coefficient function which constrains all
system components to zero.")
  (:method (problem)
    (constraint-coefficient
     (nr-of-components problem) (multiplicity problem))))

(defun identification-coefficient (master mapping)
 "A special coefficient used for identifying parts of the domain.  The
  coefficient evaluation returns the master coordinates."
 (lret ((coeff (make-instance
                '<coefficient> :name 'IDENTIFICATION
                :dimension t :demands '(:local :global)
                :eval mapping :residual nil :jacobian nil)))
   (setf (get-property coeff 'MASTER) master)))

(defun identification-coefficient-p (coeff)
  (and (eq (coefficient-name coeff) 'IDENTIFICATION)
       (get-property coeff 'MASTER)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun constraint-p (coeff)
  (member (coefficient-name coeff) '(CONSTRAINT IDENTIFICATION)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun sigma-p (coeff)
  (eq (coefficient-name coeff) 'SIGMA))

;;; fe-parameters

(concept-documentation
 "fe-parameters is an element in the demand list of coefficient functions
which is a pair consisting of :fe-parameters followed by a list of symbols
denoting the ansatz-space functions on the blackboard to be evaluated at
each integration point.  Each entry is of the form (symbol . k). which
means that a k-jet has to be evaluated.

Examples:
@code{(:fe-parameters (:solution . 1))}
@code{(:fe-parameters (:solution . 1) (:flow-field . 0))}")
 
(defun fe-parameter-p (demand)
  (and (consp demand)
       (eq :fe-parameters (car demand))))

(defun parameter-name (para)
  (etypecase para
    (symbol para)
    (list (car para))))

(defun parameter-order (para)
  (etypecase para
    (symbol 0)
    (list (cdr para))))

(defun compress-fe-parameters (paras)
  (remove-duplicates
   (loop for para in paras collect
	 (cons (parameter-name para)
	       (loop for para2 in paras
		     when (eq (parameter-name para2) (parameter-name para))
		     maximize (parameter-order para2))))
   :test #'equalp))

(defun required-fe-functions (coeffs)
  "Returns a list of finite element functions required by the coefficients
in the property list @arg{coeffs}."
  (compress-fe-parameters
   (loop for coeff in coeffs appending
	 (loop for demand in (demands coeff)
	       when (fe-parameter-p demand)
	       appending (cdr demand)))))

(defun add-fe-parameters-demand (demands new-paras)
  "Adds a list of fe-functions to the demands."
  (let ((old-paras (cdr (find-if #'fe-parameter-p demands))))
    (cons (cons :fe-parameters
		(compress-fe-parameters (union old-paras new-paras)))
	  (remove-if #'fe-parameter-p demands))))

(defun solution-dependent (coeffs)
  (find :solution (required-fe-functions coeffs) :key #'car))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; some coefficient functions
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun fx->coefficient (name func)
  "The function argument @arg{func} is transformed into a coefficient
depending on global coordinates."
  (make-instance '<coefficient> :name name :demands '(:global)
		 :eval #'(lambda (&key global &allow-other-keys)
			   (funcall func global))))

(defun fu->coefficient (name func)
  "The function argument @arg{func} is transformed into a coefficient
depending on the solution."
  (make-instance '<coefficient> :name name :demands '((:fe-parameters :solution))
		 :eval #'(lambda (&key solution &allow-other-keys)
			   (funcall func solution))))

(defun fxu->coefficient (name func)
  "The function argument @arg{func} is transformed into a coefficient
depending on position and solution."
  (make-instance '<coefficient> :name name
		 :demands '(:global (:fe-parameters :solution))
		 :eval #'(lambda (&key global solution &allow-other-keys)
			   (funcall func global solution))))

(defun function->coefficient (name func)
  (fx->coefficient name func))

(defun ensure-coefficient (name obj)
  "Returns @arg{obj} if it is a coefficient, converts @arg{obj} into a
coefficient depending on the space variable if @arg{obj} is a function;
otherwise, @arg{obj} is made into a constant coefficient."
  (cond ((typep obj '<coefficient>) obj)
	((functionp obj) (fx->coefficient name obj))
	((typep obj '<function>)
	 (make-instance '<coefficient> :name name :demands '(:global)
			:eval #'(lambda (&key global &allow-other-keys)
				  (evaluate obj global))))
	(t (constant-coefficient name obj))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; <domain-problem>
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass <domain-problem> (<problem>)
  ((domain :reader domain :initarg :domain :type <domain>)
   (components :reader components :initarg :components
	       :documentation "A list whose elements are lists of the form
\(symbol dim) or \(symbol dim type) or \(symbol subcomponents) describing
the components occuring in the pde.  Alternatively, this slot can contain a
function/dictionary mapping patches to such lists.")
   (multiplicity :reader multiplicity :initform 1 :initarg :multiplicity
		 :documentation "Multiplicity of the right-hand side.")
   (coefficients :accessor coefficients :initarg :coefficients :documentation
		 "A function mapping patches to coefficient lists."))
  ;;
  (:documentation "An instance of this class describes a problem posed on
@slot{domain} with coefficients given on each patch of the domain.  The
slot @slot{multiplicity} is a positive integer which denotes the number of
right-hand sides and solutions (e.g. when computing several eigenvectors at
once)."))

(defun component-symbol (comp)
  (if (symbolp comp) comp (first comp)))

(defun component-name (comp)
  (symbol-name (component-symbol comp)))

(defun component-length (comp)
  (cond ((symbolp comp) 1)
	((numberp (second comp)) (second comp))
	(t (count-components (second comp)))))

(defun count-components (comps)
  "Counts the total number of components."
  (reduce #'+ comps :key #'component-length))

(defun check-components (comps)
  "Checks if comps is of the correct form, see the documentation of the
slot @slot{components} in @class{domain-problem}."
  (every (lambda (comp)
	   (or (symbolp comp)
	       (and (listp comp)
		    (or (numberp (second comp))
			(check-components (second comp))))))
	  comps))

(defun extraction-information (components component)
  "If @arg{component} is in @arg{components}, a triple consisting of
position, length, and a flag is returned.  The flag is true, if the
component is a scalar."
  (loop
   for comp in components
   and off = 0 then (+ off (component-length comp))
   when (eq (component-symbol comp) component)
   do (return (values off (component-length comp) (symbolp comp)))
   when (and (listp comp) (listp (second comp)))
   do (multiple-value-bind (from length flag)
	  (extraction-information (second comp) component)
	(when from
	  (return (values from length flag))))))

(defun component-position (components comp)
  "Translates a symbol denoting a component to a position."
  (etypecase comp
    (number (values comp 1))
    (symbol
     (multiple-value-bind (from length)
	 (extraction-information components comp)
       (if from
	   (values from length)
	   (let* ((name (symbol-name comp))
		  (pos (position-if #'alpha-char-p name :from-end t)))
	     (when pos
	       (multiple-value-bind (from length)
		   (extraction-information
		    components (intern (subseq name 0 (1+ pos))
                                       (symbol-package comp)))
		 (when from
		   (let ((n (parse-integer (subseq name (1+ pos)))))
		     (when (<= 1 n length)
		       (values (+ from n -1) 1))))))))))))

				  
(defmethod initialize-instance :after ((problem <domain-problem>)
				       &key patch->coefficients &allow-other-keys)
  "Checks @slot{components}, does @slot{coefficients} setup when the
problem specification is given as a list in @arg{patch->coefficients}, and
finally memoizes @slot{coefficients}."
  (with-slots (domain coefficients) problem
    ;; component name check
    (doskel (patch domain)
      (let ((components (components-of-patch patch problem)))
	(assert (check-components components))
	(when (some (lambda (comp)
		      (let ((name (component-name comp)))
			(or (string-equal name "X")
			    (string-equal name "TIME")
			    (string-equal (subseq name 0 1) "D"))))
		    components)
	  (error "Illegal component names."))))
    ;; setup coefficients slot from old-style patch->coefficients
    (when patch->coefficients
      (assert (listp patch->coefficients))
      (setf coefficients
	    (lambda (patch)
	      (loop with classifications = (patch-classification patch domain)
		 for (id coeffs) in patch->coefficients
		 when (test-condition id classifications)
		 do (return coeffs)))))
    ;; memoize coefficients
    (setf coefficients
	  (if (slot-boundp problem 'coefficients)
	      (memoize-1 coefficients)
	      (constantly ())))))

(defmethod domain-dimension ((problem <domain-problem>))
  (dimension (domain problem)))

(defmethod describe-object :after ((problem <domain-problem>) stream)
  (doskel ((patch properties) (domain problem))
    (format stream "~&Cell ~A  [Mapping ~A]~%Properties: ~S~%Coeffs: ~A~2%"
	    patch (mapping patch) properties
	    (coefficients-of-patch patch problem))))

(defgeneric nr-of-components (problem)
  (:documentation "Returns the number of components for @arg{problem}.")
  (:method ((problem <domain-problem>))
    "Counts the number of components."
    (with-slots (components) problem
      (typecase components
	(list (reduce #'+ components :key #'component-length))
	(t nil)  ; no clear definition of number of components
	))))

(defun components-of-patch (patch problem)
  "Reader for the components of @arg{problem} on @arg{patch}."
  (with-slots (components) problem
    (etypecase components
      (list components)
      (function (funcall components patch))
      (hash-table (gethash patch components)))))

(defun (setf components-of-patch) (value patch problem)
  "Writer for the components of @arg{problem} on @arg{patch}."
  (with-slots (components) problem
    (ensure components (make-hash-table))
    (setf (gethash patch components) value)))

(defgeneric components-of-cell (cell mesh problem)
  (:documentation
   "Returns the components of @arg{problem} on @arg{cell}.")
  (:method (cell mesh problem)
    (components-of-patch (patch-of-cell cell mesh) problem)))

(defun coefficients-of-patch (patch problem)
  "Reader for the coefficients of @arg{patch} for @arg{problem}."
  (funcall (coefficients problem) patch))

(defgeneric coefficients-of-cell (cell mesh problem)
  (:documentation
   "Returns the coefficients of @arg{problem} on @arg{cell}.")
  (:method (cell mesh problem)
    "Default method, returns the coefficients of the associated patch."
    (coefficients-of-patch (patch-of-cell cell mesh) problem)))

(defun filter-applicable-coefficients (coeffs cell patch &key (constraints t))
  "Filters out the applicable coefficients for the respective cell with the
given patch."
  (remove-if-not
   (lambda (coeff)
     (and (or (eq (dimension coeff) t)
	      (= (dimension cell)
		 (or (dimension coeff) (dimension patch))))
	  (or constraints (not (constraint-p coeff)))))
   coeffs))

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

(defun classify-problem (problem)
  "Classifies a problem.  The result is written in its property list.
Unfortunately, this is rater crude at the moment.  Ideally, it could also
classify some cases of coerciveness automatically."
  (with-slots (linear-p coefficients) problem
    (setf linear-p t)
    (doskel (patch (domain problem))
      (when (solution-dependent (coefficients-of-patch patch problem))
	(setf linear-p nil)
	(return)))))

(defmethod initialize-instance :after
    ((problem <pde-problem>)
     &key (linear-p t linear-p-supplied) (coercive nil coercive-supplied)
     (classify t) &allow-other-keys)
  (when linear-p-supplied
    (setf (get-property problem 'linear-p) linear-p))
  (when coercive-supplied
    (setf (get-property problem 'coercive) coercive))
  ;; automatic classification
  (when classify (classify-problem problem)))

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

;;; Testing

(defun test-pde-problem ()
  (function->coefficient
   'test
   #'(lambda (&key global &allow-other-keys) (exp (aref global 0))))
  (make-instance '<pde-problem> :domain (n-cube-domain 1) :components '(u))
  )

;;; (test-pde-problem)
(fl.tests:adjoin-test 'test-pde-problem)