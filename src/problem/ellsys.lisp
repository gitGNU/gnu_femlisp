;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; ellsys.lisp
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;
;;; Copyright (C) 2006 Nicolas Neuss, University of Karlsruhe.
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
;;; NO EVENT SHALL THE AUTHOR, THE UNIVERSITY OF KARLSRUHE OR OTHER
;;; CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
;;; EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
;;; PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
;;; PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
;;; LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
;;; NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
;;; SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(in-package :cl-user)

(defpackage "FL.ELLSYS"
  (:use "COMMON-LISP" "FL.MACROS" "FL.UTILITIES" "FL.DEBUG"
	"FL.MATLISP" "FL.FUNCTION"
	"FL.MESH" "FL.PROBLEM")
  (:export "<ELLSYS-PROBLEM>" "NR-OF-COMPONENTS" "LINEARIZATION"
	   "ISOTROPIC-DIFFUSION" "ISOTROPIC-DIFFUSION-COEFFICIENT"
	   "DIAGONAL-REACTION-COEFFICIENT"
	   "ELLSYS-ONE-FORCE-COEFFICIENT"
	   "UNIT-VECTOR-FORCE-COEFFICIENT"
	   "ELLSYS-MODEL-PROBLEM")
  (:documentation "This package contains the problem definition of systems
of convection-diffusion-reaction equations.  The system is given in the
following form which is suited for a fixed-point iteration:

@math{-div(a(x,u_old,\nabla u_old) \nabla u)
 + div(b(x,u_old,\nabla u_old) u) +
 + c(x,u_old,\nabla u_old) \nabla u +
 + r(x,u_old,\nabla u_old) u
= f(x,u_old, \nabla u_old) 
- div(g(x,u_old, \nabla u_old))
- div(a(x,u_old,\nabla u_old) h(x,u_old, \nabla u_old)) }

where @math{u:G \to \R^N}.  Note that the last two terms are introduced in
the variational formulation and imply a natural Neumann boundary condition
@math{\derivative{u}{n} = (g+a h) \cdot n} at boundaries where no Dirichlet
constraints are posed."))

(in-package :fl.ellsys)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; <ellsys-problem>
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass <ellsys-problem> (<pde-problem>)
  ()
  (:documentation "Systems of convection-diffusion-reaction equations.  The
coefficients should be vector-valued functions in this case."))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Generation of standard ellsys problems
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun isotropic-diffusion (dim values)
  "Returns a sparse diagonal diffusion tensor with isotropic diffusion in
each component.  @arg{value} should be a vector or a number and contains
the amount of diffusion in each component."
  (let ((values (if (numberp values) (vector values) values)))
    (diagonal-sparse-tensor
     (map 'vector #'(lambda (val) (scal val (eye dim))) values))))

(defun isotropic-diffusion-coefficient (dim values)
  (constant-coefficient 'FL.ELLSYS::A (isotropic-diffusion dim values)))

(defun diagonal-reaction-coefficient (values)
  (constant-coefficient
   'FL.ELLSYS::R
    (diagonal-sparse-tensor
     (map 'vector #'(lambda (val) (if (numberp val) (scal val (ones 1)) val))
	  values))))

(defun ellsys-one-force-coefficient (nr-comps multiplicity)
  (constant-coefficient
   'FL.ELLSYS::F
   (make-array nr-comps :initial-element (ones 1 multiplicity))))

(defun unit-vector-force-coefficient (direction)
  (constant-coefficient
   'FL.ELLSYS::F
   (lret ((result (make-instance '<sparse-tensor> :rank 1)))
     (setf (tensor-ref result direction) (eye 1)))))

(defun ellsys-model-problem
    (domain components
     &key a b c d r f g h (dirichlet nil dirichlet-p)
     initial sigma evp properties derived-class)
  "Generates a rather general elliptic problem on the given domain."
  (when (numberp components)
    (setq components (list (list 'u components))))
  (setq domain (ensure-domain domain))
  (unless dirichlet-p
    (setq dirichlet
	  (let ((nr-comps (loop for comp in components summing
				(if (listp comp) (second comp) 1))))
	    (constraint-coefficient nr-comps 1))))
  (apply #'fl.amop:make-programmatic-instance
	 (remove nil (list (or derived-class '<ellsys-problem>)
			   (and initial sigma '<time-dependent-problem>)))
	 :components components
	 :properties properties
	 :domain domain :patch->coefficients
	 `((:external-boundary
	    ,(when dirichlet (list dirichlet)))
	   (:d-dimensional
	    ,(macrolet ((coefflist (&rest coeffs)
				   `(append
				     ,@(loop for sym in coeffs collect
					     `(when ,sym (list (ensure-coefficient ',sym ,sym)))))))
		       (coefflist a b c d r f g h initial sigma))))
	 (append
	  (when evp (destructuring-bind (&key lambda mu) evp
		      (list :lambda lambda :mu mu))))
	 ))

(with-memoization (:type :local :id :linearize :size 1)
  (defun linearization (reaction-jet)
    "Returns a list of two functions namely @math{u \mapsto R(u)-DR(u)u}
and @math{u \mapsto -DR(u)} which can be used directly in the
discretization as source and reaction term."
    (flet ((linearize (u)
	     (memoizing-let ((u u))
	       (dbg :disc "Evaluating reaction jet")
	       (destructuring-bind (R DR)
		   (funcall reaction-jet u)
		 (list (gemm -1.0 DR u +1.0 R)
		       (scal -1.0 DR))))))
      (list (compose #'first #'linearize)
	    (compose #'second #'linearize)))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; new problem definition support
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



;;; Testing

(defun test-ellsys ()
  (let ((problem
	 (ellsys-model-problem 
	  2 '(u v)
	  :a (isotropic-diffusion 2 #(1.0 2.0))
	  :b (vector #m((1.0) (0.0)) #m((0.0) (1.0)))
	  :f (vector #m(1.0) #m(1.0))
	  :dirichlet (constraint-coefficient 2 1))))
    (assert (components problem))
    problem)
  )

;;; (test-ellsys)
(fl.tests:adjoin-test 'test-ellsys)
