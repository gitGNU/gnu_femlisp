;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; linsolve.lisp
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

(in-package :iteration)

;;; This file customizes the iterative solving strategy defined in
;;; iterate.lisp and solve.lisp to the solving of linear systems.

(defclass <linear-solver> (<iterative-solver>)
  ()
  (:documentation "Class for linear iterative solvers."))

(defmethod initial-guess ((linsolve <linear-solver>) blackboard)
  "Ensure an initial guess and its residual."
  (with-items (&key solution matrix rhs residual residual-p) blackboard
    (unless solution
      (setq solution (make-domain-vector-for matrix (multiplicity rhs))))))

(defmethod ensure-residual ((linsolve <linear-solver>) blackboard)
  "Ensure the residual for a linear problem."
  (with-items (&key solution matrix rhs residual residual-p) blackboard
    (unless residual-p  ; ensure residual
      (unless residual
	(setq residual (make-image-vector-for matrix (multiplicity rhs))))
      (setq residual (compute-residual matrix solution rhs residual))
      (setq residual-p t))))

(defun linsolve (mat rhs &key sol res output iteration (residual-norm #'norm)
		 (threshold 1.0d-12) reduction (maxsteps 100)
		 success-if failure-if &allow-other-keys)
  "In easy situations, this function provides a simpler interface for
solving a linear system."
  (let ((result
	 (solve (make-instance
		 '<linear-solver> :iteration iteration :output output
		 :residual-norm residual-norm :success-if
		 (or success-if
		     (cons 'or
			   (remove nil
				   (list 
				    (and threshold `(< :defnorm ,threshold))
				    (and reduction `(< :reduction ,reduction))))))
		 :failure-if (or failure-if (and maxsteps `(> :step ,maxsteps))))
		(blackboard :matrix mat :rhs rhs :solution sol :residual res))))
    (values (getbb result :solution) result)))

(defparameter *lu-solver*
  (make-instance
   '<linear-solver> :iteration *lu-iteration*
   :success-if '(>= :step 1) :output nil)
  "LU decomposition without pivoting.")

(defclass <special-solver> (<iterative-solver>)
  ((solver-function :reader solver-function :initarg :solver-function :type function))
  (:documentation "If you happen to have a problem-adapted solver given as
a function, you may use this base class."))

(defmethod solve ((solver <special-solver>) &optional blackboard)
  (funcall (solver-function solver) blackboard))

;;; solver-iteration
(defclass <solver-iteration> (<linear-iteration>)
  ((solver :initarg :solver)))

(defmethod make-iterator ((solve-it <solver-iteration>) mat)
  (make-instance
   '<iterator>
   :matrix mat
   :residual-before nil
   :initialize nil
   :iterate
   #'(lambda (x b r)
       (getbb (solve (slot-value solve-it 'solver)
		     (blackboard :matrix mat :rhs b :solution x :residual r))
	      :solution))
   :residual-after nil))

;;; Problem-dependent solving

(defgeneric select-linear-solver (object blackboard)
  (:documentation "Selects a linear solver for OBJECT.  OBJECT is usually a
matrix or a linear problem with certain characteristics."))

(defmethod select-linear-solver :around (object blackboard)
  "If a solver is on the blackboard, use it."
  (let ((solver (aif (getf blackboard :solver) it (call-next-method))))
    (setf (slot-value solver 'output)
	  (getbb blackboard :output))
    solver))

(defmethod select-linear-solver (object blackboard)
  "Default method selects LU decomposition."
  *lu-solver*)

;;; Not yet used: a general linear solver
(defclass <gps-linear-solver> (<linear-solver>)
  ((subsolver :type <linear-solver>))
  (:documentation "Problem-dependent linear solving."))

(defmethod solve ((gps-ls <gps-linear-solver>) &optional blackboard)
  (unless (slot-boundp gps-ls 'subsolver)
    (setf (slot-value gps-ls 'subsolver)
	  (select-linear-solver (getbb blackboard :matrix) blackboard)))
  (solve (slot-value gps-ls 'subsolver) blackboard))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Testing: (iteration::test-linsolve)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun test-linsolve ()
  "Linear solving on different types of vector/matrix representations."
  (let ((lu *lu-iteration*)
	(jac *undamped-jacobi*)
	(gs *gauss-seidel*))
    
    ;; application to matlisp matrices
    (let ((A #m((2.0 -1.0  0.0)
		(-1.0  2.0 -1.0)
		(0.0 -1.0  2.0)))
	  (b #m((1.0) (2.0) (1.0))))
      (linsolve A b :output t :iteration lu)
      (linsolve A b :output t :iteration jac)
      (linsolve A b :output t :iteration gs)
      (linsolve A b :output t :iteration (make-instance '<sor> :omega 1.18))
      ;; new interface
      (solve *lu-solver* (blackboard :matrix A :rhs b))
      (solve (make-instance '<linear-solver> :iteration jac :output t
			    :success-if '(> :step 10))
	     (blackboard :matrix A :rhs b)))
    
    ;; application to sparse matrices
    (let* ((constantly-1 (constantly 1))
	   (b (make-instance '<sparse-vector> :key->size constantly-1))
	   (A (make-sparse-matrix
	       :row-key->size constantly-1
	       :col-key->size constantly-1
	       :keys->pattern (constantly (full-crs-pattern 1 1))))
	   (sol (make-instance '<sparse-vector> :key->size constantly-1)))
      (setf (vref (vref b 1) 0) 1.0)
      (setf (vref (vref b 2) 0) 2.0)
      (setf (vref (vref b 3) 0) 1.0)
      
      (setf (vref (mref A 1 1) 0) 2.0)
      (setf (vref (mref A 2 2) 0) 2.0)
      (setf (vref (mref A 3 3) 0) 2.0)
      (setf (vref (mref A 1 2) 0) -1.0)
      (setf (vref (mref A 2 1) 0) -1.0)
      (setf (vref (mref A 2 3) 0) -1.0)
      (setf (vref (mref A 3 2) 0) -1.0)
      (assert (eql (mref (mref A 3 2) 0 0) -1.0))
      (show (linsolve A b :output t :iteration lu))
      (show (linsolve A b :output t :iteration gs))
      (solve *lu-solver* (blackboard :matrix A :rhs b :solution sol)))
    )
  )

(fl.tests:adjoin-test 'test-linsolve)

