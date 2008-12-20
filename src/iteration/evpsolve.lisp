;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; evpsolve.lisp - Provide solvers for eigenvalue problems
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;
;;; Copyright (C) 2004- Nicolas Neuss, University of Heidelberg.
;;; Copyright (C) 2008- Nicolas Neuss, University of Karlsruhe.
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
;;; NO EVENT SHALL THE AUTHOR, THE UNIVERSITIES HEIDELBERG AND KARLSRUHE OR
;;; OTHER CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
;;; SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
;;; LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
;;; DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
;;; THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
;;; (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
;;; OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(in-package :fl.iteration)

(defclass <evp-solver> (<nonlinear-solver>)
  ()
  (:documentation "General class for a solver of eigenvalue problems."))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Wielandt vector iteration
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass <wielandt-iteration> (<evp-solver>)
  ()
  (:documentation "Wielandt iteration."))

(defmethod next-step ((wielandt <wielandt-iteration>) blackboard)
  "Simply calls the linear solver on the linearized problem."
  (with-items (&key problem solution residual residual-p
		    linearization linear-solver-failure step)
      blackboard
    (unless (slot-boundp wielandt 'linear-solver)
      (setf (slot-value wielandt 'linear-solver)
	    (?1 (lu-solver)
		(select-linear-solver linearization blackboard))))
    (copy! solution (rhs linearization))
    (handler-case
       (solve (linear-solver wielandt)
	      (blackboard :problem linearization :solution solution
			  :residual residual :residual-p residual-p))
     (arithmetic-error ()
       ;; In the case of an arithmetic error we have found an eigenvalue.
       ;; In practice, this will rarely occur.  Note that if this occurs in
       ;; step 1, the eigenvector may be completely wrong.
       (setq linear-solver-failure t)))
    (let ((mass (mass problem solution))
	  (energy (energy problem solution)))
      (with-slots (eigenvalues) problem
	(setf (aref eigenvalues 0) (/ energy mass))
	(scal! (/ (sqrt mass)) solution)
	(dbg :iter "mass=~A energy=~A eigenvalue=~A x=~A~%"
	     mass energy (aref eigenvalues 0) solution)))
    (setq residual-p nil)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; LOBPCG
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass <lobpcg> (<evp-solver>)
  ()
  (:documentation "Knyazev's LOBPCG iteration."))

(defun gram-schmidt (X B &key (tolerance 0.0))
  "Orthogonalizes the columns of X wrt to the bilinear form given by the
symmetric matrix B."
  ;; split X in vectors, do the orthogonalization, then join the vectors
  (let* ((n (ncols X))
	 (vecs (loop for k below n collect
		    (matrix-slice X :from-col k :ncols 1))))
    (flet ((bilinear (B u v)
	     (dot (m* B u) v)))
      (apply #'join :horizontal
	     (loop+ (i (u vecs))
		when
		(loop
		   ;; orthogonalize
		   (loop+ ((j (range :below i)) (v vecs)) do
		      (axpy! (- (bilinear B u v)) v u))
		   (let ((s (abs (bilinear B u u))))
		     (if (> s tolerance)
			 (return (scal! (/ (sqrt s)) u))
			 #+(or)(fill-random! u 1.0)
			 (return nil))))
		collect it)))))

(defun gram-schmidt-new (X B)
  "Unfortunately, this routine squares a bad condition number of X."
  (let* ((gram (m*-tn X (m* B X)))
	 (G (cholesky gram :inverse t :compactify t)))
    (let ((fl.matlisp::*print-matrix-pretty* t)
	  (fl.matlisp::*print-matrix* t))
    (dbg :lobpcg "gram=~%~A~%G=~%~A~%" gram G))
    (m* X G)))

(defmethod next-step ((lobpcg <lobpcg>) blackboard)
  (declare (optimize debug))
  (with-items (&key problem solution old-solution residual residual-p
		    linearization step)
      blackboard
    (with-slots (stiffness-matrix mass-matrix multiplicity) problem
      (unless (slot-boundp lobpcg 'linear-solver)
	(setf (slot-value lobpcg 'linear-solver)
	      (select-linear-solver linearization blackboard)))
      (dbg :lobpcg "A=~A, B=~B" stiffness-matrix mass-matrix)
      (let ((new-direction (copy solution)))
	(solve (linear-solver lobpcg)
	       (blackboard :problem linearization :solution new-direction
			   :residual residual :residual-p residual-p))
	(let* ((modified-new-direction
		(gram-schmidt new-direction mass-matrix))
	       (all-x
		(gram-schmidt
		 (apply #'join :horizontal
			(remove nil (list solution modified-new-direction
					  (and old-solution
					       (m- old-solution solution)))))
		 mass-matrix)))
	  (dbg-when :lobpcg
	    (format t "all-x=~%")
	    (let ((fl.matlisp::*print-matrix-pretty* t)
		  (fl.matlisp::*print-matrix* t))
	      (show all-x)))
	  ;; Rayleigh-Ritz
	  (let ((A (m*-tn all-x (m* stiffness-matrix all-x)))
		(B (m*-tn all-x (m* mass-matrix all-x))))
	    (let ((fl.matlisp::*print-matrix-pretty* t)
		  (fl.matlisp::*print-matrix* t))
	      (dbg :lobpcg "HEGV called with A=~%~A~%B=~A~%" A B))
	    (multiple-value-bind (mu Z)
		(hegv A B :V)
	      (assert (not (mzerop Z)))
	      (dbg :lobpcg "HEGV returned mu=~A" mu)
	      (shiftf old-solution solution
		      (m* all-x (matrix-slice Z :from-col 0 :ncols multiplicity)))
	      (shiftf (getbb blackboard :previous-eigenvalues)
		      (getbb blackboard :eigenvalues)
		      (copy! (subseq (store mu) 0 multiplicity)
			     (slot-value problem 'eigenvalues)))))))
      (setq residual-p nil))))

(defvar *evp-success*
  '(and (> :step 2)
    (>= :step-reduction 1.0)
    (>= (vector-last :eigenvalues) (vector-last :previous-eigenvalues))))

(defmethod select-solver ((problem <evp>) blackboard)
  (declare (ignore blackboard))
  (?2 (make-instance
       '<wielandt-iteration>
       :success-if (or (getbb blackboard :success-if)
		       '(and (> :step 2) (> :step-reduction 0.5)))
       :failure-if (or (getbb blackboard :failure-if)
		       '(and (> :step 1) (> :step-reduction 1.0) (> :defnorm 1.0e-5))))
      (make-instance
       '<lobpcg>
       :success-if *evp-success*)))

;;;; Testing

(defun test-evpsolve ()
  (let* ((n 3)
	 (evp (make-instance
	       '<ls-evp>
	       :stiffness-matrix (laplace-full-matrix n 1)
	       :mass-matrix (eye n)))
	 (bb (blackboard :problem evp :output t)))
    (solve bb)
    (slot-value evp 'eigenvalues))
  
  (let* ((n 20) (dim 2)
	 (evp (make-instance
	       '<ls-evp>
	       :stiffness-matrix (laplace-full-matrix n dim)
	       :mass-matrix (eye (expt n dim))
	       :multiplicity 10))
	 (bb (blackboard :problem evp :output 1)))
    (solve bb)
    (scal (expt (1+ n) 2) (slot-value evp 'eigenvalues)))
  )

;;;; Testing: (test-evpsolve)
(fl.tests:adjoin-test 'test-evpsolve)
