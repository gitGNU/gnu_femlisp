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

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Linear solving
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defparameter *linsolve-output*
  (list :initial
	#'(lambda (&key indentation iteration &allow-other-keys)
	    (indented-format t "Linear iteration ~A" indentation
			     (class-name (class-of iteration)))
	    (indented-format t "Step (k)  |residual|  res[k]/res[k-1]" indentation)
	    (indented-format t "~37,,,'-<~>~%" indentation))
	:after-step
	#'(lambda (i defnorm previous-defnorm &key indentation &allow-other-keys)
	    (indented-format t "~4D     ~12,5,2E" indentation i defnorm)
	    (when previous-defnorm
	      (if (zerop previous-defnorm)
		  (format t "~@12A" (if (zerop defnorm) "nan" "infinity"))
		  (format t "  ~12,5,2E" (/ defnorm previous-defnorm))))
	    (terpri))
	:final nil))

(defun safe-divide-by-zero (a b)
  (if (zerop a)
      0.0d0
      (if (zerop b)
	  :infinity
	  (/ a b))))

;;; Old interface, better use the new one below

(defun linsolve (mat rhs &rest assembly-line
		 &key sol res output iteration (residual-norm #'norm)
		 (threshold 1.0d-12) threshabs reduction (maxsteps 100)
		 success-if failure-if &allow-other-keys)
  "This function solves the linear equation A sol = rhs by using the
linear-iteration on the matrix mat."
  (let ((indentation *output-indentation*)
	(*output-indentation* (+ *output-indentation* 5))
	(threshabs (and threshabs (if (functionp threshabs)
				      (apply threshabs assembly-line)
				      threshabs))))
    (awhen (and output (getf *linsolve-output* :initial))
      (funcall it :indentation indentation :iteration iteration))
    (setq res (or res (make-column-vector-for mat (multiplicity rhs))))
    (setq sol (or sol (make-row-vector-for mat (multiplicity rhs))))
    (dbg :iter "linsolve: making iterator")
    (with-slots (initialize iterate residual-after)
      (make-iterator iteration mat)
      (dbg :iter "linsolve: starting loop")
      (loop with initial-residual = (compute-residual mat sol rhs res)
	    with defnorm0 = (funcall residual-norm initial-residual)
	    initially (when initialize (funcall initialize sol rhs res))
	    for i from 0
	    for new-sol = sol then (progn (funcall iterate new-sol rhs residual) new-sol)
	    for residual = initial-residual then
	    (cond ((eq residual-after t) residual)
		  ((null residual-after) (compute-residual mat new-sol rhs residual))
		  (t (funcall residual-after new-sol rhs residual)))
	    for defnorm = defnorm0 then (funcall residual-norm residual)
	    and previous-defnorm = nil then defnorm
	    for red-factor = (if (zerop i)
				 :undefined
				 (safe-divide-by-zero defnorm defnorm0))
	    for success-p =
	    (if success-if
		(test-condition
		 success-if :defnorm defnorm :previous-defnorm previous-defnorm
		 :initial-defnorm defnorm0 :step i)
		(or (and threshold (< defnorm threshold))
		    (and (and (or (not threshabs)
				  (< defnorm threshabs))
			      reduction)
			 (not (eq reduction :ignore))
			 (numberp red-factor)
			 (< red-factor reduction))))
	    do
	    (awhen (and output (getf *linsolve-output* :after-step))
	      (funcall it i defnorm previous-defnorm :indentation indentation))
	    until (or success-p
		      (if failure-if
			  (test-condition
			   failure-if :defnorm defnorm :previous-defnorm previous-defnorm
			   :initial-defnorm defnorm0 :step i)
			  (and maxsteps (>= i maxsteps))))
	    finally
	    (awhen (and output (getf *linsolve-output* :final))
	      (funcall it :indentation indentation))
	    (return
	      (values
	       new-sol
	       ;; we return also a status report in the form of a property list
	       (list :solution new-sol :res res :defnorm defnorm
		     :last-step-reduction (and defnorm previous-defnorm
					       (safe-divide-by-zero defnorm previous-defnorm))
		     :steps i :reduction red-factor :convergence-rate
		     (if (numberp red-factor) (expt red-factor (/ i)) red-factor)
		     :status (if success-p :success :failure))))))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Object-oriented interface
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass <linear-solver> (<solver>)
  ((iteration :accessor iteration :initarg :iteration))
  (:documentation "Linear solver class."))
  
(defmethod solve ((solver <linear-solver>) &rest rest
		  &key matrix rhs solution residual &allow-other-keys)
  "A better interface for the old linsolve function."
  (with-slots (iteration maxsteps threshold threshabs reduction success-if failure-if
			 residual-norm output)
    solver
    (apply #'linsolve matrix rhs :sol solution :res residual :iteration iteration
	   :maxsteps maxsteps :threshold threshold :threshabs threshabs
	   :reduction reduction :success-if success-if :failure-if failure-if
	   :residual-norm residual-norm :output output
	   rest)))

(defparameter *lu-solver*
  (make-instance '<linear-solver> :iteration *lu-iteration* :maxsteps 1)
  "LU decomposition without pivoting.")

(defclass <special-solver> (<solver>)
  ((solver-function :reader solver-function :initarg :solver-function :type function))
  (:documentation "If you happen to have a problem-adapted solver given as
a function, you may use this base class."))

(defmethod solve ((solver <special-solver>) &rest parameters)
  (apply (solver-function solver) parameters))

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
       (solve (slot-value solve-it 'solver)
	      :matrix mat :rhs b :solution x :residual r))
   :residual-after nil))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Testing: (test-linsolve)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun test-linsolve ()
  "Linear solving on different types of vector/matrix representations."
  (let ((lu *lu-iteration*)
	(jac *undamped-jacobi*)
	(gs *gauss-seidel*))
    
    ;; application to matlisp matrices
    (let ((A [[2.0 -1.0 0.0]' [-1.0 2.0 -1.0]' [0.0 -1.0 2.0]'])
	  (b [1.0 2.0 1.0]'))
      (linsolve A b :output t :iteration lu)
      (linsolve A b :output t :iteration jac)
      (linsolve A b :output t :iteration gs)
      (linsolve A b :output t :iteration (make-instance '<sor> :omega 1.18))
      ;; new interface
      (solve *lu-solver* :matrix A :rhs b)
      (solve (make-instance '<linear-solver> :iteration jac :output t :maxsteps 10)
	     :matrix A :rhs b))
      
    ;; application to sparse matrices
    (let* ((constantly-1 (constantly 1))
	   (b (make-instance '<sparse-vector> :key->size constantly-1))
	   (A (make-sparse-matrix
	       :row-key->size constantly-1
	       :col-key->size constantly-1
	       :keys->pattern (constantly (full-crs-pattern 1 1))))
	   (sol (make-instance '<sparse-vector> :key->size constantly-1)))
      (setf (vec-ref (vec-ref b 1) 0) 1.0d0)
      (setf (vec-ref (vec-ref b 2) 0) 2.0d0)
      (setf (vec-ref (vec-ref b 3) 0) 1.0d0)
      
      (setf (vec-ref (mat-ref A 1 1) 0) 2.0d0)
      (setf (vec-ref (mat-ref A 2 2) 0) 2.0d0)
      (setf (vec-ref (mat-ref A 3 3) 0) 2.0d0)
      (setf (vec-ref (mat-ref A 1 2) 0) -1.0d0)
      (setf (vec-ref (mat-ref A 2 1) 0) -1.0d0)
      (setf (vec-ref (mat-ref A 2 3) 0) -1.0d0)
      (setf (vec-ref (mat-ref A 3 2) 0) -1.0d0)
      (assert (eql (matrix-ref (mat-ref A 3 2) 0 0) -1.0d0))
      (print-svec (linsolve A b :output t :iteration lu))
      (print-svec (linsolve A b :output t :iteration gs))
      (print-svec (solve *lu-solver* :matrix A :rhs b :solution sol)))
    )
  )

(tests::adjoin-femlisp-test 'test-linsolve)

