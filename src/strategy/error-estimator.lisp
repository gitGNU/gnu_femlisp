;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; error-estimator.lisp - Error estimation
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

(in-package :strategy)

(defclass <error-estimator> ()
  ()
  (:documentation "An error estimator is used as first argument in the
generic functions local-estimate, global-estimate, and estimate which work
on an assembly line."))

(defgeneric local-estimate (error-estimator assembly-line)
  (:documentation "Puts a hash-table of local error contributions on the
assembly line."))

(defgeneric global-estimate (error-estimator assembly-line)
  (:documentation "Puts an estimate for the total error on the assembly
line.  Will usually need the result of the local estimate."))

(defgeneric estimate (error-estimator assembly-line)
  (:documentation "Yields both local and global estimate."))

(defmethod estimate ((errest <error-estimator>) assembly-line)
  "Executes local and global error estimation in sequence."
  (local-estimate errest assembly-line)
  (global-estimate errest assembly-line))

(defun eta->p2-vec (eta problem mesh)
  "Maps eta into a P2 ansatz function.  This is used for plotting the error
distribution.  Note that for tetrahedra, no cell contributions are shown."
  (let* ((as (make-fe-ansatz-space (lagrange-fe 2) problem mesh))
	 (result (make-ansatz-space-vector as)))
    (dohash ((cell value) eta)
      (let ((key (cell-key cell mesh)))
	(when (plusp (funcall (key->size as) key))
	  (setf (vec-ref result key) [value]))))
    result))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; standard-error-estimator
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass <standard-error-estimator> (<error-estimator>)
  ()
  (:documentation "For this form of error estimator the local estimate is
further decomposed in the steps project-solution, compute-weight-function,
and compute-local-estimate, all of them acting on the assembly line."))

(defmethod compute-error-approximant ((errest <standard-error-estimator>) assembly-line)
  assembly-line)
(defmethod compute-weight-function ((errest <standard-error-estimator>) assembly-line)
  assembly-line)
(defmethod compute-local-estimate ((errest <standard-error-estimator>) assembly-line)
  assembly-line)

(defmethod local-estimate ((errest <standard-error-estimator>) assembly-line)
  (compute-error-approximant errest assembly-line)
  (compute-weight-function errest assembly-line)
  (compute-local-estimate errest assembly-line))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; compute-error-approximant
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass <difference-with-projection> ()
  ()
  (:documentation "This is a mixin to <standard-error-estimator> which
adapts project-solution to compute the difference between the solution and
its projection to a coarser level in the quantity :solution-increment on
the assembly line."))

(defmethod compute-error-approximant ((errest <difference-with-projection>) assembly-line)
  "Computes solution-increment computed as difference solution on the
finest mesh and a projection to the next-coarser level."
  (with-items (&key ansatz-space solution interpolation projection
				 P*sol I*P*sol solution-increment)
      assembly-line
    (when (and interpolation projection)
      (setq P*sol (m* projection solution))
      (setq I*P*sol (m* interpolation P*sol))
      (setq solution-increment
	    (sparse->ansatz-space-vector (m- solution I*P*sol)
					 :ansatz-space ansatz-space)))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; compute-local-estimate
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass <global-norm> ()
  ((global-p :initform 2 :initarg :global-p))
  (:documentation "This is a mixin to <standard-error-estimator> which
leads to computation of a p-norm of the local contributions.  If global-p
is :identity a simple summation is performed."))

(defclass <global-and-local-norm> (<global-norm>)
  ((local-p :initform 2 :initarg :local-p))
  (:documentation "This is a mixin to <standard-error-estimator> which
computes a local cell norm."))

(defmethod compute-local-estimate ((errest <global-and-local-norm>) assembly-line)
  (with-items (&key mesh solution-increment eta)
    assembly-line
    (when solution-increment
      (setq eta
	    (let* ((local-p (slot-value errest 'local-p))
		   (global-p (slot-value errest 'global-p))
		   (eta (make-hash-table)))
	      (doskel (cell mesh :dimension :highest :where :surface)
		(let ((local (cell-integrate cell solution-increment
					     :key #'(lambda (x) (expt (norm x local-p) global-p)))))
		  (setf (gethash cell eta) (expt local (/ global-p)))))
	      eta)))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; global-estimate
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmethod global-estimate ((errest <global-norm>) assembly-line)
  "Computes the global error from the local error contributions."
  (with-items (&key eta global-eta)
      assembly-line
    (when eta
      (setq global-eta
	    (let ((global-p (slot-value errest 'global-p)))
	      (case global-p
		(:identity
		 (loop for val being each hash-value of eta
		       summing val))
		(t (expt (loop for val being each hash-value of eta
			       summing (expt val global-p))
			 (/ global-p)))))))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Concrete error estimators
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Projection error estimator
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass <projection-error-estimator>
    (<difference-with-projection> <global-and-local-norm> <standard-error-estimator>)
  ()
  (:documentation "Estimates the error by measuring the difference between
the solution and a projected solution in a hierarchical mesh by a certain
norm given by local-p and global-p."))

;;; compute-error-approximant: <difference-with-projection>
;;; compute-residual: none
;;; compute-weight-function: none
;;; compute-local-estimate: <global-and-local-norm>

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Duality error estimator
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; This type of error estimator computes an error weight function by
;;; solving a dual problem with the error functional as right hand side in
;;; a larger ansatz-space.  This is multiplied with the residual computed
;;; in the larger ansatz-space.

(defclass <setup-enriched-ansatz-space> ()
  ()
  (:documentation "A mixin leading to the assembly of interpolated
solution, matrix and rhs in an enriched ansatz space."))

(defmethod compute-error-approximant ((errest <setup-enriched-ansatz-space>) assembly-line)
  "Setup of interpolated solution, matrix, and rhs in the enriched ansatz
space."
  (with-items (&key enlarged-as-assembly-line problem mesh refined-cells)
    assembly-line
    (unless enlarged-as-assembly-line
      (setf enlarged-as-assembly-line (make-assembly-line)))
    (let* ((as-low (get-al assembly-line :ansatz-space))
	   (order (discretization-order (fe-class as-low)))
	   (as-high (or (get-al enlarged-as-assembly-line :ansatz-space)
			(make-fe-ansatz-space (lagrange-fe (1+ order)) problem mesh))))
      (with-items (&key low->high high->low low->high-T
				     interior-matrix interior-rhs matrix)
	enlarged-as-assembly-line
	;; to be improved later avoiding reassembly on old regions
	(setf (get-al enlarged-as-assembly-line :ansatz-space) as-high)
	(assemble-constraints as-high)
	(setf low->high (transfer-matrix as-low as-high :no-slaves nil)
	      high->low (transfer-matrix as-high as-low :no-slaves nil)
	      low->high-T (transpose low->high))
	(setf (get-al enlarged-as-assembly-line :solution)
	      (m* low->high (get-al assembly-line :solution)))
	(setf interior-matrix nil
	      interior-rhs nil)
	(fe-discretize enlarged-as-assembly-line)))))

(defclass <solve-dual-problem> ()
  ((functional :accessor functional :initarg :functional)
   (solver :reader solver :initform nil :initarg :solver))
  (:documentation "A mixin leading to the computation of the solution of
the dual problem for the given functional.  The result is put in
:error-weight-function."))

(defmethod compute-weight-function ((errest <solve-dual-problem>) assembly-line)
  "Solves a dual problem in an ansatz-space of higher order.  This error
estimator is computationally more intensive than the original problem."
  (with-items (&key mesh problem strategy refined-cells
		    dual-problem-assembly-line enlarged-as-assembly-line)
      assembly-line
    (unless dual-problem-assembly-line
      (setf dual-problem-assembly-line (make-assembly-line)))
    (unless (get-al dual-problem-assembly-line :ansatz-space)
      (let* ((dual-problem #-(or) (dual-problem problem (functional errest))
			   #+(or) problem)
	     (as-low (get-al assembly-line :ansatz-space))
	     (order (discretization-order (fe-class as-low)))
	     (as-high (make-fe-ansatz-space (lagrange-fe (1+ order))
					    dual-problem mesh)))
	(setf (get-al dual-problem-assembly-line :ansatz-space) as-high)
	(setf (get-al dual-problem-assembly-line :solution)
	      (make-ansatz-space-vector as-high))))
    ;;
    (when refined-cells
      (unless (eq refined-cells (get-al dual-problem-assembly-line :refined-cells))
	;; there has been a refinement which we did not yet track, such
	;; that we have to update the dual solution
	(setf (get-al dual-problem-assembly-line :refined-cells) refined-cells)
	(update-I-P-sol dual-problem-assembly-line)))
    (with-items (&key interior-matrix matrix interior-rhs rhs solution)
      dual-problem-assembly-line
      ;; should be improved later to avoid reassembly
      (setf interior-matrix nil interior-rhs nil)
      (cond ((self-adjoint-p problem)
	     (setf matrix nil)
	     (cond ((eq (functional errest) :load-functional)
		    ;; energy error estimate
		    (setf rhs (get-al enlarged-as-assembly-line :rhs)))
		   (t (fe-discretize dual-problem-assembly-line)))
	     (setf matrix (get-al enlarged-as-assembly-line :matrix)))
	    (t (fe-discretize dual-problem-assembly-line)))
      (setf solution  (apply #'solve (or (solver errest) (solver strategy))
			     dual-problem-assembly-line)))
    ))

;;; compute-local-estimate

(defclass <local-test-with-dual> ()
  ()
  (:documentation "A mixin testing the residual with an approximation of a
dual solution.  To localize we subtract a local projection on lower order
polynomials."))

(defmethod compute-local-estimate ((errest <local-test-with-dual>) assembly-line)
  "Evaluates the duality error estimator given the dual solution and the
approximate solution u.  Works only for the :load-functional case at the
moment, because the right-hand side is not assembled for the primal
problem."
  (with-items (&key mesh enlarged-as-assembly-line dual-problem-assembly-line)
      assembly-line
    (when (and enlarged-as-assembly-line dual-problem-assembly-line)
      (let* (#+debug(solution (get-al assembly-line :solution))
	     #+debug(mat-low (get-al assembly-line :matrix))
	     (rhs-low (get-al assembly-line :rhs))
	     (rhs-high (get-al enlarged-as-assembly-line :rhs))
	     (mat-high (get-al enlarged-as-assembly-line :matrix))
	     (dual-sol (get-al dual-problem-assembly-line :solution))
	     (low->high (get-al enlarged-as-assembly-line :low->high))
	     (high->low (get-al enlarged-as-assembly-line :high->low))
	     (low->high-T (get-al enlarged-as-assembly-line :low->high-T))
	     (Isol (get-al enlarged-as-assembly-line :solution))
	     (mat-high*Isol (m* mat-high Isol))
	     (res-high (m- rhs-high mat-high*Isol))
	     (dual-low (m* high->low dual-sol))
	     (dual-proj (m* low->high dual-low))
	     (weight (m- dual-sol dual-proj))
	     #+debug(gal-mat (sparse-m* low->high-T (sparse-m* mat-high low->high)))
	     (part1 (dot dual-low (m- (m* low->high-T rhs-high) rhs-low)))
	     #+debug(part1b (dot-abs dual-low (m- (m* low->high-T rhs-high) rhs-low)))
	     #+debug(part2 (dot (m* (m- mat-low gal-mat) solution) dual-low))
	     #+debug(part2b (dot-abs (m* (m- mat-low gal-mat) solution) dual-low))
	     (part3 (dot res-high dual-sol))
	     #+debug(part3b (dot-abs res-high dual-sol)))
	#+debug
	(format t "rhs-c=~12,5,2E:~12,5,2E, mat-c=~12,5,2E:~12,5,2E, err=~12,5,2E:~12,5,2E~%"
		part1 part1b part2 part2b part3 part3b)
	#+debug
	(format t "low=~A, high=~A, delta=~A, est=~A~%" (dot rhs-high (m* low->high solution))
		(dot rhs-high dual-sol) (dot rhs-high (m- dual-sol (m* low->high solution)))
		(dot dual-sol res-high))
	;; note: better would be (+ part1 part2 part3), but it is
	;; computationally quite expensive
	(setf (get-al assembly-line :global-estimate-guess) (+ part1 part3))
	;; set key-error
	(let ((key->error (make-hash-table))
	      (eta (make-hash-table)))
	  (let ((sum 0.0d0)
		(sum-abs 0.0d0))
	    (for-each-key
	     #'(lambda (key)
		 (let ((dot (dot (vec-ref res-high key) (vec-ref weight key)))
		       (dot-abs (dot-abs (vec-ref res-high key) (vec-ref weight key))))
		   (incf sum dot)
		   (incf sum-abs dot-abs)
		   (setf (gethash key key->error) dot-abs)))
	     res-high)
	    (setf (get-al assembly-line :global-estimate-1) sum)
	    (setf (get-al assembly-line :global-estimate-2) sum-abs))
	  ;; distribute to cells of highest dimension in some way (ignoring
	  ;; contributions of lower dimensional cells).  maybe better would
	  ;; be if the indicator would directly compute something out of this
	  ;; data.
	  (doskel (cell mesh :where :surface :dimension :highest)
	    (setf (gethash cell eta) 0.0d0)
	    (incf (gethash cell eta)
		  (abs (gethash (cell-key cell mesh) key->error)))
	    (loop for side across (boundary cell) do
		  (incf (gethash cell eta)
			(* 0.5
			   (abs (gethash (cell-key side mesh) key->error))))))
	  (setf (get-al assembly-line :error-contributions) key->error)
	  (setf (get-al assembly-line :eta) eta)
	  #+(or)
	  (format t "guess=~12,5,2E  dot=~12,5,2E  dot-abs=~12,5,2E~%"
		  (get-al assembly-line :global-estimate-guess)
		  (get-al assembly-line :global-estimate-1)
		  (get-al assembly-line :global-estimate-2))
	  )))
  ;; pass on assembly-line
  assembly-line))

(defmethod global-estimate ((errest <local-test-with-dual>) assembly-line)
  (with-items (&key global-eta global-estimate-guess)
      assembly-line
    (setf global-eta (abs global-estimate-guess))))

;;; The concrete class
(defclass <duality-error-estimator>
    (<setup-enriched-ansatz-space> <solve-dual-problem>
     <local-test-with-dual> <standard-error-estimator>)
  ()
  (:documentation "Estimates the error by testing the difference z-IPz
against the residual.  Here z is the solution of a dual problem in an
enriched finite element space."))

;;; compute-error-approximant: none
;;; compute-residual: none
;;; compute-weight-function: <solve-dual-problem>
;;; compute-local-estimate: <local-test-with-dual>

;;;; Testing

(defun test-error-estimator ()
  (make-instance '<projection-error-estimator>)
  )

(tests::adjoin-femlisp-test 'test-error-estimator)
