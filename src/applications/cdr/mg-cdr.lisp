;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; mg-cdr.lisp - Solving CDR problems with multigrid
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

(in-package :application)

(time
 (let* ((dim 3) (level 2) (order 4)
	(problem (laplace-test-problem-on-domain (n-cube-domain dim)))
	(smoother #+(or)(make-instance '<local-bgs>)
		  #-(or)(make-instance '<gauss-seidel>))
	(v-cycle (make-instance
		  '<geometric-cs>
		  ;;'<geometric-fas>
		  :gamma 1 :pre-steps 1 :pre-smooth smoother :post-steps 0
		  :fmg nil :base-level (1- level)
		  :coarse-grid-iteration *lu-iteration*
		  #+(or)(make-instance '<multi-iteration> :base smoother :nr-steps 1)))
	(mesh (uniformly-refined-hierarchical-mesh (domain problem) level))
	(fedisc (lagrange-fe order)))
   (multiple-value-bind (A b)
       (discretize-globally problem mesh fedisc)
     #+(or)(display A)
     #-(or)(nth-value 1 (linsolve A b :output t :iteration v-cycle :maxsteps 20)))))
;;(fe-value sol (make-double-vec dim 0.5d0))
     #+(or)(plot sol)
     ))))

;;; For a periodic problem we do a post-smoothing step to set the solution
;;; correct for plotting.  This should be automatized in the future, maybe
;;; by suitable :after-methods on gemm!, x+=A*y in the case of ansatz-space
;;; vectors satisfying constraints.
(time
 (let* ((dim 2) (level 4)
	(problem (make-instance
		  '<cdr-problem> :domain (n-cell-domain dim) :patch->coefficients
		  #'(lambda (patch)
		      (declare (ignore patch))
		      (list 'DIFFUSION
			    (make-<coefficient>
			     :input (make-<coefficient-input> :global t)
			     :eval #'(lambda (ci) (eye dim)))
			    'SOURCE
			    (make-<coefficient>
			     :input (make-<coefficient-input> :global t)
			     :eval #'(lambda (ci)
				       (if (>= (aref (ci-global ci) 0) 0.5) 1.0d0 -1.0d0)))))))
       (v-cycle (make-instance
		 '<geometric-cs>
		 ;;'<geometric-fas>
		 :gamma 1 :pre-steps 1 :post-steps 1
		 :fmg t :base-level 1
		 :coarse-grid-iteration
		 (make-instance '<multi-iteration> :base *gauss-seidel* :nr-steps 1))))
  (multiple-value-bind (A b)
      (problem-discretization problem :level level :order 1)
    ;;(show A)))
    ;;(show (aref (multigrid::a-vec (multilevel-decomposition v-cycle A)) 1))))

    ;;;(show (aref (multigrid::i-vec (multilevel-decomposition v-cycle A)) 1))))


    (let ((sol (linsolve A b :output t :iteration v-cycle :maxsteps 2)))
      (plot sol)
      ;;(fe-value sol (make-double-vec dim 0.5d0))
      ))))

;;; geometric V-cycle for higher-order problems
(let* ((dim 1) (level 1) (order 2)
       (problem (laplace-test-problem-on-domain (n-cube-domain dim)))
       (v-cycle (geometric-cs :base-level 1)))
  (multiple-value-bind (A b)
      (problem-discretization problem :level level :order order)
    #+(or)
    (let ((mg-data (multilevel-decomposition v-cycle A)))
      (show (aref (get-al mg-data :a-vec) 0)))
    #-(or)
    (let ((mg-data (multilevel-decomposition v-cycle A)))
      (show (aref (get-al mg-data :i-vec) 0)))
    #+(or)(plot (mesh b))
    #+(or)(linsolve A b :output t :iteration v-cycle :maxsteps 10)
    ))

(defun mg-cdr-tests ()

;;; test if CR is small enough
(let* ((dim 1) (level 3) (order 1)
       (problem (laplace-test-problem-on-domain (n-cube-domain dim)))
       (geomg (geometric-cs :gamma 1 :pre-steps 1 :post-steps 1 :base-level 1)))
  (multiple-value-bind (A b)
      (problem-discretization problem :level level :order order)
    #+(or)(linsolve A b :output t :iteration geomg :maxsteps 10)
    #+(or)(let ((interior-mat (getf (discretization-info A) :interior-matrix)))
	    (show (geomg::extend-level-matrix (extract-level interior-mat 1)
					      interior-mat :layers 0)))
    #+(or)(let ((assembly-line (multilevel-decomposition geomg A)))
      (show (aref (get-al assembly-line :a-vec) 0)))
    #-(or)
    (let ((cr (getf (nth-value 1 (linsolve A b :output t :iteration geomg :maxsteps 10))
		    :convergence-rate)))
      (assert (< cr 0.08))
      cr)))

)
