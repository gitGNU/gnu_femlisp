(in-package :fl.application)

;;; testing routines for the elahom problem

(defun initialize-elahom-calculation (dim order levels)
  (let* ((domain (n-cube-domain dim))
         (problem (elasticity-model-problem
                   domain :lambda 1.0 :mu 1.0
                          :force (ellsys-one-force-coefficient dim 1)
                          :dirichlet nil)))
    (storing
      (solve
       (blackboard
	:problem problem
	:fe-class (lagrange-fe order :nr-comps dim :type :gauss-lobatto)
	:estimator (make-instance '<projection-error-estimator>)
	:indicator (make-instance '<largest-eta-indicator> :fraction 1.0)
	:success-if `(>= :max-level ,(1- levels))
	:solver
        (let* ((smoother (make-instance '<jacobi> :damp 1.0))
                (cs (geometric-cs
                     :gamma 1 :smoother smoother :pre-steps 1 :post-steps 0
                     :coarse-grid-iteration
                     (make-instance '<multi-iteration> :base smoother :nr-steps 1)
                     :combination :additive))
                (bpx (make-instance '<cg> :preconditioner cs :restart-cycle 30)))
           (make-instance '<linear-solver>
                          :iteration bpx
                          :success-if `(> :step 2)
                          :failure-if `(> :step 100)))
	:plot-mesh nil
	:observe *stationary-fe-strategy-observe*
        :output 2)))))

(initialize-elahom-calculation 3 5 2)
