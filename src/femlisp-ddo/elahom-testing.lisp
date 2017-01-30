(in-package :femlisp-ddo-test)

;;; testing routines for the elahom problem

(defparameter *n*
  (and (mpi-initialized)
       (mpi-comm-size)))

(defparameter *mesh* nil)
(defparameter *ansatz-space* nil)
(defparameter *rhs* nil)
(defparameter *mat* nil)
(defparameter *res* nil)
(defparameter *sol* nil)
(defparameter *corr* nil)
(defparameter *diag* nil)

(defun elahom-calculation (levels &optional (base-level 0))
  (setq *mesh* nil
        *ansatz-space* nil
        *rhs* nil
        *mat* nil
        *res* nil
        *sol* nil
        *corr* nil
        *diag* nil)
  (when (mpi-initialized)
    (format t "GC~%")
    (sb-ext:gc :full t)
    (synchronize)
    (reset-distributed-objects))
  
  ;;
  (let* ((dim 3) (order 5)
         (domain (?1 (n-cell-with-ball-hole dim)
                     (n-cell-with-ball-inlay dim)))
         (problem (fl.application::elasticity-inlay-cell-problem domain))
         (mesh (make-mesh-from domain :parametric :from-domain))
         (*output-depth* 2))
    ;; start from refined mesh if desired
    (loop repeat base-level do (setf mesh (refine mesh)))
    ;; distribute mesh
    (when (mpi-initialized)
      (assert (<= *n* (mpi-comm-size)))
      (distribute-mesh mesh *n*))
    (setq *mesh* (change-class mesh '<hierarchical-mesh>))
    (setq *ansatz-space*
          (let ((fe-class (lagrange-fe order :nr-comps dim :type :uniform)))
            (make-fe-ansatz-space fe-class problem *mesh*)))
    (setq *rhs* (make-ansatz-space-vector *ansatz-space*))
    (setq *mat* (make-ansatz-space-automorphism *ansatz-space*))
    (assemble-interior *ansatz-space* :surface :rhs *rhs* :matrix *mat*)
    (setq *res* (make-ansatz-space-vector *ansatz-space*))
    (setq *sol* (make-domain-vector-for *mat*))
    (setq *corr* (make-domain-vector-for *mat*))
    (for-each-key (lambda (key) (vref *res* key) (vref *sol* key) (vref *corr* key)) *rhs*)
    (setq *diag* (femlisp-ddo::diagonal-asa *mat*)))
  
  ;;
  (ddo-performance-check)
  ;;
  (let* ((problem (problem *ansatz-space*))
         (mesh (mesh *ansatz-space*))
         (disc (fe-class *ansatz-space*))
         (dim (dimension mesh))
         saved-effective-tensor)
  
    ;; solve
    (fl.application::storing
     (solve
      (blackboard
       :ansatz-space *ansatz-space*
       :problem problem
       :fe-class disc
       :mesh mesh
       :estimator nil
       :indicator (make-instance '<uniform-refinement-indicator>)
       :success-if `(>= :max-level ,(- levels 1 base-level))
       :solver
       (let* ((smoother (make-instance '<distributed-jacobi> :damp 1.0))
              (cs (geometric-cs
                   :gamma 1 :smoother smoother :pre-steps 1 :post-steps 0
                   :coarse-grid-iteration
                   (make-instance '<multi-iteration> :base smoother :nr-steps 1)
                   :combination :additive))
              (bpx (make-instance '<cg> :preconditioner cs :restart-cycle 30)))
         (make-instance '<linear-solver>
                        :iteration bpx
                        :output t
                        :success-if `(or (zerop :defnorm) (and (> :step 2) (> :step-reduction 1.0) (< :defnorm 1.0e-7)))
                        :failure-if `(and (> :step 200) (> :step-reduction 1.0) (> :defnorm 1.0e-7))))
       :plot-mesh nil
       :observe
       (append *stationary-fe-strategy-observe*
               (list fl.strategy::*mentries-observe*)
               (list
                (list (format nil "~19@A~19@A~19@A" "A^00_00" "A^00_11" "A^10_10") "~57A"
                      #'(lambda (blackboard)
                          (let ((tensor (fl.application::effective-tensor blackboard)))
                            (setq saved-effective-tensor tensor)
                            (format nil "~19,10,2E~19,10,2E~19,10,2E"
                                    (and tensor (mref (mref tensor 0 0) 0 0))
                                    (and tensor (> dim 1) (mref (mref tensor 0 0) 1 1))
                                    (and tensor (> dim 1) (mref (mref tensor 1 0) 1 0))))))))
       :output 1)
      ))
    ;; compute the homogenized coefficient
    #+(or)
    (awhen (and effective-tensor-p
                (or saved-effective-tensor (effective-tensor *result*)))
           (format t "The effective elasticity tensor is:~%~A~%" it)
           it)
    )
  )
