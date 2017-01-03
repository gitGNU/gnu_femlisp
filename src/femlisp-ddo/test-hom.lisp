(in-package :ddo-femlisp-test)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Parallel part
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(format t "~&*** Parallel calculation:~%")

;;; Initialization/termination 

;;; (lfarm:end-kernel)
(load "../connect-to-mpi-workers.lisp")

(setq ddo::*on-controller-p* nil)

(ddo (defparameter *n*
       (if ddo::*on-controller-p*
           1
           (mpi-comm-size)))
     *n*)

(defparameter *n*
  (if ddo::*on-controller-p*
      1
      (lfarm:kernel-worker-count)))

(ddo (in-package :ddo-femlisp-test)
     *package*)

;; (ddo (princ-to-string (sb-cpu-affinity:get-cpu-affinity-mask)))
;; (ddo (progn (fl.parallel::new-kernel 6) nil))


(ddo
  (defparameter *mesh* nil)
  (defparameter *ansatz-space* nil)
  ;; (defparameter *rhs* nil)
  ;; (defparameter *mat* nil)
  ;; (defparameter *res* nil)
  ;; (defparameter *sol* nil)
  ;;(defparameter *corr* nil)
  ;; (defparameter *diag* nil)
  (when (mpi-initialized)
    (format t "GC~%")
    (sb-ext:gc :full t)
    (synchronize)
    (reset-distributed-objects))
  (let* ((dim 2) (order 1)
         (domain (?1 (n-cell-with-ball-hole dim)
                     (n-cell-with-ball-inlay dim)))
         (problem (?1 (fl.application::cdr-cell-problem
                       domain :diffusion (fl.application::ball-diffusion dim .1)
                       :corner-constraint-p nil)
                      (fl.application::elasticity-inlay-cell-problem domain)))
         (mesh (make-mesh-from domain :parametric :from-domain))
         (*output-depth* 2)
         saved-effective-tensor)
    ;; distribute mesh
    (when (mpi-initialized)
      (assert (<= *n* (mpi-comm-size)))
      (distribute-mesh mesh *n*))
    (defparameter *mesh* (change-class mesh '<hierarchical-mesh>))
    (setq *ansatz-space*
          (let ((fe-class (lagrange-fe order :nr-comps (nr-of-components problem) :type :uniform)))
           (make-fe-ansatz-space fe-class problem *mesh*)))
    ;; (setq *rhs* (make-ansatz-space-vector *ansatz-space*))
    ;; (setq *mat* (make-ansatz-space-automorphism *ansatz-space*))
    ;; (assemble-interior *ansatz-space* :surface :rhs *rhs* :matrix *mat*)
    ;; (setq *res* (make-ansatz-space-vector *ansatz-space*))
    ;; (setq *sol* (make-domain-vector-for *mat*))
    ;; (setq *corr* (make-domain-vector-for *mat*))
    ;; (for-each-key (lambda (key) (vref *res* key) (vref *sol* key) (vref *corr* key)) *rhs*)
    ;; (setq *diag* (ddo-femlisp::diagonal-asa *mat*))
    ))

;;; (ddo- (print (mapcar #'midpoint (cells-of-highest-dim *mesh*))))

(ddo

  (let* ((problem (problem *ansatz-space*))
         (mesh (mesh *ansatz-space*))
         (disc (fe-class *ansatz-space*))
         (dim (dimension mesh))
         (ncomps (nr-of-components problem))
         (order (discretization-order *ansatz-space*))
         (levels 4))
    
    ;; solve
    (storing
      (solve
       (blackboard
        :ansatz-space *ansatz-space*
        :problem problem
        :fe-class disc
        :mesh mesh
        :estimator nil
        :indicator (make-instance '<uniform-refinement-indicator>)
        :success-if `(>= :max-level ,(1- levels))
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
                         :success-if `(or (zerop :defnorm) (and (> :step 2) (> :step-reduction 1.0) (< :defnorm 1.0e-8)))
                         :failure-if `(and (> :step 100) (> :step-reduction 1.0) (> :defnorm 1.0e-8))))
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
                                     (and tensor (> dim 1) (> ncomps 1) (mref (mref tensor 0 0) 1 1))
                                     (and tensor (> dim 1) (> ncomps 1) (mref (mref tensor 1 0) 1 0))))))))
        :output 1)
       ))
    ;; compute the homogenized coefficient
    (or saved-effective-tensor
        (fl.application::effective-tensor *result*))
    ))


(ddo (fl.mesh::substance-boundary-cells (mesh *ansatz-space*)))
(ddo- (fl.plot::plot (mesh *ansatz-space*)))
(ddo
  (let ((blackboard fl.application::*result*))
    (with-items (&key mesh strategy refinement-table refined-cells)
        blackboard
      (indicate (slot-value strategy 'fl.strategy::indicator) blackboard)
      (setf refined-cells
            (let ((ht refinement-table))
              (nth-value 1 (refine mesh :indicator #'(lambda (cell)
                                                       (nth-value 1 (gethash cell ht)))))))
      (assert (not (skel-empty-p refined-cells)))
      ;; update solution and interpolation and projection operators
      (fl.strategy::update-I-P-sol blackboard))))

(ddo
  (with-items (&key matrix) *result*
    (defparameter *corr* (make-domain-vector-for matrix))))
(ddo-
  (with-items (&key interior-rhs rhs)
      *result*
    (print rhs)
    (print interior-rhs)
    ;;(copy! rhs *corr*)
    ;;(make-consistent *corr*)
    ))

(ddo (dbg-on :constraints))
(ddo
  (fe-discretize *result*))
(ddo
  (with-items (&key interior-rhs rhs) *result*
    (copy! rhs *corr*)
    (make-consistent *corr*)))

(ddo-
  ;;(show (?1 *rhs* (getbb *result* :rhs))))
  (print (sparse-vector-to-list
          (?2 *corr* (getbb *result* :rhs)))))

(ddo- (fl.plot:plot (cells-on-level *mesh* 1)))

(ddo-
  (with-items (&key interpolation projection rhs solution)
      fl.application::*result*
    (show rhs)))
(fl.iteration::terminate-p)

(ddo-
  (with-items (&key mesh)
      *result*
    (print (surface-cells-of-highest-dim mesh))))
(fl.iteration::terminate-p)

(ddo- (check *mesh*))
(fl.plot:plot *mesh*)
(ddo- (display-ht (fl.mesh::skeleton-substance *mesh*)))

(ddo- (skeleton-boundary *mesh*))
(ddo (refinement-interface *mesh*))
(ddo (nr-of-cells *mesh*))
 
*res*
(ddo (dbg-on :partition-graph))

(describe *mesh*)
(setf *graph* (replace-nodes-by-numbers (mesh-graph *mesh*)))
(partition-graph 2)
(ddo- (asdf:req :ddo-femlisp))
(let ((*graph* (mesh-graph *mesh*)))
  (partition-graph 2))




(ddo- (dbg t "~A" (cells-of-highest-dim *mesh*)))
(ddo (distributed-data))
(elasticity-interior-effective-coeff-demo
     (elasticity-inlay-cell-problem (n-cell-with-ball-hole 3))
     :order 5 :levels 1 :plot nil :output 2)

(time
 (ddo
   (defparameter *my-mesh* nil)
   (defparameter *ansatz-space* nil)
   (defparameter *rhs* nil)
   (defparameter *mat* nil)
   (defparameter *res* nil)
   (defparameter *sol* nil)
   (format t "GC~%")
   (sb-ext:gc :full t)
   (synchronize)
   (reset-distributed-objects)
   (defparameter *n*
     (and (mpi-initialized) (mpi-comm-size)))
   ;;(dbg-on :ddo-refine)
   ;;(dbg-on :distribute)
   (let ((dim 2))
     (setq *my-mesh*
           (uniform-mesh-on-box-domain (n-cube-domain dim) (make-fixnum-vec dim 2)))
     (when (mpi-initialized)
       (assert (<= *n* (mpi-comm-size)))
       (ddo-femlisp::distribute-mesh *my-mesh* *n*)))
   (change-class *my-mesh* '<hierarchical-mesh>)
   (format t "Refining~%")
   (loop repeat 4 do (refine *my-mesh*))
   (format t "Ready with refinement~%")
   (setq *ansatz-space*
         (let* ((order 1)
                (problem
                  (?1 (fl.cdr::cdr-model-problem
                       (domain *my-mesh*)
                       :dirichlet nil :reaction 1.0
                       :source (lambda (x)
                                 (reduce #'* x :key (lambda (x) (* x (- 1.0 x))))))
                      (ellsys-model-problem
                       (domain *my-mesh*) '((u 1))
                       :a (diagonal-sparse-tensor (vector (eye (dimension *my-mesh*))))
                       :r (diagonal-sparse-tensor (diag #(10.0)))
                       :f (lambda (x)
                            (vector (ensure-matlisp 1.0))))))
                (fe-class (lagrange-fe order :nr-comps 1 :type :uniform)))
           (make-fe-ansatz-space fe-class problem *my-mesh*)))
   (setq *rhs* (make-ansatz-space-vector *ansatz-space*))
   (setq *mat* (make-ansatz-space-automorphism *ansatz-space*))
   (assemble-interior *ansatz-space* :surface :rhs *rhs* :matrix *mat*)
   (setq *res* (make-ansatz-space-vector *ansatz-space*))
   (setq *sol* (make-domain-vector-for *mat*))
   (defparameter *corr* (make-domain-vector-for *mat*))
   (defparameter *diag* (ddo-femlisp::diagonal-asa *mat*))
   ))

(ddo
  (let* ((smoother (make-instance 'ddo-femlisp::<distributed-jacobi>
                                  :damp 0.5))
         (cs (fl.geomg::geometric-cs
              :gamma 1 :smoother smoother :pre-steps 1 :post-steps 0
              :base-level 0
              :coarse-grid-iteration
              (make-instance 'fl.iteration::<multi-iteration>
                             :base smoother :nr-steps 1)
              :combination :additive))
         (bpx (make-instance 'ddo-femlisp::<distributed-cg>
                             :preconditioner cs :restart-cycle 30)))
    (fl.iteration::linsolve
     *mat* *rhs* :sol *sol* :res *res*
     :iteration bpx :output t
     :maxsteps 20)))
