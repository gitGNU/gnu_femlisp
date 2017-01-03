(in-package :ddo-femlisp-test)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Parallel part
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(format t "~&*** Parallel calculation:~%")

;;; Initialization/termination 

;;; (lfarm:end-kernel)
(load "../connect-to-mpi-workers.lisp")

(ddo (defparameter *n* (mpi-comm-size))
     *n*)
(defparameter *n* (lfarm:kernel-worker-count))

(ddo (in-package :ddo-femlisp-test)
     *package*)

;; (ddo (princ-to-string (sb-cpu-affinity:get-cpu-affinity-mask)))
;; (ddo (progn (fl.parallel::new-kernel 6) nil))

(setq ddo::*on-controller-p* nil)

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
   (loop repeat 2 do (refine *my-mesh*))
   (format t "Ready with refinement~%")
   (setq *ansatz-space*
         (let* ((order 5)
                (problem
                  (?1 (fl.cdr::cdr-model-problem
                       (domain *my-mesh*)
                       :dirichlet nil
                       :diffusion (eye (dimension *my-mesh*))
                       :reaction 1.0
                       :source (lambda (x)
                                 (reduce #'* x :key (lambda (x) (* x (- 1.0 x))))))
                      (ellsys-model-problem
                       (domain *my-mesh*) '((u 1))
                       :a (diagonal-sparse-tensor (vector (eye (dimension *my-mesh*))))
                       :r (diagonal-sparse-tensor (vector #m(10.0)))
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

;;; BPX

(ddo
  (let* ((smoother (make-instance 'ddo-femlisp::<distributed-jacobi>))
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
     :maxsteps 10)))

;;; geometric multigrid

(ddo
  (let* ((smoother (?2 (make-instance 'fl.iteration::<gauss-seidel>)
                       (make-instance (?1 'ddo-femlisp::<distributed-jacobi>
                                          'fl.iteration::<jacobi>)
                                      :damp 0.5)))
         (cs (fl.geomg::geometric-cs
              :gamma 1 :smoother smoother :pre-steps 2 :post-steps 0
              :base-level 0
              :coarse-grid-iteration
              (make-instance 'fl.iteration::<multi-iteration>
                             :base smoother :nr-steps 30))))
    (fl.iteration::linsolve
     *mat* *rhs* :sol *sol* :res *res*
     :iteration cs :output t
     :maxsteps 20)))

(show *mat*)
(ddo- (dbg t "~A" (get-property *sol* :distributed-p)))
(ddo (Asdf:req :ddo-femlisp))
(ddo (unintern (intern "PARALLEL-NORM"  (find-package :ddo-femlisp-test)) (find-package :ddo-femlisp-test)))

;;; preconditioned CG iteration

(ddo (let* ((jac (make-instance 'ddo-femlisp::<distributed-jacobi>))
            (cg (make-instance 'ddo-femlisp::<distributed-cg>
                              :preconditioner jac)))
       (fl.iteration::linsolve *mat* *rhs* :sol *sol* :res *res*
                                           :iteration cg :output t
                                           :maxsteps 10)))
(ddo- (show *sol*))

;;; Jacobi iteration by hand

(loop for k from 0 below 10 do
  (ddo (fl.iteration::compute-residual *mat* *sol* *rhs* *res*))
  (ddo (format t "Residual: ~A ~A~%"
               (ddo-femlisp::l2-norm *res*)
               (ddo-femlisp::parallel-norm *res* *corr*)))
  (ddo (ddo-femlisp::jacobi-correction *sol* 0.5 *corr* *diag* *res*)))

;;; Jacobi iteration with the Femlisp solver classes

(ddo
 (let ((jac (make-instance 'ddo-femlisp::<distributed-jacobi> :damp 0.5)))
   (fl.iteration::linsolve *mat* *rhs*
                           :sol *sol* :res *res*
                           :iteration jac :output t
                           :maxsteps 1000)))


;;; CG iteration

(ddo (let ((cg (make-instance 'ddo-femlisp::<distributed-cg>
                              :damp 0.5)))
       (fl.iteration::linsolve *mat* *rhs* :sol *sol* :res *res*
                                           :iteration cg :output t
                                           :maxsteps 10)))

;;; Multigrid iteration

(ddo (let* ((smoother (make-instance 'ddo-femlisp::<distributed-jacobi> :damp 0.5))
            (cs (fl.geomg::geometric-cs
                 :gamma 1 :smoother smoother :pre-steps 4 :post-steps 4
                 :base-level 0
                 :coarse-grid-iteration
                 (make-instance 'fl.iteration::<multi-iteration>
                                :base smoother :nr-steps 10))))
       (fl.iteration::linsolve
        *mat* *rhs* :sol *sol* :res *res*
                    :iteration cs :output t
                    :maxsteps 20)))


       (defparameter *iterator* (fl.iteration::make-iterator cs *mat*))))

(ddo (defparameter *mat2*
       (lret ((l-mat (fl.discretization::make-ansatz-space-automorphism (ansatz-space *mat*))))
         (fl.discretization::assemble-interior *ansatz-space* :refined :matrix l-mat :level 0))))

(ddo- (let ((data (slot-value *iterator* 'fl.multigrid::mg-data)))
        (loop for a across (getbb data :a-vec)
              do
                 (dbg t "~A ~D ~A~%" a (fl.matlisp::nrows a) (get-property a :distributed-p)))))


(format t "Residual: ~8,2g~%" 1.2341243e-15)

(ddo (asdf:req :ddo-femlisp))
(ddo (load "ddo-discretize"))
(ddo (dbg-on :merger))
(ddo (get-property *sol* :distributed-p))
(ddo- (show *diag*))
(ddo- (show *mat*))
(ddo- (show *rhs*))
(ddo- (show *res*))
(ddo- (show *sol*))
  (format t (ddo- (show *sol*))
(ddo (aso-make-distributed *diag*))
(ddo
  (copy! *rhs* *res*)
  (gemm! 1.0 *mat* *sol* 1.0 *res*)
  )
(ddo- (show *res*))

(ddo (distributed-data))
(ddo (nr-of-levels *my-mesh*))
(ddo (nr-of-cells *my-mesh* :highest))
(ddo (x<-0 *rhs*)
     (x<-0 *mat*))

(loop repeat 3 do
  (time
   (ddo
     (time
      (assemble-interior *ansatz-space* :surface :rhs *rhs* :matrix *mat*))
     )))

;;; Serial calculation

(asdf:req :ddo-femlisp)
(in-package :ddo-femlisp-test)
(format t "~&*** Serial calculation:~%")
(defparameter *my-mesh* nil)
(defparameter *ansatz-space* nil)
(defparameter *rhs* nil)
(defparameter *mat* nil)
(format t "GC~%")
(sb-ext:gc :full t)
(defparameter *n* 1)
(let ((dim 3))
  (setq *my-mesh*
        (uniform-mesh-on-box-domain (n-cube-domain dim) (make-fixnum-vec dim 2))))
(change-class *my-mesh* '<hierarchical-mesh>)
(format t "Refining~%")
(loop repeat 0 do (refine *my-mesh*))
(format t "Ready with refinement~%")
(setq *ansatz-space*
      (let* ((order 5)
             (problem (ellsys-model-problem
                       (domain *my-mesh*) '((u 1))
                       :r (diagonal-sparse-tensor (vector #m(1.0)))
                       :f (lambda (x)
                            (vector (ensure-matlisp 1.0)))))
             (fe-class (lagrange-fe order :nr-comps 1 :type :uniform)))
        (make-fe-ansatz-space fe-class problem *my-mesh*)))
(setq *rhs* (make-ansatz-space-vector *ansatz-space*))
(setq *mat* (make-ansatz-space-automorphism *ansatz-space*))
(loop repeat 1 do
  (sb-ext:gc :full t)
  (sb-sprof:reset)
  (sb-sprof:start-profiling)
  (format t "Discretizing~%")
  (time
   (assemble-interior *ansatz-space* :surface :rhs *rhs* :matrix *mat*))
  (sb-sprof:stop-profiling)
  (sb-sprof:report :type :flat :sort-by :cumulative-samples :max 200)
      )

(ddo *features*)
(ddo (load "packages.lisp"))
(ddo (load "ddo.lisp"))
(ddo (load "ddo-discretize.lisp"))
(ddo (load "ddo-mesh.lisp"))
(ddo (asdf:req :ddo-femlisp))

(ddo (for-each-entry-and-key 
      (lambda (entry rk ck)
        (format t "~A ~A ~A ~A ~A ~A~%"
                entry rk ck (distributed-p entry)
                (distributed-p (representative rk)) (distributed-p (representative ck)))
        (force-output))
      *mat*))

(ddo (distributed-p *mat*))
(ddo (asdf:req :ddo-femlisp))


(ddo (setq fl.debug::*dbg-ids* nil))
(ddo (dbg-on :merger))
(ddo (dbg-on :ddo-discretize))
(ddo fl.debug::*dbg-ids*)
(ddo (get-property *rhs* :distributed-p))

#+(or)
(ddo (loop repeat 3 do (refine *my-mesh*))
     (nr-of-surface-cells *my-mesh*))

(loop repeat 2 do
  (time
   (ddo
     (let* ((order 5)
            (problem (ellsys-model-problem
                      (domain *my-mesh*) '((u 1))
                      :r (diagonal-sparse-tensor (vector #m(1.0)))
                      :f (lambda (x)
                           (vector (ensure-matlisp (float #I"sin(2*pi*norm(x))" 1.0))))))
            (fe-class (lagrange-fe order :nr-comps 1 :type :uniform))
            (ansatz-space (make-fe-ansatz-space fe-class problem *my-mesh*)))
       (setq *rhs* (make-ansatz-space-vector ansatz-space))
       (assemble-interior ansatz-space :surface :rhs *rhs*)
       (nr-of-entries *rhs*)))))

#|

;;; For interactive work

;;; reset
(ddo (setq *my-mesh* nil
           *rhs* nil)
     (sb-ext:gc :full t)
     (synchronize)
     (reset-distributed-objects))

(ddo (mpi-comm-rank))
(ddo (dbg-on :distribute))
(ddo (distributed-data))
(ddo (synchronize))
(ddo (reset-distributed-objects))
(ddo (sb-ext:gc :full t))
(ddo (load "ddo-discretize.lisp"))
(ddo (load "ddo.lisp"))

|#

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Finishing
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(lfarm:end-kernel)
(sb-ext:quit)

