(in-package :ddo-femlisp-test)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Parallel part
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(format t "~&*** Parallel calculation:~%")

;;; Initialization/termination 

;;; (lfarm:end-kernel)
;;; (load "../connect-to-mpi-workers.lisp")

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

;; (ddo- (load "ddo-femlisp/mpi-affinity.lisp"))
;; (ddo (princ-to-string (sb-cpu-affinity:get-cpu-affinity-mask)))
;; (ddo (progn (fl.parallel::end-kernel) (fl.parallel::new-kernel 6) nil))
;; (ddo fl.parallel::*kernel*)

(ddo (dbg-off :communication)
     (dbg-off :distribute)
     (dbg-on :log-times)
     (dbg-off :synchronization-time-in-approximation-step)
     (setq ddo::*debug-show-data* nil)
     (setq ddo::*communicate-with-all* t)
     #+(or)(progn (fl.parallel::end-kernel) (fl.parallel::new-kernel 1 nil) nil)
     (setq ddo::*report-ranks* '(0)))

(ddo- (Asdf:req :ddo-femlisp))

(ddo (elahom-calculation 1))
(ddo (elahom-calculation 2))

(ddo (elahom-calculation 2 1)
     (elahom-calculation 3 1)
     (elahom-calculation 4 1))

(ddo (elahom-calculation 3))
(ddo (progn (fl.parallel::end-kernel) (fl.parallel::new-kernel 6) nil)
     (elahom-calculation 2 1)
     (elahom-calculation 3 1)
     (elahom-calculation 4 1))

(ddo (elahom-calculation 3))

#|

(ddo- (Asdf:req :ddo-femlisp))

(ddo (with-items (&key mesh) *result*
       (mapcar #'midpoint (cells-of-highest-dim (cells-on-level mesh 0)))))

(defmethod mextreme ((mat standard-matrix))
  (let (min max)
    (for-each-entry (lambda (entry)
                      (if min
                          (setf min (min min entry)
                                max (max max entry))
                          (setf min entry
                                max entry)))
                    mat)
    (cons min max)))

(ddo (with-items (&key solution) *result*
       (let ((extremes (fe-extreme-values solution)))
         (print (m+ (car extremes) (cdr extremes)))))
(setq fl.matlisp::*print-matrix-pretty* t)
(setq *print-matrix* 10)

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
    (dbg-on :update-i-p-sol)
    (fl.strategy::update-I-P-sol blackboard)))
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
|#
