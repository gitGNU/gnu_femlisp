(in-package :ddo-femlisp-test)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Parallel part
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(format t "~&*** Parallel calculation:~%")

;;; Initialization/termination 

;;; (lfarm:end-kernel)
(load "connect-to-mpi-workers.lisp")

(ddo (defparameter *n* (mpi-comm-size))
     *n*)
(ddo (in-package :ddo-femlisp-test)
     *package*)

(ddo
  (princ-to-string (sb-cpu-affinity:get-cpu-affinity-mask)))

(ddo
  (progn (fl.parallel::new-kernel 6) nil))
(ddo
  (progn (fl.parallel::end-kernel) nil))

(time
 (ddo
   (defparameter *my-mesh* nil)
   (defparameter *ansatz-space* nil)
   (defparameter *rhs* nil)
   (defparameter *mat* nil)
   (format t "GC~%")
   (sb-ext:gc :full t)
   (synchronize)
   (reset-distributed-objects)
   (defparameter *n* (mpi-comm-size))
   ;;(dbg-on :ddo-refine)
   ;;(dbg-on :distribute)
   (let ((nr-workers *n*))
     (when (< (mpi-comm-rank) nr-workers)
       (let ((dim 3))
         (setq *my-mesh*
               (uniform-mesh-on-box-domain (n-cube-domain dim) (make-fixnum-vec dim 2))))
       (ddo-femlisp::distribute-mesh *my-mesh* nr-workers)))
   (change-class *my-mesh* '<hierarchical-mesh>)
   (format t "Refining~%")
   (loop repeat 3 do (refine *my-mesh*))
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
   ))

(ddo (distributed-data))
(ddo (nr-of-levels *my-mesh*))
(ddo (nr-of-cells *my-mesh* :highest))
(ddo (x<-0 *rhs*)
     (x<-0 *mat*))

(ddo (dbg-on :memoize))
(ddo (dbg-on :distribute))

(loop repeat 3 do
  (time
   (ddo
     (time
      (assemble-interior *ansatz-space* :surface :matrix *mat* :rhs *rhs*))
     )))

(ddo (asdf:req :ddo-femlisp))

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
(loop repeat 3 do (refine *my-mesh*))
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
(loop repeat 4 do
  (format t "Discretizing~%")
  (time
   (assemble-interior *ansatz-space* :surface :rhs *rhs* :matrix *mat*)))

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

