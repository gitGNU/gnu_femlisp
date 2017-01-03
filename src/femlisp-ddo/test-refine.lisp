(in-package :ddo-femlisp-test)

;;; (lfarm:end-kernel)
(load "connect-to-mpi-workers.lisp")
#+(or)
(setf lfarm:*kernel*
      (lfarm:make-kernel
       '(
         ("norton" 20000)
         ("norton" 20001)
         ("norton" 20002)
         ("norton" 20003))))

;;; number of workers
(defparameter *n* (lfarm:kernel-worker-count))
(ddo (defparameter *n* (mpi-comm-size)) *n*)

(ddo (load "packages.lisp"))
;;; (ddo (asdf:load-system :ddo-femlisp))  ; Problem: asdf is not working in parallel
(ddo (in-package :ddo-femlisp-test))

(ddo (mpi-comm-rank))
(ddo (dbg-on :distribute))
(ddo (distributed-data))
(ddo (synchronize))
(ddo (reset-distributed-objects))
(ddo (sb-ext:gc :full t))
(ddo (load "ddo-refine.lisp"))
(ddo (load "ddo.lisp"))

(fl.plot:plot (refine (skeleton (n-cube 3))))
(ddo (midpoint (first (cells-of-dim *my-mesh* 1))))
(cells-of-dim (refine (skeleton (n-cube 2))) 0)

;;; local
(defparameter *my-mesh*
  (uniform-mesh-on-box-domain (n-cube-domain 2) #(2 2)))
(defparameter *my-refined-meshes*
  (make-array 0 :fill-pointer t))
(vector-push-extend *my-mesh* *my-refined-meshes*)

(time
 (loop repeat 8 do
   (let ((mesh (vector-last *my-refined-meshes*)))
     (vector-push-extend (refine mesh) *my-refined-meshes*))))

(doskel (cell *my-mesh*)
  (make-distributed-object cell '(0)))
(synchronize)

(ddo (dbg-off :distribute))

;;; remote

(ddo (defparameter *my-mesh* nil)
     (defparameter *my-refined-meshes* nil))

(ddo
  (let ((nr-workers *n*))
    (when (< (mpi-comm-rank) nr-workers)
      (setq *my-mesh*
            (uniform-mesh-on-box-domain (n-cube-domain 2) #(2 2)))
      (ddo-femlisp::distribute-mesh *my-mesh* nr-workers))))

(ddo (distributed-data))

(time
 (ddo
   (setq *my-refined-meshes*
         (make-array 0 :fill-pointer t))
   (vector-push-extend *my-mesh* *my-refined-meshes*)
   (time
    (loop repeat 8 do
      (let ((mesh (vector-last *my-refined-meshes*)))
        (vector-push-extend (refine mesh) *my-refined-meshes*))))))

(ddo (nr-of-cells (vector-last *my-refined-meshes*)))
(vector-last nil)
(ddo (room))

(ddo
  (when (< (mpi-comm-rank) 2)
    (nr-of-cells (vector-last *my-refined-meshes*))))

(ddo (defparameter *my-doubly-refined-mesh*
       (refine *my-refined-mesh*)))

;;; reset
(ddo (setq *my-mesh* nil
           *my-refined-meshes* nil)
     (sb-ext:gc :full t)
     (synchronize))
