(in-package :fl.dictionary)

(defun extract-complete-subgraph (objects get-dependents &optional result)
  "Extract a complete subgraph which contains the vertices in objects.
The result is a hash-table mapping vertices to lists of depending vertices.
If the hash-table in @arg{result} is given it is assumed to contain a
partial result of this operation."
  (lret ((result (or result (make-hash-table))))
    (dolist (object objects)
      (let ((dependents (coerce (funcall get-dependents object) 'list)))
        (acond ((gethash object result)
                (assert (set-equal dependents it)))
               (t (setf (gethash object result)
                        dependents)
                  (extract-complete-subgraph dependents get-dependents result)))))))

(defun invert-graph (graph)
  "Invert the directed graph given in @arg{graph} as a hash-table."
  (lret ((result (make-hash-table)))
    (maphash
     (lambda (object dependents)
       (dolist (dependent dependents)
         (pushnew object (gethash dependent result))))
     graph)))

(defun get-dependents (node graph)
  "Get the dependents of @arg{node} in @arg{graph}."
  (gethash node graph))

(defclass parallel-heap (waitqueue-mixin)
  ((available-objects
    :initform (make-instance 'sorted-hash-table :insertion-order :heap))
   (blocked-objects :initform (make-hash-table))
   (done-objects :initform ())
   (dependents :initarg :dependents
               :documentation "A function returning all dependents of an object."))
  (:documentation "This class implements a heap for traversing a set of objects in parallel in such a way that the operation does not conflict with other threads working on the same set.  While the parallel operation works on some object, all dependents are blocked.  Furthermore, after the operation on an object is finished its dependents are priorized for being operated on for optimizing cache use."))

(defun make-parallel-heap (objects dependents)
  (lret ((heap (make-instance 'parallel-heap :dependents dependents)))
    (with-slots (available-objects) heap
      ;; Maybe: Change insertion order temporarily to :fifo?
      (loop for object in objects do
        (setf (dic-ref available-objects object) t)))))

(defun take-object (parallel-heap)
  (with-mutual-exclusion (parallel-heap)
    (with-slots (available-objects blocked-objects) parallel-heap
      (fl.parallel::wait
       parallel-heap
       :while (_ (dic-empty-p available-objects))
       :finish (_ (dic-empty-p blocked-objects))
       :perform
       (_ 
         (lret ((object (fl.dictionary::dic-pop available-objects)))
           ;; block neighborhood
           (with-slots (blocked-objects dependents) parallel-heap
             (let ((dependents (funcall dependents object)))
               (dbg :parallel-heap "Incrementing neighbors(~A) =~%~A~%"
                    object dependents)
               (loop for neighbor in dependents do
                 (cond ((dic-ref available-objects neighbor)
                        (setf (gethash neighbor blocked-objects) 1)
                        (dic-remove neighbor available-objects))
                       ((gethash neighbor blocked-objects)
                        (incf (gethash neighbor blocked-objects)))))))))))))

(defun drop-object (object parallel-heap)
  (with-mutual-exclusion (parallel-heap)
    (with-slots (available-objects blocked-objects done-objects dependents)
        parallel-heap
      ;; unblock neighborhood
      (let ((dependents (funcall dependents object)))
        (dbg :parallel-heap "Decrementing neighbors(~A) =~%~A~%"
             object dependents)
        (loop for neighbor in dependents do
          (whereas ((count (gethash neighbor blocked-objects)))
            (cond ((> count 1)
                   (setf (gethash neighbor blocked-objects)
                         (the fixnum (1- count))))
                  (t (remhash neighbor blocked-objects)
                     (setf (dic-ref available-objects neighbor) t))))))
      (push object done-objects))
    (fl.parallel::notify parallel-heap)))

(defgeneric popper (collection)
  (:documentation "Returns a function which accesses the given collection in parallel.
Calling this function yields a new argument list together with a second value
which is either NIL or a finalizer that should be called when working on
that argument list is finished.")
  (:method ((heap parallel-heap))
    (lambda ()
      (let ((object (take-object heap)))
        (values object
                (lambda () (drop-object object heap))))))
  (:method ((queue parqueue))
    (lambda () (dequeue queue))))

(defun multi-worker (popper worker)
  (loop (multiple-value-bind (object finalizer)
            (funcall popper)
          (unless object (return))
          (awhen worker (funcall it object))
          (awhen finalizer (funcall it)))))

(defmacro process-in-parallel (collection arglist &body body)
  (with-gensyms (popper worker)
  `(let ((,popper (popper ,collection)))
     (flet ((,worker ,arglist ,@body))
       (if lparallel::*kernel*
           (fl.parallel::broadcast-task (_ (multi-worker ,popper (function ,worker))))
           (multi-worker ,popper (function ,worker)))))))

;;;; Testing

(defun test-parallel-heap ()
  #+(or)
  (let* ((A (fl.matlisp:laplace-sparse-matrix 10))
         (heap (make-parallel-heap
                (range< 1 5) (curry #'fl.matlisp:matrix-row A))))
    (let ((o1 (take-object heap))
          (o2 (take-object heap)))
      (print o1)
      (print o2)
      (drop-object o1 heap)
      (drop-object o2 heap))
    (loop for object = (take-object heap) while object do
      (print object)
      (drop-object object heap)))
  ;;; this is OK for the iterative solver, not for the discretization!
  #+(or)
  (let ((mesh (uniform-mesh-on-box-domain (n-cube-domain 2) #(2 2))))
    (let ((heap (make-parallel-heap (cells-of-dim mesh nil)
                                    (_ (cdr (coerce (subcells _) 'list))))))
      (print (take-object heap))
      (print (take-object heap))
      (print (take-object heap))
      (print (take-object heap))
      heap
      ))

  
  )
