(in-package :fl.discretization)

(defun make-forest (objects get-children &optional result)
  "Returns a forest containing OBJECTS for the relation GET-CHILDREN.
The result is a hash-table mapping objects to their children.  E.g., for the
boundary relation the hash-table maps n-dimensional cells to their
n-1-dimensional sides."
  (lret ((result (or result (make-hash-table))))
    (dolist (object objects)
      (let ((children (coerce (funcall get-children object) 'list)))
        (acond ((gethash object result)
                (assert (set-equal children it)))
               (t (setf (gethash object result)
                        children)
                  (make-forest children get-children result)))))))

(defun invert-forest (forest)
  "Invert the forest.  The result is a hash-table mapping a child to its parents.  In the case of a forest derived by the boundary relation this results in a mapping of n-dimensional cells to their n+1-dimensional coboundaries."
  (lret ((result (make-hash-table)))
    (maphash
     (lambda (object children)
       (dolist (child children)
         (pushnew object (gethash child result))))
     forest)))

(defun get-branches (tree forest)
  (gethash tree forest))

(defclass assembly-heap (waitqueue-mixin)
  ((available-cells
    :initform (make-instance 'sorted-hash-table :insertion-order :heap))
   (blocked-cells :initform (make-hash-table))
   (done-cells :initform ())
   (children :initarg :children)
   (parents :initarg :parents)))

(defun assembly-heap (cells)
  (lret* ((children (make-forest cells #'boundary))
          (parents (invert-forest children))
          (heap (make-instance 'assembly-heap
                               :children children
                               :parents parents)))
    (with-slots (available-cells) heap
      ;; Maybe: Change insertion order temporarily to :fifo?
      (loop for cell in cells do
        (setf (dic-ref available-cells cell) t)))))

(defun mark-neighborhood (object mark forest neighborhood)
  (let ((object-mark (gethash object neighborhood)))
    (unless object-mark
      (setf (gethash object neighborhood) mark))
    (unless (eql object-mark mark)
      (dolist (child (get-branches object forest))
        (mark-neighborhood child mark forest neighborhood)))))

(defun recursive-neighborhood (cells mark forest &optional (neighborhood (make-hash-table)))
  (dolist (cell cells)
    (mark-neighborhood cell mark forest neighborhood))
  neighborhood)
  
(defun take-cell (assembly-heap)
  (with-mutual-exclusion (assembly-heap)
    (with-slots (available-cells blocked-cells) assembly-heap
      (fl.parallel::wait
       assembly-heap
       :while (_ (dic-empty-p available-cells))
       :finish (_ (dic-empty-p blocked-cells))
       :perform
       (_ 
         (lret ((cell (fl.dictionary::dic-pop available-cells)))
           (with-slots (blocked-cells children parents) assembly-heap
             (let ((neighborhood-1 (recursive-neighborhood (list cell) 1 children)))
               (dic-for-each-key
                (lambda (neighbor)
                  (when (dic-ref available-cells neighbor)
                    (pushnew cell (gethash neighbor blocked-cells ()))
                    (dic-remove neighbor available-cells)))
                (recursive-neighborhood (keys neighborhood-1) 2 parents neighborhood-1))))))))))

(defun drop-cell (cell assembly-heap)
  (with-mutual-exclusion (assembly-heap)
    (with-slots (available-cells blocked-cells done-cells children parents)
        assembly-heap
      ;; unblock union of neighborhood-1 and neighborhood-2
      (let ((neighborhood-1
              (recursive-neighborhood (list cell) 1 children)))
          (dic-for-each-key
           (lambda (neighbor)
             (awhen (gethash neighbor blocked-cells)
               (let ((blocks (delete cell it)))
                 (cond (blocks (setf (gethash neighbor blocked-cells) blocks))
                       (t (remhash neighbor blocked-cells)
                          (setf (dic-ref available-cells neighbor) t))))))
           (recursive-neighborhood (keys neighborhood-1) 2 parents neighborhood-1)))
      (push cell done-cells))
    (fl.parallel::notify assembly-heap)))

(defgeneric popper (collection)
  (:documentation "Returns a function which accesses the given collection in parallel.
Calling this function yields a new argument list together with a second value
which is either NIL or a finalizer that should be called when working on
that argument list is finished.")
  (:method ((heap assembly-heap))
    (lambda ()
      (let ((cell (take-cell heap)))
        (values cell
                (lambda () (drop-cell cell heap))))))
  (:method ((queue parqueue))
    (lambda () (dequeue queue))))

(defun multi-worker (popper worker)
  (loop (multiple-value-bind (object finalizer)
            (funcall popper)
          (unless object (return))
          (awhen worker (funcall it object))
          (awhen finalizer (funcall it)))))

(defmacro parallel-processing (collection arglist &body body)
  (with-gensyms (popper worker)
  `(let ((,popper (popper ,collection)))
     (flet ((,worker ,arglist ,@body))
       (if lparallel::*kernel*
           (fl.parallel::broadcast-task (_ (multi-worker ,popper #',worker)))
           (multi-worker ,popper #',worker))))))

;;;; Testing

(defun test-assembly-heap ()
  (let* ((dim 2) (level 1)
         (mesh (uniformly-refined-mesh (n-cube-domain dim) level))
         (cells (cells-of-dim mesh dim)))
    (let ((heap (assembly-heap cells)))
      (loop for cell = (take-cell heap) while cell do
        (print cell)
        (drop-cell cell heap))))

  (time
   (let* ((dim 2) (level 2)
          (mesh (uniformly-refined-mesh (n-cube-domain dim) level))
          (cells (cells-of-dim mesh dim)))
     (let ((heap (assembly-heap cells)))
       (multi-worker (popper heap)
                     (_ (sleep 0.5) (print _))))))

  (time
   (let* (;lparallel:*kernel*
          (dim 2) (level 2)
          (mesh (uniformly-refined-mesh (n-cube-domain dim) level))
          (cells (cells-of-dim mesh dim)))
     (parallel-processing (assembly-heap cells)
         (cell)
       (sleep 0.5)
       (print cell))))
  )
