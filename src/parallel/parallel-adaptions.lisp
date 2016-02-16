(in-package :fl.parallel)

(defun call-with-workers (worker producer &key (parallel t))
  (if (and parallel *kernel* (not *worker-id*))
      (let ((channel (lparallel:make-channel))
            (count 0))
        (flet ((work-on (&rest args)
                 (apply #'lparallel:submit-task channel
                        (lambda (&rest args) (apply worker args) nil)
                        args)
                 (incf count)))
          (funcall producer #'work-on))
        (loop repeat count do (lparallel:receive-result channel)))
      (funcall producer worker)))

(defmacro with-workers ((work &rest keyword-args &key &allow-other-keys)
                        &body body)
  "This macro distributes work generated in body with calling the locally
bound function @function{work-on} on some arguments to several working
threads which call @arg{func} on those arguments."
  (with-gensyms (worker)
    `(call-with-workers
      ,work
      (lambda (,worker)
        (flet ((work-on (&rest args)
                 (apply ,worker args)))
          ,@body))
      ,@keyword-args)))

(defmacro with-accumulators ((name initial-element reduce-op) &body body)
  "This macro sets up an array of accumulators for each worker initialized
with initial-element.  Within each worker the symbol @arg{name} can be used
for accessing the private accumulator.
After the body was executed (probably containing parallel operations),
the result array is reduced using @arg{reduce-op}."
  (with-gensyms (accumulators)
  `(if (and *kernel* (not *worker-id*))
       (let ((,accumulators (make-array (kernel-worker-count) :initial-element ,initial-element)))
         (symbol-macrolet ((,name (aref ,accumulators *worker-id*)))
           ,@body)
         (reduce ,reduce-op ,accumulators))
       (lret ((,name ,initial-element)) ,@body))))

;;; replace
(let ((lock (bordeaux-threads:make-recursive-lock)))
  (defmethod dbg (id format-string &rest args)
    (dbg-when id
      (bordeaux-threads:with-recursive-lock-held (lock)
        (format *debug-io* "~&~?~%" format-string args)
        (force-output *debug-io*)))))

(defun test-parallel-adaptions ()
  (time
   (with-workers ((lambda (k) (print (* k k))))
     (loop for i below 10 do (work-on i))))
  (loop for i below 1000 do
    (assert (= 45
               (with-accumulators (sum 0 #'+)
                 (with-workers ((lambda (x)
                                  (incf sum x)))
                   (loop for i below 10 do (work-on i)))))))
  (with-accumulators (sum 0 #'+)
    (pwork (lambda () (incf sum 1))))
  
  (let ((x #(1 2 3))
        (y #(2 3 4)))
    (with-accumulators (sum 0 #'+)
      (with-workers
          ((lambda (xc i)
             (incf sum (* xc (aref y i)))))
        (loop for i from 0 and xc across x do
          (work-on xc i)))))
  (dbg-on :mp)
  (with-workers ((lambda (k)
                   (mp-dbg "Working on ~D~%" k)))
    (work-on 1)
    (work-on 2))
  (call-with-workers (lambda () (print fl.parallel::*worker-id*))
                     (lambda (worker) (funcall worker)))

  (call-with-workers (lambda (k) (mp-dbg "Work on ~D" k))
                     (lambda (worker)
                       (flet ((work-on (&rest args)
                                (apply worker args)))
                         (work-on 1)
                         (work-on 2))))

  (with-workers ((lambda (k) (mp-dbg "Work on ~D" k)))
    (Work-on 1)
    (work-on 2))
  (with-workers ((lambda (k)
                   (mp-dbg "id=~A k=~D" *worker-id* k)
                   (with-workers ((lambda (l)
                                    (mp-dbg "Working on ~D and ~D ~%" k l)))
                     (work-on 3)
                     (work-on 4))))
    (work-on 1)
    (work-on 2))
  )

  
