(in-package :fl.parallel)

;;; replace
(defmacro with-workers ((work) &body body)
  "This macro distributes work generated in body with calling the locally
bound function @function{work-on} on some arguments to several working
threads which call @arg{func} on those arguments."
  (with-gensyms (channel worker count)
    `(if *kernel*
         (let ((,worker ,work)
               (,channel (lparallel:make-channel))
               (,count 0))
           (flet ((work-on (&rest args)
                    (apply #'lparallel:submit-task ,channel
                           (lambda (&rest args) (apply ,worker args) nil)
                           args)
                    (incf ,count)))
             ,@body
             (loop repeat ,count do (lparallel:receive-result ,channel))))
         (let ((,worker ,work))
           (flet ((work-on (&rest args)
                    (apply ,worker args)))
             ,@body)))))

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
  )

  
