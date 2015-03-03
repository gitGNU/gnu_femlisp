(in-package :fl.parallel)

;;; Redefine memoization for making it thread safe

(defmacro with-memoization ((&key (type :global) size id (debug t) (test ''equal)) &body body)
  "Sets up a memoization environment consisting of a table, and a captured
symbol @symbol{memoizing-let} memoizing its body depending on the
arguments.  Example of usage:
@lisp
  (with-memoization (:type :local :id 'test)
    (defun test (n)
      (memoizing-let ((k (* 2 n)))
	(sleep 1)
	(* k k))))
@end lisp
If @arg{type} is :global, the table is thread-safe and the same for all
threads, if @arg{type} is :local, it is special for each thread."
  (when (eq type :local) (assert id))
  (with-gensyms (mutex table key key-exprs value-body func foundp value)
    `(let ((,mutex (bordeaux-threads:make-recursive-lock))
	   (,table ,(ecase type
                           (:global `(fl.utilities::memoization-table ,size ,test))
                           (:local
                            `(ensure (getf *thread-local-memoization-table* ,id)
                                     (fl.utilities::memoization-table ,size ,test))))))
	;; the memoizing-let symbol is captured, 
       (macrolet ((memoizing (,func)
                    `(lambda (&rest ,',key)
		       (bordeaux-threads:with-recursive-lock-held (,',mutex)
			 (multiple-value-bind (,',value ,',foundp)
			     (dic-ref ,',table ,',key)
			   (if ,',foundp
			       ,',value
			       (progn
				 ,',(when debug `(dbg :memoize "Memoizing: id=~A, key=~A" ,id ,key))
				 (setf (dic-ref ,',table ,',key)
				       (apply ,,func ,',key))))))))
                  (memoizing-let (,key-exprs &body ,value-body)
                    (let ((,key-exprs
                           (mapcar (lambda (entry)
                                     (if (symbolp entry)
                                         (list entry entry)
                                         entry))
                                   ,key-exprs)))
                      `(let ,,key-exprs
                         ;; declares should be inserted here from value-body
                         (bordeaux-threads:with-recursive-lock-held (,',mutex)
                           (let ((,',key (list ,@(mapcar #'car ,key-exprs))))
                             (multiple-value-bind (,',value ,',foundp)
                                 (dic-ref ,',table ,',key)
                               (if ,',foundp
                                   ,',value
                                   (progn
                                     ,',(when debug `(dbg :memoize "Memoizing: id=~A, key=~A" ,id ,key))
                                     (setf (dic-ref ,',table ,',key)
                                           (progn ,@,value-body)))))))))))
         ,@body))))

;;; replace
(defmacro with-workers ((work) &body body)
  "This macro distributes work generated in body with calling the locally
bound function @function{work-on} on some arguments to several working
threads which call @arg{func} on those arguments."
  (with-gensyms (channel worker count)
    `(let ((,worker ,work)
           (,channel (lparallel:make-channel))
           (,count 0))
       (flet ((work-on (&rest args)
                (apply #'lparallel:submit-task ,channel ,worker args)
                (incf ,count)))
	 ,@body
         (loop repeat ,count do (lparallel:receive-result ,channel))))))

;;; replace

(let ((lock (bordeaux-threads:make-recursive-lock)))
  (defmethod dbg (id format-string &rest args)
    (dbg-when id
      (bordeaux-threads:with-recursive-lock-held (lock)
        (format *debug-io* "~&~?~%" format-string args)
        (force-output *debug-io*)))))

(defun test-parallel-adaptions ()
  (let* ((count 0)
         (f (with-memoization (:type :global :id 'test)
              (memoizing (lambda (k)
                           (sleep 1.0)
                           (incf count)
                           (* k k))))))
    (pwork f #((1) (1)))
    (time (funcall f 1))
    count)
  
  (time
   (with-workers ((lambda (k) (* k k)))
     (loop for i below 100000 do (work-on i))))
  )

  
