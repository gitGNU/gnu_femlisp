(in-package :fl.parallel)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Very simple wrapping of a lock around some data
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass mutex-object ()
  ((lock :initform (bordeaux-threads:make-recursive-lock "my-lock"))
   (data :initarg :data)))

(defun mutex-wrap (object)
  (make-instance 'mutex-object :data object))

(defmacro accessing-exclusively (bindings &body body)
  (if bindings
      (destructuring-bind ((name object) . rest)
          bindings
        `(let ((,name ,object))
           (bordeaux-threads:with-recursive-lock-held ((slot-value ,name 'lock))
             (with-slots ((,name data)) ,name
               (accessing-exclusively ,rest ,@body)))))
      `(locally ,@body)))

(defmacro with-mutex (object &body body)
  "Needs to be called unfortunately around push/incf/setf."
  `(bordeaux-threads:with-recursive-lock-held ((slot-value ,object 'lock))
     ,@body))

