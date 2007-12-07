;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; multiprocessing.lisp - routines for multiprocessing
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;
;;; Copyright (C) 2003-2006 Nicolas Neuss, University of Heidelberg.
;;; All rights reserved.
;;; 
;;; Redistribution and use in source and binary forms, with or without
;;; modification, are permitted provided that the following conditions are
;;; met:
;;; 
;;; 1. Redistributions of source code must retain the above copyright
;;; notice, this list of conditions and the following disclaimer.
;;; 
;;; 2. Redistributions in binary form must reproduce the above copyright
;;; notice, this list of conditions and the following disclaimer in the
;;; documentation and/or other materials provided with the distribution.
;;; 
;;; THIS SOFTWARE IS PROVIDED ``AS IS'' AND ANY EXPRESS OR IMPLIED
;;; WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
;;; MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN
;;; NO EVENT SHALL THE AUTHOR, THE UNIVERSITY OF HEIDELBERG OR OTHER
;;; CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
;;; EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
;;; PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
;;; PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
;;; LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
;;; NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
;;; SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defpackage "FL.MULTIPROCESSING"
  (:nicknames "FL.MP")
  (:use "COMMON-LISP" "FL.UTILITIES" "FL.MACROS" "FL.PORT" "FL.DEBUG")
  (:export
   ;; threading
   "MAKE-THREAD" "MAKE-MUTEX" "WITH-MUTEX"
   "MAKE-WAITQUEUE" "CONDITION-WAIT" "CONDITION-NOTIFY"
   
   "MUTEX-MIXIN" "WITH-LOCK"
   "WITH-MUTUAL-EXCLUSION"

   "LOCKED-REGION-MIXIN"
   "WITH-REGION-DO" "WITH-REGION"
   "MPQUEUE" "PARQUEUE"
   "*THREAD-LOCAL-MEMOIZATION-TABLE*"
   "*NUMBER-OF-THREADS*" "*CHUNK-SIZE*"
   "WITH-WORKERS" "WORK-ON")
  (:documentation
   "This package provides an interface for allowing parallel execution in
Femlisp.  It abstracts also from the underlying Lisp's threading features.
Unfortunately, these features somewhat differ (e.g. the waitqueue concept
in SBCL and gates in Allegro), although one can usually achieve the same
effects."))

(in-package "FL.MULTIPROCESSING")

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; A implementation-independent interface
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; threads

(defun make-thread (func &key name initial-bindings)
  "Execute @arg{func} as a thread with name @arg{name} within dynamic
bindings given by @arg{initial-bindings}."
  (when initial-bindings
    (warn "Not supported by SBCL."))
  #+allegro (mp:process-run-function
	     (list :name name :initial-bindings initial-bindings)
	     func)
  #+cmu (mp:make-process func :name name :initial-bindings initial-bindings)
  #+sb-thread (sb-thread:make-thread func :name name)
  #-(or allegro cmu sb-thread) nil
  )

(defun thread-yield (&optional absolutely)
  "Allows other processes to run.  If @arg{absolutely} is T, then other
processes have to run."
  (declare (ignorable absolutely))
  #+allegro (mp:process-allow-schedule)
  #+cmu (mp:process-yield)
  #-(or allegro cmu) nil
  )

(defun current-thread ()
  "Returns the current thread."
  #+allegro system:*current-process*
  #+cmu (mp:current-process)
  #+sb-thread sb-thread:*current-thread*
  #-(or allegro cmu sb-thread) (list :main-thread)
  )

(defun thread-name (thread)
  "Returns the name of the current thread."
  #+allegro (mp:process-name thread)
  #+cmu (mp:process-name thread)
  #+sb-thread (sb-thread:thread-name thread)
  #-(or allegro cmu sb-thread) nil
  )

(defun list-all-threads ()
  "Returns a list of all threads"
  #+allegro mp:*all-processes*
  #+cmu (mp:all-processes)
  #+sb-thread (sb-thread:list-all-threads)
  #-(or allegro cmu sb-thread) (list :main-thread)
  )

(defun terminate-thread (thread)
  #+allegro (mp:process-kill thread)
  #+cmu (mp:destroy-process thread)
  #+sb-thread (sb-thread:terminate-thread thread)
  #-(or allegro cmu sb-thread)
  ;; no threading: thus only the main-thread which must not be killed
  (assert (null thread))
  )
  
;;; mutices

(defun make-mutex (&key (name "mutex"))
  (declare (ignorable name))
  #+allegro (mp:make-process-lock :name name)
  #+cmu (mp:make-lock name)
  #+sb-thread (sb-thread:make-mutex :name name)
  #-(or allegro cmu sb-thread) nil
  )

(defmacro with-mutex ((mutex) &body body)
  (declare (ignorable mutex))
  #+allegro `(mp:with-process-lock (,mutex) ,@body)
  #+cmu `(mp:with-lock-held (,mutex) ,@body)
  #+sb-thread `(sb-thread:with-recursive-lock (,mutex) ,@body)
  #-(or allegro cmu sb-thread) `(locally ,@body))

;;; waitqueues (following the SBCL interface)

(defun make-waitqueue ()
  "Generates a waitqueue."
  #+allegro (mp:make-gate nil)
  #+sb-thread (sb-thread:make-waitqueue)
  #-(or allegro sb-thread) nil 
  )

(defun condition-wait (waitqueue mutex)
  "Registers on the waitqueue, releases mutex, and waits for a notification
on the waitqueue."
  (declare (ignorable waitqueue mutex))
  #+allegro
  (progn
    (mp:close-gate waitqueue)
    (mp:process-unlock mutex)
    (mp:process-wait "waiting for gate" #'mp:gate-open-p waitqueue)
    (mp:process-lock mutex))
  #+sb-thread (sb-thread:condition-wait waitqueue mutex)
  #-(or allegro sb-thread) (thread-yield)
  )

(defun condition-notify (waitqueue &optional (n 1))
  "Notifies on the waitqueue."
  (declare (ignorable waitqueue n))
  #+allegro (mp:open-gate waitqueue)
  #+sb-thread
  (if (eq n t)
      (sb-thread:condition-broadcast waitqueue)
      (sb-thread:condition-notify waitqueue n))
  #-(or allegro sb-thread) (thread-yield)
  )

;;; nicer OO interface

(defclass mutex-mixin ()
  ((mutex :reader mutex :initform (make-mutex) :documentation
	  "A mutex for excluding access of other threads."))
  (:documentation "A mixin which adds a mutex to every instance of the
class."))

(defmacro with-mutual-exclusion ((obj) &body body)
  "Execute @arg{body} on the waitqueue @arg{obj} without other threads
interfering."
  `(with-mutex ((mutex ,obj))
    ,@body))

(defclass waitqueue-mixin (mutex-mixin)
  ((waitqueue :reader waitqueue :initform (make-waitqueue)))
  (:documentation "Waitqueue mixin."))

(defgeneric wait (waitqueue-object &key while until perform)
  (:documentation "Waits on @arg{waitqueue-object} while @arg{while} is
satisfied or until @arg{until} is satisfied.  After this, if it is given,
the function @perform is called.")
  (:method ((wq waitqueue-mixin) &key while until perform)
    (with-mutual-exclusion (wq)
      (loop while (or (and while (funcall while wq))
		      (and until (not (funcall until wq))))
	    do (condition-wait (waitqueue wq) (mutex wq)))
      (awhen perform (funcall it)))))

(defgeneric notify (waitqueue-object &optional n)
  (:documentation "Notifies @arg{waitqueue-object}.")
  (:method ((wq waitqueue-mixin) &optional (n t))
    ;;(mp-dbg "notifying ~A" wq)
    (condition-notify (waitqueue wq) n)))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; locked-region
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass locked-region-mixin (waitqueue-mixin)
  ((locked-region :reader locked-region :initform (make-hash-table))))

(defmethod lock-region ((object locked-region-mixin) keys)
  "Should not be used externally.  Use WITH-REGION instead."
  (with-mutual-exclusion (object)
    (let ((table (locked-region object)))
      (loop+ ((key keys)) do
	 (awhen (gethash key table)
	   (assert (eq it (current-thread))))
	 (setf (gethash key table)
	       (current-thread))))))

(defmethod unlock-region ((object locked-region-mixin) keys)
  (with-mutual-exclusion (object)
    (let ((table (locked-region object)))
      (loop+ ((key keys)) do
	 (remhash key table)))))

(defparameter *count* (cons 0 0))
(defmethod with-region-do ((object locked-region-mixin) keys perform)
  (with-mutual-exclusion (object)
    (wait object
	  :until
	  (lambda (object)
	    (let ((table (locked-region object)))
	      (lret ((result (notany (lambda (key)
				       (aand (gethash key table)
					     (not (eq it (current-thread)))))
				     keys)))
		(if result
		    (incf (car *count*))
		    (incf (cdr *count*)))))))
    (lock-region object keys))
  (funcall perform)
  (unlock-region object keys)
  (notify object))

(defmacro with-region ((object keys) &body body)
  `(with-region-do ,object ,keys
		   (lambda () ,@body)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Multithreaded debugging
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defvar *debug-mutex* (make-mutex))

(defmethod dbg :around (id format-string &rest args)
  "Serializing of debug output."
  (declare (ignore id format-string args))
  (with-mutex (*debug-mutex*)
    (call-next-method)))

(defun mp-dbg (format &rest args)
  (funcall #'dbg :mp "~A:~%~?" (current-thread) format args))

#+(or)
(format t "~A:~%~?" 'THREAD "working on ~A" '((1 2 3)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass multithread (waitqueue-mixin)
  ((threads :accessor threads :initform ())
   (work :initarg :work :documentation "Function doing the actual work.")))

(defmethod initialize-instance :after ((mt multithread) &key number-of-threads &allow-other-keys)
  (with-slots (waitqueue mutex threads work) mt
    (with-mutual-exclusion (mt)
      (setf threads
	    (loop repeat number-of-threads collect
		  (make-thread
		   ;; worker function
		   (lambda ()
		     (funcall work)
		     ;; epilogue after doing the work
		     (with-mutual-exclusion (mt)
		       (setf threads (remove (current-thread) threads))
		       (notify mt)))
		   :name "my-worker"))))))

#+(or)
(defun execute-in-parallel (work &key (number-of-threads *number-of-threads*))
  (let ((mt (make-instance 'multithread :work work
			   :number-of-threads number-of-threads)))
    (wait mt :while #'threads)))

#+(or)
(loop for i from 1 upto 5 do
      (format t "~R thread~:P~%" i)
      (time (execute-in-parallel
	     (lambda ()
	       (loop repeat (expt 2 27) sum 1))
	     :number-of-threads i)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; thread-safe queue
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass parqueue (queue waitqueue-mixin)
  ((finished-p :accessor finished-p :initform nil :documentation
	       "Indicator for no-more-input-allowed."))
  (:documentation "A thread-safe queue waiting for input."))

(defmethod enqueue :around (obj (pq parqueue))
  (assert (not (finished-p pq)))
  (mp-dbg "enqueueing ~A in ~A" obj pq)
  (with-mutual-exclusion (pq)
    (call-next-method)
    (notify pq) ; notify of change!
    nil
    ))

(defmethod dequeue :around ((pq parqueue))
  (with-mutual-exclusion (pq)
    (loop
     (mp-dbg "trying to dequeue ~A (finished: ~A, empty: ~A, top: ~A)"
	     pq (finished-p pq) (emptyp pq) (car (fl.utilities::head pq)))
     (cond ((or (finished-p pq) (not (emptyp pq)))
	    (multiple-value-bind (value emptyp)
		(call-next-method)
	      (notify pq)			; notify of change
	      (mp-dbg "dequeued value=~A, emptyp=~A" value emptyp)
	      (return-from dequeue (values value emptyp))))
	   (t (mp-dbg "waiting on parqueue ~A" pq)
	      (condition-wait (waitqueue pq) (mutex pq))))
     )))

(defmethod finish ((pq parqueue))
  (with-mutual-exclusion (pq)
    (mp-dbg "finishing queue ~A" pq)
    (setf (finished-p pq) t)
    (notify pq)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; work-group
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass work-group (waitqueue-mixin)
  ((threads :accessor threads :initform ())
   (tasks :reader tasks :initform (make-instance 'parqueue))
   (work :reader work :initarg :work))
  (:documentation "A class representing a group of threads working on a
queue."))

(defparameter *thread-local-memoization-table* nil
  "Dynamic variable specialized for each worker thread separately.  Should
not be used globally.")

(defmethod add-worker ((group work-group))
  (with-mutual-exclusion (group)
    (push (make-thread
	   (lambda ()
	     (mp-dbg "starting...")
	     (let ((*thread-local-memoization-table* ()))
	       (loop
		(multiple-value-bind (args emptyp)
		    (dequeue (tasks group))
		  (when emptyp (return))
		  (mp-dbg "working on ~A" args)
		  (apply (work group) args)
		  (mp-dbg "finished working on ~A." args)
		  (thread-yield)))
	       (mp-dbg "...done")
	       (remove-worker group (current-thread))))
	   :name "femlisp-worker")
	  (threads group))))

#+(or)
(defmethod add-worker ((group work-group) &optional multiple-args)
  (with-mutual-exclusion (group)
    (push (make-thread
	   (lambda ()
	     (mp-dbg "starting...")
	     (let ((*thread-local-memoization-table* ()))
	       (loop
		(multiple-value-bind (args emptyp)
		    (dequeue (tasks group))
		  (when emptyp (return))
		  (mp-dbg "working on ~A" args)
		  (if multiple-args
		      (dolist (args args)
			(apply (work group) args))
		      (apply (work group) args))
		  (mp-dbg "finished working on ~A." args)
		  (thread-yield)))
	       (mp-dbg "...done")
	       (remove-worker group (current-thread))))
	   :name "femlisp-worker")
	  (threads group))))

(defmethod remove-worker ((group work-group) thread)
  "Removes @arg{thread} from @arg{group}."
  (assert (finished-p (tasks group)))
  (with-mutual-exclusion (group)
    (assert (member thread (threads group)))
    (mp-dbg " --- removes itself.")
    (setf (threads group)
	  (delete thread (threads group)))
    (notify group)))

(defmethod send-task ((group work-group) task)
  "Send a task to the work-group."
  (enqueue task (tasks group)))

(defvar *number-of-threads* nil  ; (optimal-thread-number)
  "The number of threads in which Femlisp tries to split the work for some
computationally intensive tasks.  If NIL, no threading is used.")

;; (setq *number-of-threads* 2)

(defvar *chunk-size* nil
  "NIL means no argument grouping, a number means that the arguments are
  grouped in chunks of this size.")

(defun execute-in-parallel (work &key arguments distribute
			    (number-of-threads *number-of-threads*))
  "Executes work in parallel distributed on @arg{number-of-threads}
threads.  The arguments for work are taken from the argument-lists in
@arg{arguments}.  Alternatively, @arg{distribute} may contain a function
which generates the argument-lists by calling the function
@function{work-on}.  If group-size is non-nil, it is used as a chunk-size.
This is important, if the amount of work on each argument list is not
large."
  (if number-of-threads
      ;; parallel execution
      (let ((group (make-instance 'work-group :work work)))
	(loop repeat number-of-threads do (add-worker group))
	(if distribute
	    (funcall distribute (lambda (&rest arglist) (send-task group arglist)))
	    (loop for arglist in arguments do (send-task group arglist)))
	(finish (tasks group))
	(wait group :while #'threads))
      ;; sequential execution
      (if distribute
	  (funcall distribute work)
	  (loop for arglist in arguments
	       do (apply work arglist)))))

#+(or)
(defun execute-in-parallel (work &key arguments distribute
			    (number-of-threads *number-of-threads*)
			    (chunk-size *chunk-size*))
  "Executes work in parallel distributed on @arg{number-of-threads}
threads.  The arguments for work are taken from the argument-lists in
@arg{arguments}.  Alternatively, @arg{distribute} may contain a function
which generates the argument-lists by calling the function
@function{work-on}.  If chunk-size is non-nil, it is used as a size.  This
is important, if the amount of work on each argument list is not large."0
  (if number-of-threads
      ;; parallel execution
      (let ((group (make-instance 'work-group :work work)))
	;; generate worker threads
	(loop repeat number-of-threads do
	     (add-worker group chunk-size))
	;; generate work
	(let ((chunk ())
	      (current-size 0))
	  (flet ((collect-and-send (&rest args)
		   (cond ((null chunk-size)
			  (send-task group args))
			 ((= current-size chunk-size)
			  (send-task group chunk)
			  (setf chunk () current-size 0))
			 (t (push args chunk)
			    (incf current-size)))))
	    (if distribute
		(funcall distribute #'collect-and-send)
		(mapc #'collect-and-send arguments))
	    (when chunk
	      (send-task group chunk))))
	(finish (tasks group))
	(wait group :while #'threads))
      ;; sequential execution
      (if distribute
	  (funcall distribute work)
	  (loop for arglist in arguments
	       do (apply work arglist)))))

(defun femlisp-workers ()
  (remove "femlisp-worker" (list-all-threads)
	  :test-not #'string= :key #'thread-name))

(defun terminate-workers ()
  (dolist (thread (femlisp-workers))
    (terminate-thread thread)))

(defmacro with-workers ((work) &body body)
  "This macro distributes work generated in body with calling the locally
bound function @function{work-on} on some arguments to several working
threads which call @arg{func} on those arguments."
  (with-gensyms (send-work)
    `(execute-in-parallel
      ,work
      :distribute
      (lambda (,send-work)
	(flet ((work-on (&rest args) (apply ,send-work args)))
	 ,@body)))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Find optimal threading parameter for this architecture
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun optimal-thread-number ()
  (let ((count
	 (nth-value 1 (fl.utilities::measure-time-repeated
		       (lambda () (loop repeat 10000 sum 1))
		       0.5))))
    (loop for nr-threads from 1 upto 8
	  for previous-time = nil then time
	  for time =
	  (fl.utilities::measure-time
	   (lambda ()
	     (execute-in-parallel
	      (lambda () (loop repeat (* 10000 count) sum 1))
	      :number-of-threads nr-threads
	      :arguments (make-list nr-threads :initial-element ())))
	   1 t)
	  until (and previous-time (> time (* previous-time 1.1)))
	  finally (return (1- nr-threads)))))

;;(optimal-thread-number)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Thread-global and thread-local cacheing of results
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; in general.lisp

;;; Testing
(defun femlisp-multiprocessing-tests ()

  (dbg-on :mp)
  (let ((mutex (make-mutex))
	(result ())
	(pq (make-instance 'parqueue)))
    (flet ((worker ()
	     (loop for obj = (dequeue pq) while obj
		   do (with-mutex (mutex)
			(push obj result)))
	     (format t "~A done.~%" (current-thread))))
      (make-thread #'worker :name "femlisp-worker")
      (make-thread #'worker :name "femlisp-worker")
      (make-thread #'worker :name "femlisp-worker"))
    (loop for k from 1 upto 10 do
	  (enqueue k pq))
    (finish pq)
    (sleep 0.5)
    result)

  (assert (null (femlisp-workers)) () "workers remained - A")
  ;; (terminate-workers)

  (time
   (let ((n 10000)
	 (q (make-instance 'parqueue)))
     (let ((*number-of-threads* 4))
       (with-workers ((lambda (k)
			(enqueue k q)))
	 (loop for k below n do (work-on k))))
     :done))
  
  (sleep 0.5)
  (assert (null (femlisp-workers)) () "workers remained - B")
  ;; (terminate-workers)

  (let ((k 100000)
	(mutex (make-mutex))
	(sum 0))
    (declare (ignorable mutex))
    (time (print (loop repeat k summing 1)))
    (time
     (let ((*number-of-threads* 4))
       (with-workers
	   ((lambda (k)
	      (declare (type fixnum k))
	      (let ((s (loop repeat k summing 1)))
		(with-mutex (mutex)
		  (incf sum s)))))
	 (loop repeat *number-of-threads* do (work-on k)))
       (format t "Sum: ~A" sum))))
  
  (sleep 0.5)
  (assert (null (femlisp-workers)) () "workers remained - C")
  ;;; (terminate-workers)

  ;; no tasks for 3 threads
  (time (execute-in-parallel
	 (lambda () (error "should not be called"))
	 :arguments ()
	 :number-of-threads 3))

  (sleep 0.5)
  (assert (null (femlisp-workers)) () "workers remained - D")
  ;; (terminate-workers)

  (let ((n (expt 2 27)))
    (loop for i from 1 upto 5 do
	  (format t "~R thread~:P~%" i)
	  (time (execute-in-parallel
		 (lambda (k)
		   (loop repeat k sum 1))
		 :arguments (make-list i :initial-element (list n))
		 :number-of-threads i))))

  (sleep 0.5)
  (assert (null (femlisp-workers)) () "workers remained - E")
  ;; (terminate-workers)

  (let ((n (expt 2 27)))
    (loop for i from 1 upto 5 do
	 (format t "~R thread~:P~%" i)
	 (let ((*number-of-threads* i))
	   (time (with-workers ((lambda (k) (loop repeat k sum 1)))
		   (loop repeat i do (work-on n)))))))

  (sleep 0.5)
  (assert (null (femlisp-workers)) () "workers remained - F")
  ;; (list-all-threads)
  ;; (terminate-workers)
  
  #+(or)
  (let ((x (make-instance 'mrsw-mixin)))
    (with-read-access (x)
      (print 'hi))
    (with-readwrite-access (x)
      (print 'ho)))
  
  (dbg-off :mp)
  )

;;; (femlisp-multiprocessing-tests)