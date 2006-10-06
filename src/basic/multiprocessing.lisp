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

   "*THREADS*"
   "MUTEX-MIXIN" "WITH-LOCK"
   "WITH-MUTUAL-EXCLUSION"
   "MPQUEUE" "PARQUEUE"
   "*THREAD-LOCAL-MEMOIZATION-TABLE*"
   "WITH-WORKERS" "WORK-ON")
  (:documentation
   "This package provides an interface for allowing parallel execution in
Femlisp.  It abstracts also from the underlying Lisp's threading features.
Unfortunately, these features somewhat differ (e.g. the waitqueue concept
in SBCL and gates in Allegro), although one can usually achieve the same
effects."))

(in-package "FL.MULTIPROCESSING")

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Reaching a common interface
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Threads
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(eval-when (:compile-toplevel :load-toplevel :execute)
  (defvar *threads* nil
    "The number of threads in which Femlisp tries to split the work for some
computationally intensive tasks.  If NIL, no threading is used.")
  )


(defun make-thread (func &key name initial-bindings)
  "Execute @arg{func} as a thread with name @arg{name} within dynamic
bindings given by @arg{initial-bindings}."
  (declare (ignorable name initial-bindings))
  #+allegro (mp:process-run-function
	     `(:name ,name :initial-bindings ,initial-bindings)
	     func)
  #+cmu (mp:make-process func :name name :initial-bindings initial-bindings)
  #+sb-thread (sb-thread:make-thread func :name name)
  )

(defun thread-yield (&optional absolutely)
  "Allows other processes to run.  If @arg{absolutely} is T, then other
processes have to run."
  (declare (ignorable absolutely))
  #+allegro (mp:process-allow-schedule)
  #+cmu (mp:process-yield)
  )

(defun current-thread ()
  "Returns the current thread."
  #+allegro system:*current-process*
  #+cmu (mp:current-process)
  #+sb-thread sb-thread:*current-thread*
  )
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Mutex
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun make-mutex (&key (name "mutex"))
  (declare (ignorable name))
  #+allegro (mp:make-process-lock :name name)
  #+cmu (mp:make-lock name)
  #+sb-thread (sb-thread:make-mutex :name name)
  )

(defmacro with-mutex ((mutex) &body body)
  (declare (ignorable mutex))
  (if *threads*
      (list 'progn
	    #+allegro `(mp:with-process-lock (,mutex) ,@body)
	    #+cmu `(mp:with-lock-held (,mutex) ,@body)
	    #+sb-thread `(sb-thread:with-mutex (,mutex) ,@body))
      `(progn ,@body)))

;;; Mutex mixin - should maybe be a metaclass?

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

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Waitqueue
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; SBCL-specific performance enhancement

(defun make-waitqueue ()
  "Generates a waitqueue."
  #+sb-thread (sb-thread:make-waitqueue)
  )

(defun condition-wait (waitqueue mutex)
  "Registers on the waitqueue, releases mutex, and waits for a notification
on the waitqueue."
  (declare (ignorable waitqueue mutex))
  #+sb-thread (sb-thread:condition-wait waitqueue mutex)
  )

(defun condition-notify (waitqueue)
  "Notifies on the waitqueue."
  (declare (ignorable waitqueue))
  #+sb-thread (sb-thread:condition-notify waitqueue)
  )

(defclass waitqueue-mixin (mutex-mixin)
  ((waitqueue :reader waitqueue :initform (make-waitqueue)))
  (:documentation "Waitqueue mixin."))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; MRSW (multiple-read-single-write) access (untested)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass mrsw-mixin (waitqueue-mixin)
  ((readcount :reader readcount :initform 0 :documentation
	      "Number of read accesses.")
   (writecount :reader writecount :initform 0 :documentation
	       "Number of (desired) write accesses, should be 0 or 1."))
  (:documentation "Waitqueue mixin."))

(defmacro with-read-access ((obj) &body body)
  "Execute @arg{body} while registered for reading @arg{obj}."
  `(with-slots (mutex waitqueue readcount writecount) ,obj
    (with-mutex (mutex)
      (loop until (zerop writecount) do
	    (condition-wait waitqueue mutex))
      (incf readcount)
      ,@body
      (decf readcount))))

(defmacro with-readwrite-access ((obj) &body body)
  "Execute @arg{body} while registered for writing @arg{obj}."
  `(with-slots (mutex waitqueue readcount writecount) ,obj
    (with-mutex (mutex)
      (loop until (zerop writecount) do
	    (condition-wait waitqueue mutex))
      (incf writecount)
      (loop until (zerop readcount) do
	    (condition-wait waitqueue mutex))
      (incf readcount))
    ,@body
    (with-mutex (mutex)
      (decf readcount)
      (decf writecount))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; thread-safe queue
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass parqueue (queue waitqueue-mixin)
  ((finished-p :reader finished-p :initform nil :documentation
	       "Indicator for end-of-input."))
  (:documentation "A thread-safe queue waiting for input."))

(defmethod enqueue :around (obj (pq parqueue))
  (declare (ignorable obj))
  (with-mutual-exclusion (pq)
    (call-next-method)
    (condition-notify (waitqueue pq))))

(defmethod dequeue :around ((pq parqueue))
  (loop do
	(when (finished-p pq)
	  (with-mutual-exclusion (pq)
	    (return (call-next-method))))
	(with-mutual-exclusion (pq)
	  (whereas ((obj (call-next-method)))
	    (return obj))
	  (condition-wait (waitqueue pq) (mutex pq)))
	(thread-yield)))

(defmethod finish-queue ((parqueue parqueue))
  (with-mutual-exclusion (parqueue)
    (setf (slot-value parqueue 'finished-p)  t)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Working-group
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass working-group (waitqueue-mixin)
  ((threads :reader threads :initform ())
   (queue :reader queue :initform (make-instance 'parqueue)))
  (:documentation "A class representing a group of threads working on a
queue."))

(defun add-worker (group thread)
  "Adds @arg{thread} as worker to @arg{group}."
  (with-mutual-exclusion (group)
    (dbg :mp "Adding worker ~A" thread)
    (with-slots (threads waitqueue) group
      (push thread threads)
      (condition-notify waitqueue))))

(defun remove-worker (group thread)
  "Removes @arg{thread} from @arg{group}."
  (with-mutual-exclusion (group)
    (with-slots (threads waitqueue queue) group
      (assert (finished-p queue))
      (assert (member thread threads))
      (dbg :mp "Removing worker ~A" thread)
      (setf threads (delete thread threads))
      (condition-notify waitqueue))))

(defun send-work (group work)
  "Send work to the group.  @arg{work} is an argument list for the worker
functions."
  (with-mutual-exclusion (group)
    (enqueue (queue group) work)))

(defun wait-until-finished (group)
  (loop while (threads group) do
	(with-slots (waitqueue mutex) group
	  (with-mutex (mutex)
	    (condition-wait waitqueue mutex))
	  (thread-yield))))

(defparameter *thread-local-memoization-table* nil
  "Dynamic variable specialized for each worker thread separately.  Should
not be used globally.")

(defmacro with-workers ((func) &body body)
  "This macro distributes work generated in body with calling the locally
bound function @function{work-on} on some arguments to several working
threads which call @arg{func} on those arguments."
  (if *threads*
      (with-gensyms (group thunk)
	`(let ((,group (make-instance 'working-group)))
	  (flet ((work-on (&rest args) (enqueue args (queue ,group)))
		 (,thunk ()
		   (let ((*thread-local-memoization-table*
			  (make-hash-table :test 'equalp)))
		     (loop for args = (dequeue (queue ,group)) while args do
			   (dbg :mp "Handling ~A" args)
			   (apply ,func args)
			   (thread-yield))
		     (remove-worker ,group (current-thread)))))
	    (loop repeat *threads* do
		  (add-worker ,group (make-thread #',thunk :name "worker")))
	    ,@body
	    (finish-queue (queue ,group))
	    (wait-until-finished ,group))))
      `(flet ((work-on (&rest args) (apply ,func args)))
	,@body)
      ))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Thread-global and thread-local cacheing of results
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; in general.lisp

;;; Testing
(defun femlisp-multiprocessing-tests ()
  (let ((q (make-instance 'parqueue)))
    (enqueue 1 q)
    (dequeue q))

  (let ((k 100000000)
	(mutex (make-mutex))
	(sum 0))
    (declare (ignorable mutex))
    (time (print (loop for n to k summing 1)))
    (time
     (let ((*threads* 10))
       (with-workers
	   ((lambda (k)
	      (declare (type fixnum k))
	      (let ((s (loop for n to k summing 1)))
		(with-mutex (mutex)
		  (incf sum s)))))
	 (loop for i below *threads* do (work-on k)))
       (format t "Sum: ~A" sum))))

  (let ((x (make-instance 'mrsw-mixin)))
    (with-read-access (x)
      (print 'hi))
    (with-readwrite-access (x)
      (print 'ho)))
  )

;;; (femlisp-multiprocessing-tests)