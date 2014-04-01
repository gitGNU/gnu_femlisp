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
   "MAKE-THREAD" "CURRENT-THREAD" "THREAD-NAME"
   "LIST-ALL-THREADS" "TERMINATE-THREAD"
   
   "MAKE-MUTEX" "WITH-MUTEX"
   "MUTEX-MIXIN" "WITH-LOCK"
   "WITH-MUTUAL-EXCLUSION"

   "WAITQUEUE-MIXIN" "MAKE-WAITQUEUE" "CONDITION-WAIT" "CONDITION-NOTIFY"
   
   "MP-DBG"
   "LOCKED-REGION-MIXIN"
   "WITH-REGION"
   "MPQUEUE" "PARQUEUE" "WAIT" "NOTIFY"
   "*THREAD-LOCAL-MEMOIZATION-TABLE*"
   "*NUMBER-OF-THREADS*"
   "WITH-WORKERS" "WORK-ON" "FEMLISP-WORKERS"

   "WITH-PARALLEL-NETWORK" "PARCELL" "NAMED-PARCELL" "VALUE"
   )
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
  #+lispworks (mp:process-run-function name nil func)
  #+sb-thread (sb-thread:make-thread func :name name)
  #+scl (thread:thread-create func :name name)
  #-(or allegro cmu scl lispworks sb-thread) nil
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
  #+lispworks mp:*current-process*
  #+sb-thread sb-thread:*current-thread*
  #+scl thread:*thread*
  #-(or allegro cmu lispworks sb-thread scl)
  (load-time-value '(:main-thread))
  )

(defun thread-name (thread)
  "Returns the name of the current thread."
  #+allegro (mp:process-name thread)
  #+(or cmu lispworks) (mp:process-name thread)
  #+sb-thread (sb-thread:thread-name thread)
  #+scl (thread:thread-name thread)
  #-(or allegro cmu lispworks sb-thread scl) nil
  )

(defun list-all-threads ()
  "Returns a list of all threads"
  #+allegro mp:*all-processes*
  #+cmu (mp:all-processes)
  #+lispworks (mp:list-all-processes)
  #+sb-thread (sb-thread:list-all-threads)
  #+scl
  (lret ((result ()))
    (thread:map-over-threads
     (lambda (thread) (push thread result))))
  #-(or allegro cmu lispworks sb-thread scl)
  (list :main-thread)
  )

(defun terminate-thread (thread)
  #+allegro (mp:process-kill thread)
  #+cmu (mp:destroy-process thread)
  #+lispworks (mp:process-kill thread)
  #+sb-thread (sb-thread:terminate-thread thread)
  #+scl (thread:destroy-thread thread)
  #-(or allegro cmu lispworks sb-thread scl)
  ;; no threading: thus only the main-thread which must not be killed
  (assert (null thread))
  )
  
;;; mutices

(defun make-mutex (&key (name "mutex"))
  (declare (ignorable name))
  #+allegro (mp:make-process-lock :name name)
  #+cmu (mp:make-lock name)
  #+lispworks (mp:make-lock :name name)
  #+sb-thread (sb-thread:make-mutex :name name)
  #+scl (thread:make-lock name :type :recursive)
  #-(or allegro cmu lispworks sb-thread scl) nil
  )

(defmacro with-mutex ((mutex) &body body)
  (declare (ignorable mutex))
  #+allegro `(mp:with-process-lock (,mutex) ,@body)
  #+cmu `(mp:with-lock-held (,mutex) ,@body)
  #+lispworks `(mp:with-lock (,mutex) ,@body)
  #+sb-thread `(sb-thread:with-recursive-lock (,mutex) ,@body)
  #+scl `(thread:with-lock-held (,mutex) ,@body)
  #-(or allegro cmu lispworks sb-thread scl) `(locally ,@body))

;;; waitqueues (following the SBCL interface)

(defun make-waitqueue ()
  "Generates a waitqueue."
  #+allegro (mp:make-gate nil)
  #+sb-thread (sb-thread:make-waitqueue)
  #+scl (thread:make-cond-var)
  #-(or allegro sb-thread scl) nil
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
  #+scl (thread:cond-var-wait waitqueue mutex)
  #-(or allegro sb-thread scl) (thread-yield)
  )

(defun condition-notify (waitqueue &optional (n 1))
  "Notifies on the waitqueue."
  (declare (ignorable waitqueue n))
  #+allegro (mp:open-gate waitqueue)
  #+sb-thread
  (if (eq n t)
      (sb-thread:condition-broadcast waitqueue)
      (sb-thread:condition-notify waitqueue n))
  #+scl
  (if (eq n t)
      (thread:cond-var-signal waitqueue)
      (thread:cond-var-broadcast waitqueue))
  #-(or allegro sb-thread scl) (thread-yield)
  )

;;; nicer OO interface

(defclass mutex-mixin ()
  ((mutex :reader mutex :initform (make-mutex) :documentation
	  "A mutex for excluding access of other threads."))
  (:documentation "A mixin which adds a mutex to every instance of the
class."))

(defmethod mutex (obj)
  "Ordinary objects are not mutex-protected."
  nil)

(defmacro with-mutual-exclusion ((obj) &body body)
  "Execute @arg{body} on the waitqueue @arg{obj} without other threads
interfering."
  (with-gensyms (do)
    `(flet ((,do () ,@body))
       (aif (mutex ,obj)
            (with-mutex (it) (,do))
            (,do)))))

(defclass waitqueue-mixin (mutex-mixin)
  ((waitqueue :reader waitqueue :initform (make-waitqueue)))
  (:documentation "Waitqueue mixin."))

(defgeneric wait (waitqueue-object &key while until perform)
  (:documentation "Waits on @arg{waitqueue-object} while @arg{while} is
satisfied or until @arg{until} is satisfied.  After this, if it is given,
the function @arg{perform} is called.")
  (:method (obj &key &allow-other-keys)
    "The default method for objects does not wait at all.")
  (:method ((wq waitqueue-mixin) &key while until perform)
    (with-mutual-exclusion (wq)
      (loop
         do (condition-wait (waitqueue wq) (mutex wq))
         when (or (not (and while until))
                  (aand while (not (funcall it wq)))
                  (aand until (funcall it wq)))
         do (loop-finish))
      (awhen perform (funcall it wq)))))

(defgeneric notify (waitqueue-object &optional n)
  (:documentation "Notifies @arg{waitqueue-object}.")
  (:method (obj &optional n)
    "The default method does nothing."
    (declare (ignore n)))
  (:method ((wq waitqueue-mixin) &optional (n t))
    ;;(mp-dbg "notifying ~A" wq)
    (condition-notify (waitqueue wq) n)))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; locked-region
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass locked-region-mixin (waitqueue-mixin)
  ((locked-region :reader locked-region)))

(defmethod initialize-instance :after ((self locked-region-mixin)
                                       &key (lock-test 'eql) &allow-other-keys)
  (setf (slot-value self 'locked-region)
        (make-hash-table :test lock-test)))

(defun lock-region (object keys)
  "Should not be used externally.  Use WITH-REGION instead."
  (with-mutual-exclusion (object)
    (let ((table (locked-region object)))
      (loop+ ((key keys)) do
	 (awhen (gethash key table)
	   (assert (eq it (current-thread))))
	 (setf (gethash key table)
	       (current-thread))))))

(defun unlock-region (object keys)
  "Should not be used externally.  Use WITH-REGION instead."
  (with-mutual-exclusion (object)
    (let ((table (locked-region object)))
      (loop+ ((key keys)) do
	 (remhash key table)))))

(defgeneric perform-with-locked-region (object keys perform)
  (:documentation
   "Perform @arg{perform} while @arg{keys} are locked in @arg{object}.")
  (:method (object keys perform)
    "No locking if @arg{object} is not a @class{locked-region-mixin}."
    (funcall perform))
  (:method ((object locked-region-mixin) keys perform)
    (with-mutual-exclusion (object)
      (wait object
            :until
            (lambda (object)
              (let ((table (locked-region object)))
                (lret ((result (notany (lambda (key)
                                         (aand (gethash key table)
                                               (not (eq it (current-thread)))))
                                       keys)))))))
      (lock-region object keys))
    (funcall perform)
    (unlock-region object keys)
    (notify object)))

(defmacro with-region ((object keys) &body body)
  `(perform-with-locked-region ,object ,keys
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

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; thread-safe queue
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass parqueue (queue waitqueue-mixin)
  ((finished-p :accessor finished-p :initform nil :documentation
	       "Indicator for no-more-input-allowed.")
   (maximal-size :initform nil :initarg :maximal-size)
   (current-size :initform 0))
  (:documentation "A thread-safe queue waiting for input."))

(defgeneric finish (object)
  (:method ((pq parqueue))
    (with-mutual-exclusion (pq)
      (mp-dbg "finishing queue ~A" pq)
      (setf (finished-p pq) t)
      (notify pq))))

(defmethod enqueue :around (obj (pq parqueue))
  (with-slots (finished-p maximal-size current-size) pq
    (assert (not finished-p))
    (mp-dbg "enqueueing ~A in ~A" obj pq)
    (with-mutual-exclusion (pq)
      (loop while (and maximal-size (>= current-size maximal-size)) do
           (mp-dbg "queue ~A full - waiting" pq)
           (condition-wait (waitqueue pq) (mutex pq)))
      (call-next-method)
      (incf current-size)
      (notify pq) ; notify of change!
      nil
      )))

(defmethod dequeue :around ((pq parqueue))
  (with-slots (finished-p current-size) pq
    (with-mutual-exclusion (pq)
      (loop
         (mp-dbg "trying to dequeue ~A (finished: ~A, empty: ~A, top: ~A)"
                 pq finished-p (emptyp pq) (car (fl.utilities::head pq)))
         (cond ((or finished-p (not (emptyp pq)))
                (multiple-value-bind (value emptyp)
                    (call-next-method)
                  (unless emptyp (decf current-size))
                  (notify pq)           ; notify of change
                  (mp-dbg "dequeued value=~A, emptyp=~A" value emptyp)
                  (return-from dequeue (values value emptyp))))
               (t (mp-dbg "waiting on parqueue ~A" pq)
                  (condition-wait (waitqueue pq) (mutex pq))))
         ))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; pipeline
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defvar *pipeline-length* 10)
(defvar *distribute-workload-p* nil)

(defun execute-in-pipeline (distribute &rest jobs)
  "Every function in @arg{jobs} produces arguments for the next one with the
function @function{call-next}.  Every job is run in a separate thread."
  (let ((joins (loop repeat (length jobs) collecting
                    (make-instance 'parqueue
                                   :maximal-size *pipeline-length*))))
    (let ((first-jobs (butlast jobs))
          (collector (car (last jobs)))
          (first-joins (butlast joins))
          (last-join (car (last joins))))
      (mapc #'add-pipeline-worker
            (cons distribute first-jobs)
            (cons nil first-joins)
            joins)
      ;; collector
      (loop (multiple-value-bind (args emptyp)
                (dequeue last-join)
              (when emptyp (return))
              (apply collector args)
              (thread-yield))))))

(defvar *input-queue* nil
  "The input queue of a worker job in a pipeline.")

(defvar *output-queue* nil
  "The output queue of a worker job in a pipeline.")

(defun hand-over (&rest args)
  "Hands over the arguments to the next worker job."
  (enqueue args *output-queue*))

(defun add-pipeline-worker (job from to)
  (make-thread
   (lambda ()
     (let ((*input-queue* from)
           (*output-queue* to))
       (mp-dbg "starting...")
       (let ((*thread-local-memoization-table* ()))
         (declare (special *thread-local-memoization-table*))
         (if (null from)
             (funcall job)  ; distributor
             (loop
                (multiple-value-bind (args emptyp)
                    (dequeue from)
                  (when emptyp
                    (return))
                  (mp-dbg "working on ~A" args)
                  (apply job args)
                  (mp-dbg "finished working on ~A" args)
                  (thread-yield)))
                )
         (mp-dbg "...done")
         (awhen to (finish it))
         #+(or)(remove-worker group (current-thread)))))
   :name "pipeline-worker"))

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

(defun add-worker (group)
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
		  (apply (work group) (mklist args))
		  (mp-dbg "finished working on ~A." args)
		  (thread-yield)))
	       (mp-dbg "...done")
	       (remove-worker group (current-thread))))
	   :name "femlisp-worker")
	  (threads group))))

(defun remove-worker (group thread)
  "Removes @arg{thread} from @arg{group}."
  (assert (finished-p (tasks group)))
  (with-mutual-exclusion (group)
    (assert (member thread (threads group)))
    (mp-dbg " --- removes itself.")
    (setf (threads group)
	  (delete thread (threads group)))
    (notify group)))

(defun send-task (group task)
  "Send a task to the work-group."
  (enqueue task (tasks group)))

(defvar *number-of-threads* nil  ; (optimal-thread-number)
  "The number of threads in which Femlisp tries to split the work for some
computationally intensive tasks.  If NIL, no threading is used.")

;; (setq *number-of-threads* nil)
;; (setq *number-of-threads* 2)

(defun execute-in-parallel (work &key arguments distribute
			    (number-of-threads *number-of-threads*))
  "Executes work in parallel distributed on @arg{number-of-threads}
threads.  The arguments for work are taken from the argument-lists in
@arg{arguments}.  Alternatively, @arg{distribute} may contain a function
which generates the argument-lists by calling the function
@function{work-on}."
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

(defun femlisp-workers ()
  (remove-if-not
   (lambda (thread)
     (member (thread-name thread) '("femlisp-worker" "pipeline-worker" "stencil-worker")
             :test #'string=))
   (list-all-threads)))

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
	(flet ((work-on (&rest args)
		 (declare (type function ,send-work))
		 (apply ,send-work args)))
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

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; For testing
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; the following does almost no consing (with SBCL, at least) and can
;;; therefore be efficiently parallelized

(defun mandelbrot-iteration (c)
  (declare (ftype (function (complex double-float) fixnum))
           (optimize speed))
  (loop for k of-type fixnum from 1 below 1000
       for z of-type (complex double-float) = c then
       (+ (* z z) c)
       until (> (abs z) 100.0)
       finally (return k)))

(defun mandelbrot-box (x1 x2 y1 y2 nx ny &optional result-p)
  (declare (type double-float x1 x2 y1 y2)
           (type fixnum nx ny))
  (let ((dx (/ (- x2 x1) nx))
        (dy (/ (- y2 y1) ny))
        (result (make-array (list (1+ nx) (1+ ny)) :element-type 'fixnum))
        (sum 0))
    (declare (type double-float dx dy)
             (type fixnum sum))
    (loop for kx upto nx
       for x of-type double-float = x1 then (+ x1 dx) do
       (loop for ky upto ny
          for y of-type double-float = y1 then (+ y1 dy) do
            (incf sum
                  (setf (aref result kx ky)
                        (mandelbrot-iteration (complex x y))))))
    (if result-p result sum)))

(defun simple-consing (n)
  (loop repeat n do
       (loop repeat 100 collect nil)))

(defun speedup-test (func)
  (loop for i from 1 upto 8 do
       (format t "~R thread~:P~%" i)
       (let ((*number-of-threads* i))
         (time (with-workers (func)
                 (loop repeat i do (work-on)))))))

;;; Testing
(defun femlisp-multiprocessing-tests ()

  (dbg-on :mp)
  #-scl  ; does not work on SCL
  (let ((mutex (make-mutex))
	(result ())
	(pq (make-instance 'parqueue)))
    (flet ((worker ()
	     (loop for obj = (dequeue pq) while obj
		   do (with-mutex (mutex)
			(push obj result)))
	     (mp-dbg t "~A done.~%" (current-thread))))
      (make-thread #'worker :name "femlisp-worker")
      (make-thread #'worker :name "femlisp-worker")
      (make-thread #'worker :name "femlisp-worker")
      )
    (loop for k from 1 upto 10 do
	  (enqueue k pq))
    (finish pq)
    (sleep 0.5)
    result)
  (dbg-off :mp)

  (assert (null (femlisp-workers)) () "workers remained - A")
  ;; (terminate-workers)
  
  (time
   (let ((n 100000)
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

  (speedup-test (let ((n (expt 2 27))) (_ (loop repeat n sum 1))))
  (speedup-test (_ (mandelbrot-box 0.0 10.0 0.0 10.0 100 100)))
  (speedup-test (_ (simple-consing 500000)))

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

  (time
   (lret ((delay 0.5)
          (result ()))
     (execute-in-pipeline
      (_ (loop for i below 10 do
              (hand-over i)))
      (_ (sleep delay)
         (hand-over (expt _ 2)))
      (_ (sleep delay)
         (hand-over (isqrt _)))
      (_ (push _ result)))))
  
  ;; (list-all-threads)
  ;; (terminate-workers)
  (dolist (thread (list-all-threads))
    (when (string= (thread-name thread) "worker-1")
      (terminate-thread thread)))

  (dbg-off :mp)
  )


;;; (femlisp-multiprocessing-tests)

