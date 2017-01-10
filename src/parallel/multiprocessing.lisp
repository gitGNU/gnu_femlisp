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

(in-package :fl.parallel)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Multithreaded printing
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defvar *print-lock* (make-recursive-lock "PRINT-LOCK"))

(defmacro with-atomic-output (&body body)
  "If output of a process is desired to be atomic
wrap it in @arg{with-atomic-output}."
  `(with-recursive-lock-held (*print-lock*)
     ,@body))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Multithreaded debugging
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defvar *debug-lock* (make-recursive-lock "DEBUG-LOCK"))

(defmethod dbg :around (id format-string &rest args)
  "Serializing of debug output."
  (declare (ignore id format-string args))
  (with-recursive-lock-held (*debug-lock*)
    (call-next-method)))

(defun mp-dbg (id format &rest args)
  (funcall #'dbg id "~A:~%~?" (current-thread) format args))

;;; Some CLOS mixins for mutex and waitqueues 

(defclass mutex-mixin ()
  ((mutex :reader mutex :initform (bordeaux-threads:make-recursive-lock) :documentation
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
            (bordeaux-threads:with-recursive-lock-held (it) (,do))
            (error "No object with mutex")))))

(defclass waitqueue-mixin (mutex-mixin)
  ((waitqueue :reader waitqueue :initform (make-condition-variable)))
  (:documentation "Waitqueue mixin."))

(defgeneric wait (waitqueue-object &key while until finish perform)
  (:documentation "Waits on @arg{waitqueue-object} while @arg{while} is
satisfied or until @arg{until} is satisfied.  After a successful waiting,
the function given in @arg{perform} is called.")
  (:method (obj &key &allow-other-keys)
    "The default method for objects does not wait at all."
    nil)
  (:method ((wq waitqueue-mixin) &key while until perform finish)
    (with-mutual-exclusion (wq)
      (loop while (or (aand while (funcall it wq))
		      (aand until (not (funcall it wq))))
	    do
               (when (aand finish (funcall it))
                 (mp-dbg :mp "Finishing wait...")
                 (notify wq)
                 (return nil))
               (mp-dbg :mp "Waiting...")
               (condition-wait (waitqueue wq) (mutex wq))
            finally
               (mp-dbg :mp "Performing after waiting...")
               (return (awhen perform (funcall it wq)))))))

(defgeneric notify (waitqueue-object)
  (:documentation "Notifies @arg{waitqueue-object}.")
  (:method (obj) "The default method does nothing.")
  (:method ((wq waitqueue-mixin))
    ;;(mp-dbg :mp "notifying ~A" wq)
    (condition-notify (waitqueue wq))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; locked-region
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass locked-region-mixin (waitqueue-mixin)
  ((locked-region :reader locked-region
                  :initform (make-hash-table :test 'equal))))

(defun lock-region (object keys)
  "Should not be used externally.  Use WITH-REGION instead."
  (with-mutual-exclusion (object)
    (wait object
          :until
          (lambda (object)
            (let ((table (locked-region object)))
              (notany (lambda (key)
                        (aand (gethash key table)
                              (not (eq it (current-thread)))))
                      keys))))
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
    (lock-region object keys)
    (funcall perform)
    (unlock-region object keys)
    (notify object)))

(defmacro with-region ((object keys) &body body)
  `(perform-with-locked-region ,object ,keys
		   (lambda () ,@body)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; thread-safe queue
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass parqueue (fl.utilities::queue waitqueue-mixin)
  ((finished-p :accessor finished-p :initform nil :documentation
	       "Indicator for no-more-input-allowed.")
   (maximal-size :initform nil :initarg :maximal-size)
   (current-size :initform 0))
  (:documentation "A thread-safe queue waiting for input."))

(defgeneric finish (object)
  (:method ((pq parqueue))
    (with-mutual-exclusion (pq)
      (mp-dbg :mp "finishing queue ~A" pq)
      (setf (finished-p pq) t)
      (notify pq))))

(defmethod enqueue :around (obj (pq parqueue))
  (with-mutual-exclusion (pq)
    (with-slots (finished-p maximal-size current-size) pq
      (assert (not finished-p))
      (mp-dbg :mp "enqueueing ~A in ~A" obj pq)
      (loop while (and maximal-size (>= current-size maximal-size)) do
        (mp-dbg :mp "queue ~A full - waiting" pq)
        (condition-wait (waitqueue pq) (mutex pq)))
      (call-next-method)
      (incf current-size)
      (notify pq) ; notify of change!
      nil
      )))

(defmethod dequeue :around ((pq parqueue))
  (with-mutual-exclusion (pq)
    (with-slots (finished-p current-size) pq
      (loop
        (mp-dbg :mp "trying to dequeue ~A (finished: ~A, empty: ~A, top: ~A)"
                pq finished-p (emptyp pq) (car (fl.utilities::head pq)))
        (cond ((or finished-p (not (emptyp pq)))
               (multiple-value-bind (value emptyp)
                   (call-next-method)
                 (unless emptyp (decf current-size))
                 (notify pq)           ; notify of change
                 (mp-dbg :mp "dequeued value=~A, emptyp=~A" value emptyp)
                 (return-from dequeue (values value emptyp))))
              (t (mp-dbg :mp "waiting on parqueue ~A" pq)
                 (condition-wait (waitqueue pq) (mutex pq))))
         ))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Multithreaded pool
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass parpool (mutex-mixin)
  ((table :type hash-table
          :documentation "A table mapping keys to a list of parpool entries")
   (ppes :type hash-table
         :documentation "A table mapping an object to its ppe")
   (test :initform 'eql :initarg :test)
   #+(or)(generator :initform (required-argument) :initarg :generator))
  (:documentation "A parallel pool implemented as a shared hash-table.
The slot pool contains a hash-table, the test slot contains the hash-table test, the generator contains a function which is called if the pool is empty, i.e. all objects are used and generates a new object for the given argument."))

(defstruct (parpool-entry (:conc-name ppe-))
  object refcount key)

(defmethod initialize-instance :after ((pool parpool) &key &allow-other-keys)
  (with-slots (table test ppes) pool
    (setf table (make-hash-table :test test))
    (setf ppes (make-hash-table :test 'eq))))

(defmethod describe-object :after ((ppool parpool) stream)
  (format stream "~A contains the following:~%" ppool)
  (maphash (lambda (key value)
             (format stream "~A ->~%~{ ~A~}~%" key value))
           (slot-value ppool 'table))
  (format stream "Furthermore, the following objects are registered:~%")
  (maphash (lambda (key value)
             (format stream "~A -> ~A~%" key value))
           (slot-value ppool 'ppes)))

(defun register-in-pool (pool key object)
  "Note that the object is only registered, but not put in the pool!"
  (with-mutual-exclusion (pool)
    (with-slots (table ppes) pool
      (let ((ppe (make-parpool-entry :object object :key key)))
        (let (*print-array*)
          (dbg :parpool "Registering new ppe ~A for object ~A under key=~A"
               ppe object key))
        (assert (not (gethash object ppes)))
        (setf (gethash object ppes) ppe)))))

(defun get-from-pool (pool key)
  (with-mutual-exclusion (pool)
    (with-slots (table) pool
      (whereas ((ppe (pop (gethash key table nil))))
        (lret ((object (ppe-object ppe)))
          (let ((*print-array*))
            (dbg :parpool "Getting from pool ppe ~A for object ~A under key=~A"
                 ppe object key)))))))

(defun put-back-in-pool (pool object)
  (with-mutual-exclusion (pool)
    (with-slots (ppes table) pool
      (let ((ppe (gethash object ppes)))
        (assert ppe)
        (let ((key (ppe-key ppe)))
          (assert (not (member ppe (gethash key table) :test 'eq)))
          (let (*print-array*)
            (dbg :parpool "Putting back in pool ppe ~A for object ~A under key=~A"
               ppe object key))
          (push ppe (gethash key table)))))))

(defun set-refcount (pool object count)
  (with-mutual-exclusion (pool)
    (with-slots (ppes) pool
      (let ((ppe (gethash object ppes)))
        (setf (ppe-refcount ppe) count)))))

(defun decrease-refcount (pool object &optional (decrement 1))
  (with-mutual-exclusion (pool)
    (with-slots (ppes) pool
      (let ((ppe (gethash object ppes)))
        (when (zerop (decf (ppe-refcount ppe) decrement))
          (setf (ppe-refcount ppe) nil)
          (put-back-in-pool pool object))))))

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
       (mp-dbg :mp "starting...")
       ;; (let ((*thread-local-memoization-table* ()))
       ;;   (declare (special *thread-local-memoization-table*))
       (if (null from)
           (funcall job)  ; distributor
             (loop
                (multiple-value-bind (args emptyp)
                    (dequeue from)
                  (when emptyp
                    (return))
                  (mp-dbg :mp "working on ~A" args)
                  (apply job args)
                  (mp-dbg :mp "finished working on ~A" args)
                  (thread-yield)))
                )
       (mp-dbg :mp "...done")
       (awhen to (finish it))
       #+(or)(remove-worker group (current-thread))))
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

(defun add-worker (group)
  (with-mutual-exclusion (group)
    (push (make-thread
	   (lambda ()
	     (mp-dbg :mp "starting...")
	     ;; (let ((*thread-local-memoization-table* ()))
             (loop
               (multiple-value-bind (args emptyp)
                   (dequeue (tasks group))
                 (when emptyp (return))
                 (mp-dbg :mp "working on ~A" args)
                 (apply (work group) (mklist args))
                 (mp-dbg :mp "finished working on ~A." args)
                 (thread-yield)))
             (mp-dbg :mp "...done")
             (remove-worker group (current-thread)))
	   :name "femlisp-worker")
	  (threads group))))

(defun remove-worker (group thread)
  "Removes @arg{thread} from @arg{group}."
  (assert (finished-p (tasks group)))
  (with-mutual-exclusion (group)
    (assert (member thread (threads group)))
    (mp-dbg :mp " --- removes itself.")
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
     (member (thread-name thread) '("femlisp-worker" "pipeline-worker" "stencil-worker" "lparallel")
             :test #'string=))
   (all-threads)))

(defun terminate-workers ()
  (dolist (thread (femlisp-workers))
    (destroy-thread thread)))

(defmacro with-femlisp-workers ((work) &body body)
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
  "Old code.  The optimal thread number depends on the kind of application and cannot really be determined in this simple way.  Probably better is to use the introspection in @file{parallel.lisp}."
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
  (declare (ftype (function ((complex double-float)) fixnum))
           (optimize speed))
  (?1
   ;; Version 1
   (loop for k of-type fixnum from 1 below 1000
         and z of-type (complex double-float) = c
               then (+ (* z z) c)
         until (> (abs z) 100.0d0)
         finally (return k))
   ;; Version 2
   (labels ((mi1 (z n)
              (declare (type (complex double-float) z)
                       (type fixnum n))
              (if (or (>= n 1000) (> (abs z) 100.0d0))
                  n
                  (mi1 (+ (* z z) c) (1+ n)))))
     (mi1 c 1))
   ;; Version 3
   (named-let mi1 ((z c) (n 1))
     (declare (type (complex double-float) z)
              (type fixnum n))
     (if (or (>= n 1000) (> (abs z) 100.0d0))
         n
         (mi1 (+ (* z z) c) (1+ n))))
   ))

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
          and x of-type double-float = x1 then (+ x dx) do
            (loop for ky upto ny
                  and y of-type double-float = y1 then (+ y dy) do
                    (incf sum
                          (setf (aref result kx ky)
                                (mandelbrot-iteration (complex x y))))))
    (if result-p result sum)))

(defun histogram (array)
  "Calculates a Histogram of Mandelbrot iteration numbers."
  (lret* ((max (loop for k below (array-total-size array)
                     maximize (row-major-aref array k)))
          (result (make-array (1+ max) :initial-element 0)))
    (loop for k below (array-total-size array)
          do (incf (aref result (row-major-aref array k))))))

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
  (let ((lock (make-recursive-lock))
	(result ())
	(pq (make-instance 'parqueue)))
    (flet ((worker ()
	     (loop for obj = (dequeue pq) while obj
		   do (with-recursive-lock-held (lock)
			(push obj result)))
	     (mp-dbg :mp "Done.~%")))
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
   (flet ((test (i)
            (+ i (* i i))))
     (let ((test2 (lambda (&rest args)
                    (apply #'test args))))
       (with-femlisp-workers (test2)
         (let ((n 100))
           (loop for i below n do
             (work-on i)))))))
  (time
   (let ((n 10)
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
	(lock (make-recursive-lock))
	(sum 0))
    (declare (ignorable lock))
    (time (print (loop repeat k summing 1)))
    (time
     (let ((*number-of-threads* 4))
       (with-workers
	   ((lambda (k)
	      (declare (type fixnum k))
	      (let ((s (loop repeat k summing 1)))
		(with-recursive-lock-held (lock)
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
  (mandelbrot-box -2.0d0 1.0d0 -1.0d0 1.0d0 2 2)
  (time (mandelbrot-box -2.0d0 1.0d0 -1.0d0 1.0d0 1536 1024))
  (histogram (mandelbrot-box -2.0d0 1.0d0 -1.0d0 1.0d0 256 256 t))
  (speedup-test (_ (mandelbrot-box 0.0 10.0 0.0 10.0 100 100)))
  (speedup-test (_ (simple-consing 500000)))

  (sleep 0.5)
  (assert (null (femlisp-workers)) () "workers remained - F")
  ;; (all-threads)
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
  
  ;; (all-threads)
  ;; (terminate-workers)
  (dolist (thread (all-threads))
    (when (string= (thread-name thread) "worker-1")
      (destroy-thread thread)))
  (dolist (thread (all-threads))
    (when (string= (thread-name thread) "lparallel")
      (destroy-thread thread)))

  (dbg-off :mp)

  (let ((pool (make-instance 'parpool)))
    (let ((o1 (vector 1))
          (o2 (vector 1 1)))
      (register-in-pool pool 1 o1)
      (register-in-pool pool 1 o2)
      (let ((o (get-from-pool pool 1)))
        (format t "~&Got ~A~%" o)
        (set-refcount pool o 4)
        (describe pool)
        (decrease-refcount pool o)
        (decrease-refcount pool o)
        (decrease-refcount pool o)
        (decrease-refcount pool o)
        (describe pool)
        )))

  )

;;; (femlisp-multiprocessing-tests)
