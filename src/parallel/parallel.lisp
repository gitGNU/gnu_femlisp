;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Parallel utilities
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(in-package :fl.parallel)

;;;; CPU info

(defun get-cpuinfo ()
  (with-open-file (stream "/proc/cpuinfo")
    (loop for line = (read-line stream nil) while line
       collect (cl-ppcre:split "\\t+: " line))))

;;; Cachesize

(defun get-cachesize ()
  (let* ((size (second (assoc "cache size" (get-cpuinfo) :test #'string=)))
         (sep-pos (position #\Space size))
         (number (parse-integer (subseq size 0 sep-pos)))
         (unit (subseq size (1+ sep-pos))))
    (* number
       (stringcase unit
         ("KB" 1024)
         (t (error "Unknown unit: ~A" unit))))))

#+(or)
(defun calculate-effective-cachesize ()
  "Calculates effective cachesize in bytes."
  (let ((l 30))
    (flet ((test (k)
             (let* ((n (expt 2 k))
                    (count (expt 2 (- l k)))
                    (x (make-double-vec n))
                    (y (make-double-vec n)))
               (declare (optimize speed (safety 0)))
               (tic)
               (loop repeat count
                  do
                    (replace x y))
                (tic))))
      (loop for k from 10 below 30
         for previous = nil then next
         and next = (test k)
         until (and previous (> next (* 1.05 previous)))
         do (format t "~2D ~A~%" k next)
         finally (return (expt 2 (+ k 1)))))))


(defvar *cachesize* (get-cachesize))

(unless (<= (expt 2 19) *cachesize* 32e6)
  (warn "Possible problem: cachesize not in standard range"))

;;; Processors

(defstruct (procinfo (:conc-name pi-) (:type list))
  processor
  physical-id
  core-id
  apicid)
           
(defun get-processors ()
  (let (processors current-processor)
    (loop for entry in (get-cpuinfo)
         when (= (length entry) 2) do
         (destructuring-bind (id value) entry
           (stringcase id
             ("processor"
              (when current-processor
                (push current-processor processors))
              (setf current-processor
                    (make-procinfo :processor (parse-integer value))))
             ("physical id"
              (setf (pi-physical-id current-processor)
                    (parse-integer value)))
             ("core id"
              (setf (pi-core-id current-processor)
                    (parse-integer value)))
             ("apicid"
              (setf (pi-apicid current-processor)
                    (parse-integer value)))))
         finally (when current-processor
                   (push current-processor processors)))
    (reverse processors)))

(defun allowed-processors ()
  #+sb-cpu-affinity
  (sb-cpu-affinity:with-cpu-affinity-mask (mask)
    (loop for proc in (get-processors)
          when (sb-cpu-affinity:cpu-affinity-p (pi-processor proc) mask)
            collect proc))
  #-sb-cpu-affinity
  (get-processors))

(defun get-processors-without-hyperthreading ()
  (let ((procs ()))
    (loop for proc in (allowed-processors)
       unless (member proc procs
                      :test (lambda (p1 p2)
                              (and (eql (second p1) (second p2))
                                   (eql (third p1) (third p2)))))
       do (push proc procs))
    (reverse procs)))

;;; (get-processors-without-hyperthreading)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; New lparallel interface
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defvar *worker-id* nil
  "Local ID set for each worker thread when it is created.")

(defun new-kernel (&optional nr-threads)
  (let ((available-procs (get-processors-without-hyperthreading)))
    (ensure nr-threads (length available-procs))
    (end-kernel)
    (when (plusp nr-threads)
      (let ((count -1))
        (flet ((my-worker-context (worker-loop)
                 (let ((*worker-id* (incf count)))
                   ;; set cpu affinity if possible
                   #+sb-cpu-affinity
                   (let* ((proc (nth *worker-id* available-procs))
                          (id (pi-processor proc)))
                     (sb-cpu-affinity:with-cpu-affinity-mask (mask :save t)
                       (sb-cpu-affinity:clear-cpu-affinity-mask mask)
                       (setf (sb-cpu-affinity:cpu-affinity-p id mask) t)))
                   ;; enter the worker loop; return when the worker shuts down
                   (funcall worker-loop))))
          (setf *kernel*
                (make-kernel nr-threads
                             :bindings '((*worker-id* . nil)
                                         (*thread-local-memoization-table* . nil))
                             :context #'my-worker-context)))))))

;;;(pwork (_ (sb-cpu-affinity:get-cpu-affinity-mask)))

#+(or)
(let* ((proc (nth 3 (get-processors)))
       (id (pi-processor proc)))
  id)
;;; start a pool of workers
;;; (fl.parallel::new-kernel 2)
;;; (end-kernel)

(defun pwork (function &optional arguments)
  "Distribute a task to each lparallel worker and wait until all finish.
Arguments may be a vector of size equal or smaller than (kernel-worker-count).
All results are collected in a vector of the same size."
  (assert (not *worker-id*))
  (ensure *kernel* (new-kernel))
  (let* ((nr-workers (kernel-worker-count))
         (arguments (or arguments
                        (make-array nr-workers :initial-element nil)))
         ;; in case fewer arguments than threads have been provided
         (n (min (length arguments) nr-workers))
         (channel (make-channel))
         (from-workers (make-queue))
         (to-workers (make-queue))
         (results (make-array n :initial-element nil)))
    (when (< n (length arguments))
      (error "More work than workers.  Systematic error?"))
    (loop repeat n
          do (submit-task channel (lambda ()
                                    (push-queue t from-workers)
                                    (pop-queue to-workers)
                                    (cons *worker-id*
                                          (apply function
                                                 (aref arguments *worker-id*))))))
    (loop repeat n do (pop-queue from-workers))
    ;; at this point every worker should wait for work
    ;; and the next line sets them doing the work
    (loop repeat n do (push-queue t to-workers))
    ;; and collect the results
    (loop repeat n
          for (id . result) = (receive-result channel) do
            (setf (aref results id) result))
    results))

;;; Avoids many debugger frames popping up, but eliminates
;;; debugging an error at the bottom level 
(setq *debug-tasks-p* nil)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; parallel pools
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun thread-local-memoization-pool (&key (test 'equal))
  (let ((n (if *kernel*
               (kernel-worker-count)
               1)))
    (coerce (loop repeat n collect (make-hash-table :test test)) 'vector)))

(defun thread-local-memoize (pool func args)
  (let ((table (aref pool *worker-id*)))
    (acond ((gethash args table)
            (values it t))
           (t (dbg :memoize "Memoizing for ~A" args)
              (values (setf (gethash args table)
                            (apply func args))
                      nil)))))

(defun test-femlisp-parallel ()
  (get-cpuinfo)
  (get-cachesize)
  (get-processors-without-hyperthreading)
  (pwork (_ *worker-id*))

  (time
   (let* ((channel (make-channel))
          (n 1000)
          (result (make-array n))
          (result-list ()))
     (flet ((test (i)
              (format t "Test: ~D~%" i)
              (force-output t)
              ;;(sleep 1.0)
              (setf (aref result i)
                    i)))
       (flet ((test2 (&rest args)
                (apply #'test args)))
         (loop for i below n do
           (submit-task channel #'test2 i))
         (loop repeat n do
           (push (receive-result channel) result-list))
         (set-difference (coerce result 'list)
                         result-list)
         ))))

  (ignore-errors (pwork (lambda (x) (* x x)) #((1 2) (2))))
  (flet ((f (x) (print (* x x))))
    (let ((channel (make-channel)))
      (loop for i below 10 do
        (submit-task channel #'f i))))
    
  (pwork (lambda (x) (* x x)) #((1) (2)))
  #+(or)
  (let ((n 1000))
    (let ((a1 (fl.matlisp:ones n))
          (a2 (fl.matlisp:ones n))
          (a3 (fl.matlisp:ones n))
          (b1 (fl.matlisp:ones n))
          (b2 (fl.matlisp:ones n))
          (b3 (fl.matlisp:ones n)))
      (time
       (loop for (x y z) across
             (vector (list a1 a2 a3) (list b1 b2 b3)) do
          (fl.matlisp:gemm-tn! 1.0 x y 1.0 z)))
      (time
       (pwork
        (lambda (x y z)
          (fl.matlisp:gemm-tn! 1.0 x y 1.0 z))
        (vector (list a1 a2 a3) (list b1 b2 b3))))))
    
  (pwork (lambda (x) (* x x)) #((1) (2)))
  )
