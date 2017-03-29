;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Parallel utilities
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(in-package :fl.parallel)

;;;; CPU info

(defun get-cpuinfo ()
  (whereas ((cpuinfo (probe-file "/proc/cpuinfo")))
    (with-open-file (stream cpuinfo)
      (loop for line = (read-line stream nil) while line
            collect (cl-ppcre:split "\\t+: " line)))))

;;; Cachesize

(defun calculate-effective-cachesize ()
  "Calculates effective cachesize in bytes."
  (let ((l 27))
    (flet ((test (k)
               (let* ((n (expt 2 k))
                      (count (expt 2 (- l k)))
                      (x (make-double-float-array n))
                      (y (make-double-float-array n)))
                 (declare (optimize speed (safety 0)))
                 (measure-time
                  (lambda () (replace x y)) count))))
      (format t "Measuring effective cache size (this may take some time)...~%")
      (loop for k from 10 below 22
            for previous = nil then next
            and next = (test k)
            do
               (format t "~2D ~A~%" k next)
               (force-output)
               (when (and previous (> next (* 1.2 previous)))
                 (loop-finish))
            finally (return (expt 2 (+ k 3)))))))

(defun get-cachesize ()
  (aif (aand (get-cpuinfo)
             (assoc "cache size" it :test #'string=))
       (let* ((size (second it))
              (sep-pos (position #\Space size))
              (number (parse-integer (subseq size 0 sep-pos)))
              (unit (subseq size (1+ sep-pos))))
         (* number
            (stringcase unit
              ("KB" 1024)
              (t (error "Unknown unit: ~A" unit)))))
    (calculate-effective-cachesize)))

(defvar *cachesize* (get-cachesize))

(unless (<= (expt 2 19) *cachesize* 32e6)
  (warn "Possible problem: cachesize not in standard range"))

;;; Processors

(defstruct (procinfo (:conc-name pi-) (:type list))
  cpu
  core
  socket
  node)

(defun get-processors ()
  (whereas ((output (or (ignore-errors (fl.port:run-program-output "lscpu" '("-p")))
                        ;; the following is a rather unsafe kludge for
                        ;; getting around an SBCL problem occuring
                        ;; when calling subprocesses when the
                        ;; available memory is small.
                        #+(or) 
                        (awhen (probe-file #p"femlisp:external;lscpu-e-output")
                          (with-open-file (stream it)
                            (loop for line = (read-line stream nil)
                                  while line collect line))))))
    (loop for line in output while line
          unless (eql (aref line 0) #\#)
          collect
          (mapcar #'parse-integer
                  (take 4 (cl-ppcre:split "," line))))))

(defun allowed-processors ()
  (if (member :cl-cpu-affinity *features*)
      (progn
        #+cl-cpu-affinity
        (cl-cpu-affinity:with-cpu-affinity-mask (mask)
          (loop for proc in (get-processors)
                when (cl-cpu-affinity:cpu-affinity-p (pi-cpu proc) mask)
                  collect proc)))
      (get-processors)))

(defun get-workers ()
  (group-by #'pi-core (allowed-processors)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; New lparallel interface
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defvar *worker-id* nil
  "Local ID set for each worker thread when it is created.")

(defun new-kernel (&optional nr-threads (set-affinity-p t))
  (let* ((available-workers (get-workers))
         (max-workers (length available-workers)))
    (ensure nr-threads (if available-workers
                           max-workers
                           1))
    (end-kernel)
    (when set-affinity-p
      (unless (member :cl-cpu-affinity *features*)
        (warn "No CPU pinning possible on this architecture.")
        (setq set-affinity-p nil)))
    (when (< nr-threads 1)
      (error "Number of worker threads must be positive"))
    (if (null available-workers)
        (when set-affinity-p
          (warn "No worker information available, so no CPU pinning is possible")
          (setq set-affinity-p nil))
        (when (> nr-threads max-workers)
          (warn "~D workers required, but only ~D CPUs are available, which may slow down calculations.
~:[~;  Furthermore no CPU pinning is possible in this situation.~]"
                nr-threads max-workers set-affinity-p)
          (setq set-affinity-p nil)))
    (when (= nr-threads 1)
      (warn "Only one worker thread is required, which does not really help with performance.
~:[~;  We also do not do CPU pinning in this situation, but use all available CPUs.~]"
            set-affinity-p)
      (setq set-affinity-p nil))
    (let ((count -1))
      (flet ((my-worker-context (worker-loop)
               (let ((*worker-id* (incf count)))
                 ;; set cpu affinity
                 (when set-affinity-p
                   (funcall (intern "SET-CPU-AFFINITIES" (find-package "CL-CPU-AFFINITY"))
                            (mapcar #'pi-cpu (nth *worker-id* available-workers))))
                 ;; enter the worker loop; return when the worker shuts down
                 (funcall worker-loop))))
        (setf *kernel*
              (make-kernel nr-threads
                           :bindings '((*worker-id* . nil)
                                       (*thread-local-memoization-table* . nil))
                           :context #'my-worker-context))))))

;;; starting and ending a pool of workers
;;; (fl.parallel::new-kernel)
;;; (end-kernel)

;;;(pwork (_ (cl-cpu-affinity:cpu-affinity-mask-string)))

(defun pwork (function &optional arguments)
  "Distribute a task to each lparallel worker and wait until all finish.
Arguments may be a vector of size equal or smaller than (kernel-worker-count).
All results are collected in a vector of the same size."
  (assert (not *worker-id*))
  (assert *kernel* () "No kernel")
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
;;; debugging an error inside threads
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
  (let ((table (aref pool (or *worker-id* 0))))
    (acond ((gethash args table)
            (values it t))
           (t (dbg :memoize "Memoizing for ~A" args)
              (values (setf (gethash args table)
                            (apply func args))
                      nil)))))

(defun test-femlisp-parallel ()
  (get-cpuinfo)
  (get-cachesize)
  (get-processors)
  (allowed-processors)
  (loop for proc in (get-processors) do
    (format t "~A~%" proc))
  (get-workers)
  (pwork (_ *worker-id*))
  #+cl-cpu-affinity (pwork (_ (cl-cpu-affinity:cpu-affinity-mask-string)))

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
  ;; a useful test which is uncommented because FL.MATLISP is not
  ;; available at this time
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
