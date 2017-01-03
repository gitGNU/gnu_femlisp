(in-package :ddo-femlisp-test)

(defparameter *n*
  (mpi-comm-size))

(format t "CPU affinity: ~A~%" (princ-to-string (sb-cpu-affinity:get-cpu-affinity-mask)))
(force-output)

;;(progn #+(or) (fl.parallel::end-kernel) (fl.parallel::new-kernel) nil)

#+(or)
(pwork (_ (format t "CPU affinity: ~A~%" (princ-to-string (sb-cpu-affinity:get-cpu-affinity-mask)))
          (force-output)))

(dbg-off :communication)
(dbg-off :distribute)
(dbg-on :log-times)
(dbg-on :synchronization-time-in-approximation-step)
(setq *debug-show-data* nil)

;;; actual calculation
(format t "~D- starting calculation~%" (mpi-comm-rank)) (force-output)
(elahom-calculation 4 1)

(sb-ext:exit)
