(in-package :femlisp-ddo)

(defmethod fl.strategy::approximate :around ((fe-strategy <stationary-fe-strategy>) blackboard)
  (let ((*communication-real-time* 0.0)
        (*synchronization-real-time* 0.0)
        (*communication-size* 0))
    (call-next-method)
    (dbg-when :synchronization-time-in-approximation-step
      (format t "~&Synchronization total time: ~F seconds~%"
              *synchronization-real-time*)
      (format t "~&Communication total time: ~F seconds~%"
              *communication-real-time*)
      (format t "~&Communication size: ~D numbers~%"
              *communication-size*))))

(defvar *iter-report-ranks* t
  "List of ranks for which reports are shown.  T means all ranks.")

(defmethod fl.iteration::iter-format :around
    (stream control-string &rest args)
  "Switch off reporting for some processors."
  (declare (ignore args)) 
  (let ((rank (and (mpi-initialized) (mpi-comm-rank))))
    (when (or (null rank)
              (eq *iter-report-ranks* T)
              (member rank *iter-report-ranks*))
          (call-next-method))))

  
  
