(in-package :ddo-femlisp)

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

