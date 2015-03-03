;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; parcells.lisp - parallel computation in a network
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;
;;; Copyright (C) 2009 Nicolas Neuss, University of Karlsruhe.
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

(defclass network (work-group)
  ()
  (:documentation "Improvement of work-group.  Will maybe replace it later
  on."))

(defparameter *network* nil
  "This variable is rebound in WITH-PARALLEL-NETWORK to a local computing
  network.  A value of NIL means no parallel calculation.")

(defgeneric value (obj)
  (:method ((obj t))
    "All objects but cells are self-evaluating."
    obj))

(defclass parcell (mutex-mixin)
  ((network :reader network :initarg :network)
   (state :accessor state :initform :waiting
          :type (member :waiting :computing :ready))
   (value :accessor value)
   (name :reader name :initform nil :initarg :name)
   (action :reader action :initarg :action :documentation
           "Action guaranteed to compute the value.")
   (inputs :accessor inputs :initform ()
           :documentation "Missing input parcells for computing the value.")
   (outputs :accessor outputs :initform ()
            :documentation "Parcells to which information about a computed value
            has to be propagated.")))

(defmethod value :before ((cell parcell))
  (wait (network cell) :until (_ (value-computed-p cell))))

(defun connect (depends on)
  (with-mutual-exclusion (depends)
    (with-mutual-exclusion (on)
      (pushnew depends (outputs on))
      (pushnew on (inputs depends)))))

(defun disconnect (depends on)
  (with-mutual-exclusion (depends)
    (with-mutual-exclusion (on)
      (deletef depends (outputs on))
      (deletef on (inputs depends)))))

(defun parcell-p (obj)
  (typep obj 'parcell))

(defgeneric value-computed-p (obj)
  (:method (obj) t)
  (:method ((cell parcell))
    (eq (state cell) :ready)))

(defgeneric waiting-p (obj)
  (:method (obj) nil)
  (:method ((cell parcell))
    (eq (state cell) :waiting)))

(defgeneric computable (obj)
  (:method (obj) nil)
  (:method ((cell parcell))
    (with-mutual-exclusion (cell)
      (and (eq (state cell) :waiting)
           (every #'value-computed-p (inputs cell))))))

(defun propagate (cell)
  "If it is computable, push the cell to the tasks to be handled."
  (with-mutual-exclusion (cell)
    (unless (computable cell)
      (return-from propagate))
    (enqueue cell (tasks (network cell)))
    (setf (state cell) :computing)))

(defun handle-cell (cell)
  (let ((network (network cell)))
    (mp-dbg "~&Performing action from parcell ~A~%" (name cell))
    (let ((result (funcall (action cell))))
      (with-mutual-exclusion (cell)
        (setf (slot-value cell 'value) result)
        (setf (state cell) :ready)
        (mp-dbg "~&Computed: ~A.~%" result)))
    (loop for dest in (outputs cell) do
         (disconnect dest cell)
         (propagate dest))
    (notify network)))

(defun make-parcell (&key (name (gensym)) (network *network*) action inputs)
  (lret ((cell (make-instance 'parcell :network network :name name :action action)))
    (cond (network
           (with-mutual-exclusion (cell)
             (loop for source in inputs
                when (parcell-p source) do
                (with-mutual-exclusion (source)
                  (unless (value-computed-p source)
                    (connect cell source))))
             (propagate cell)))
          (t ; no parallel calculation in this case
           (setf (slot-value cell 'value) (funcall action)
                 (slot-value cell 'state) :ready)))))

(defmacro named-parcell (name inputs &body body)
  "Generate a parcell depending on the given input parameters (which can be
parcells again)."
  (with-gensyms (args)
    `(let ((,args (list ,@(mapcar #'second inputs))))
       (make-parcell :name ,name
                     :action (lambda ()
                               (apply (lambda ,(mapcar #'first inputs) ,@body)
                                      (mapcar #'value ,args)))
                     :inputs ,args))))

(defmacro parcell (inputs &body body)
  "Generate an unnamed parcell."
  `(named-parcell nil ,inputs ,@body))

(defun execute-in-parallel-network (func &key number-of-threads)
  (if number-of-threads
      (let ((*network* (make-instance 'network :work #'handle-cell)))
        (loop repeat number-of-threads do
             (add-worker *network*))
        (unwind-protect (funcall func)
          (finish (tasks *network*))
          (wait *network* :while #'threads)))
      (funcall func)))

(defmacro with-parallel-network ((&key number-of-threads) &body body)
  `(execute-in-parallel-network
    (lambda () ,@body) :number-of-threads ,number-of-threads))

(defun calculate-local-matrix-slowly (cell)
  (declare (ignore cell))
  (sleep 1.0)
  1.0)

(defun sum-long-calculations (n)
  (let ((result 0.0))
    (loop for i below n
       for local-mat = (parcell () (calculate-local-matrix-slowly i))
       do (setf result (parcell ((lmat local-mat)
                                 (result result))
                         (+ lmat result))))
    (value result)))


;;;; Testing: (test-parcells)

(defun test-parcells ()

  (with-parallel-network (:number-of-threads nil)
    
    (parcell () (calculate-local-matrix-slowly 1))
    
    (time
     (let* ((c1 (make-parcell :name 'c1 :action (_ (calculate-local-matrix-slowly 1))))
            (c2 (make-parcell :name 'c2 :action (_ (calculate-local-matrix-slowly 2))))
            (c3 (make-parcell :name 'c3 :action (_ (+ (value c1) (value c2)))
                                        :inputs (list c1 c2)))
            (c4 (make-parcell :name 'c4 :action (_ (calculate-local-matrix-slowly 2))))
            (c5 (make-parcell :name 'c5 :action (_ (+ (value c3) (value c4)))
                                        :inputs (list c3 c4))))
       (value c5)))

    (time 
     (let* ((c1 (parcell ()
                  (calculate-local-matrix-slowly 1)))
            (c2 (parcell ()
                  (calculate-local-matrix-slowly 2)))
            (c3 (parcell ((c1 c1)
                          (c2 c2))
                  (+ c1 c2))))
       (value c3)))
    
    (time (sum-long-calculations 10))
    )
  
  (assert (null (femlisp-workers)))
  
  )
