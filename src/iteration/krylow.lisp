;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; krylow.lisp
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;
;;; Copyright (C) 2003 Nicolas Neuss, University of Heidelberg.
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

(in-package :fl.iteration)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Gradient method
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass <gradient-method> (<linear-iteration>)
  ()
  (:documentation "Gradient-method.  Better use CG."))

(defmethod make-iterator ((linit <gradient-method>) mat)
  "Iterator for the gradient method."
  (let ((p (make-domain-vector-for mat))
	(a (make-image-vector-for mat)))
    (make-instance
     '<iterator>
     :matrix mat
     :initialize nil
     :residual-before t
     :iterate
     #'(lambda (x b r)
	 (declare (ignore b))
	 (copy! r p)
	 (gemm! 1.0 mat p 0.0 a)
	 (let ((lam (/ (dot r p) (dot a p))))
	   (axpy! lam p x)))
     :residual-after nil)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Preconditioned conjugate gradient iteration
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass <cg> (<linear-iteration>)
  ((preconditioner :initform nil :initarg :preconditioner)
   (restart-cycle :reader restart-cycle :initform nil :initarg :restart-cycle))
  (:documentation "Preconditioned conjugate gradient iteration"))

(defmethod make-iterator ((cg <cg>) mat)
  "Standard method for the preconditioned conjugate-gradient iteration."
  (declare (optimize debug))
  (let* ((precond (aand (slot-value cg 'preconditioner)
			(make-iterator it mat))))
    (when (and precond (slot-value precond 'residual-after))
      (error "This preconditioner does not work, because the application
here wants to keep residual and rhs intact."))
    (let ((p (make-domain-vector-for mat))
	  (a (make-image-vector-for mat))
	  (w (make-image-vector-for mat))
	  (q (and precond (make-domain-vector-for mat)))
	  (alpha 0.0)
          (count 0))
      (assert (and p a w))
      (with-slots (initialize iterate) precond
        (flet ((restart (r)
                 (cond
                   (precond
                    (copy! r w)
                    (when initialize
                      (funcall initialize p w w))
                    (funcall iterate p w w))
                   (t (copy! r p)))
                 (setq alpha (dot p r))))
          (make-instance
           '<iterator>
           :matrix mat
           :residual-before t
           :initialize
           #'(lambda (x b r)
               (declare (ignore x b))
               (restart r))
           :iterate
           #'(lambda (x b r)
               (declare (ignore b))
               (when (aand (plusp count)
                           (restart-cycle cg)
                           (zerop (mod count it)))
                 (restart r))
               (unless (zerop alpha)
                 (gemm! 1.0 mat p 0.0 a)
                 (let* ((beta (dot a p))
                        (lam (/ alpha beta)))
                   (axpy! lam p x)
                   (axpy! (- lam) a r)
                   (let ((q (cond (precond
                                   (copy! r w) (x<-0 q)
                                   (funcall iterate q w w)
                                   q)
                                  (t (copy r)))))
                     (let ((new-alpha (dot q r)))
                       (scal! (/ new-alpha alpha) p)
                       (m+! q p)
                       (setq alpha new-alpha))
                     )))
               ;; restart procedure
               (unless (aand (restart-cycle cg)
                             (zerop (mod (incf count) it)))
                 (list :residual-after t)))
           :residual-after t))))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; BiCGStab
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass <bi-cgstab> (<linear-iteration>)
  ((preconditioner :initform nil :initarg :preconditioner))
  (:documentation "Preconditioned Bi-CGStab iteration"))

(defmethod make-iterator ((it <bi-cgstab>) mat)
  "Standard method for the preconditioned conjugate-gradient iteration."
  (let ((precond (whereas ((preconditioner (slot-value it 'preconditioner)))
		   (make-iterator preconditioner mat))))
    (with-slots (iterate) precond
      (make-instance
       '<iterator>
       :matrix mat
       :residual-before t
       :iterate
       (let (r~ p p^ v s s^ rho alpha omega)
	 #'(lambda (x b r)
	     (declare (ignore b))
	     (unless r~ (setf r~ (copy r)))
	     (let ((rho~ (dot r~ r)))
	       (when (zerop rho~)
		 (error "Failure"))
	       (cond
		 (p (axpy! (- omega) v p)
		    (scal! (* (/ rho~ rho) (/ alpha omega)) p)
		    (m+! r p))
		 (t (setf p (copy r))))
	       (setf rho rho~))
	     (ensure p^ (make-analog p))
	     (cond (precond (x<-0 p^) (funcall iterate p^ p p))
		   (t (copy! p p^)))
	     (setf v (m* mat p^))
	     (setf alpha (/ rho (dot r~ v)))
	     (setf s (axpy (- alpha) v r))
	     (ensure s^ (make-analog s))
	     (cond (precond (x<-0 s^) (funcall iterate s^ s s))
		   (t (copy! s s^)))
	     (let ((tee (m* mat s^)))
	       ;; check of (norm s^) dropped
	       (setf omega (/ (dot tee s) (dot tee tee)))
	       ;; update of x
	       (axpy! alpha p^ x)
	       (axpy! omega s^ x)
	       (copy! s r)
	       (axpy! (- omega) tee r)
	       ;; convergence check dropped (is done outside)
	       )
             (list :residual-after t)))
       :residual-after t
       ))))

;;;; Testing

(defun test-krylow ()
  (time (let* ((n 32)
	       (A (fl.matlisp::five-point-stencil-matrix n n))
	       (b (ones (* n n) 1)))
	  (linsolve A b :output t
		    :iteration (make-instance '<cg>))))
  (let ((A #m((2.0 -1.0  0.0)
	      (-1.0  2.0 -1.0)
	      ( 0.0 -1.0  2.0)))
	(b #m((1.0) (2.0) (1.0))))
    (linsolve A b :output t :iteration (make-instance '<cg> :preconditioner (make-instance '<jacobi>)))
    (linsolve A b :output t :iteration (make-instance '<cg>))
    (setf (mref A 2 1) 7.0)
    (linsolve A b :output t :iteration (make-instance '<cg>))  ; diverges
    (linsolve A b :output t :iteration (make-instance '<bi-cgstab> :preconditioner (make-instance '<jacobi>))) ; converges
  )
  )

;;; (test-krylow)
(fl.tests:adjoin-test 'test-krylow)
