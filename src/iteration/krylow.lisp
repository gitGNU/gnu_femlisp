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

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Gradient method
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass <gradient-method> (<linear-iteration>)
  ()
  (:documentation "Gradient-method.  Better use CG."))

(defmethod make-iterator ((linit <gradient-method>) mat)
  "Gradient-method via matlisp methods."
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
	 (copy! r p) (copy! p a)
	 (setq a (m* mat a))
	 (let ((lam (/ (dot r p) (dot a p))))
	   (axpy! lam p x)))
     :residual-after nil)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; CG method without preconditioning
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass <cg> (<linear-iteration>)
  ()
  (:documentation "CG iteration.  Think about using preconditioning,
i.e. PCG."))

(defmethod make-iterator ((linit <cg>) mat)
  "Standard method for CG iteration.  Works if matlisp methods are defined for
the given matrix-vector representation."
  (let ((p (make-domain-vector-for mat))
	(a (make-image-vector-for mat))
	(defnorm 0.0))
    (make-instance
     '<iterator>
     :matrix mat
     :residual-before t
     :initialize
     #'(lambda (x b r)
	 (declare (ignore x b))
	 (copy! r p)
	 (setq defnorm (dot r r)))
     :iterate
     #'(lambda (x b r)
	 (declare (ignore b))
	 (unless (zerop defnorm)
	   (copy! p a)
	   (setq a (m* mat a))
	   (let ((lam (/ defnorm (dot a p))))
	     (axpy! lam p x)
	     (axpy! (- lam) a r)
	     (let ((new-defnorm (dot r r)))
	       (scal! (/ new-defnorm defnorm) p)
	       (m+! r p)
	       (setq defnorm new-defnorm)))))
     :residual-after t)))

(defparameter *standard-cg* (make-instance '<cg>)
  "CG iteration.")

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; PCG iteration
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass <pcg> (<linear-iteration>)
  ((preconditioner :initform nil :initarg :preconditioner))
  (:documentation "Preconditioned CG iteration"))

(defparameter *standard-pcg* (make-instance '<pcg>)
  "PCG iteration.")

(defmethod make-iterator ((pcg <pcg>) mat)
  "Standard method for PCG iteration.  Works if matlisp methods are defined
for the given matrix-vector representation."
  (let ((p (make-domain-vector-for mat))
	(a (make-image-vector-for mat))
	(q (make-domain-vector-for mat))
	(alpha 0.0)
	(precond (whereas ((preconditioner (slot-value pcg 'preconditioner)))
		   (make-iterator preconditioner mat))))
    (with-slots (initialize iterate) precond
      (make-instance
       '<iterator>
       :matrix mat
       :residual-before t
       :initialize
       #'(lambda (x b r)
	   (declare (ignore x b))
	   (cond
	     (precond
	      (when initialize (funcall initialize p r r))
	      (funcall iterate p r r))
	     (t (copy! r p)))
	   (setq alpha (dot p r)))
       :iterate
       #'(lambda (x b r)
	   (declare (ignore b))
	   (unless (zerop alpha)
	     (copy! p a)
	     (setq a (m* mat a))
	     (let* ((beta (dot a p))
		    (lam (/ alpha beta)))
	       (axpy! lam p x)
	       (axpy! (- lam) a r)
	       (cond (precond (x<-0 q) (funcall iterate q r r))
		     (t (copy! r q)))
	       (let ((new-alpha (dot q r)))
		 (scal! (/ new-alpha alpha) p)
		 (m+! q p)
		 (setq alpha new-alpha)))))
       :residual-after t))))