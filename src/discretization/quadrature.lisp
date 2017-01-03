;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; quadrature.lisp
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

(in-package :fl.discretization)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; integration points <ip>
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defstruct (<ip> (:conc-name ip-))
  (weight (required-argument) :type double-float)
  (coords (required-argument) :type double-vec))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; integration rule
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass <integration-rule> ()
  ((cell :initarg :cell :type <cell>)
   (order :initarg :order :type fixnum)
   (points :reader integration-points :initarg :points)
   (weights :reader integration-weights :initarg :weights)))

;;; Later on, integration rules with different qualities could be
;;; registered and chosen on demand.

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; integration of polynomials
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun integrate-monomial (mono dim)
  (declare (type list mono))
  (declare (type fixnum dim))
  (/ (factorial mono)
     (factorial (+ dim (reduce #'+ mono)))))

(defun integrate-over-reference-simplex (poly n)
  "Integrates a polynomial over the reference simplex.  This is done by
splitting into monomials and computing
 \int_S(n) x^alpha = (alpha1! alpha2! ...) / (n+alpha1+...)!"
  (loop for (coeff . mono) in (split-into-monomials (coefficients poly))
	summing (* coeff (integrate-monomial mono n))))

(defun integrate-monomial-over-simplex-product (mono dims)
  (reduce #'* (mapcar #'(lambda (factor-mono n) (integrate-monomial factor-mono n))
		      (split-by-length
                       (append mono (make-list (- (reduce #'+ dims) (length mono))
                                               :initial-element 0))
                       dims)
		 dims)))

(defun integrate-over-reference-product-cell (poly factor-dims)
  (loop for c&m in (split-into-monomials (coefficients poly))
	summing (* (car c&m)
		  (integrate-monomial-over-simplex-product (cdr c&m) factor-dims))))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; computation of quadrature formulae
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; Jacobi polynomials for alpha, beta, n
(defun jacobi-polynomial (alpha beta s)
  (let* ((a+b+2 (+ alpha beta 2))
	 (J_0 (make-polynomial (list 1)))
	 (J_1 (make-polynomial (list (/ (- alpha beta) a+b+2) 1))))
    (cond
      ((= s 0) J_0)
      ((= s 1) J_1)
      (t (labels ((new-jacobi (J_n-2 J_n-1 n)
		    (let* ((2n+a+b (+ (* 2 n) alpha beta))
			   (n-1 (1- n))
			   (c_n (/ (* 4 n-1 (+ n-1 alpha) (+ n-1 beta) (+ n-1 alpha beta))
				   (* (- 2n+a+b 1) (- 2n+a+b 2) (- 2n+a+b 2) (- 2n+a+b 3))))
			   (J_n (axpy (- c_n) J_n-2
				      (poly* (make-polynomial (list (/ (* (+ alpha beta) (- alpha beta))
								       (* 2n+a+b (- 2n+a+b 2))) 1))
					     J_n-1))))
		      (if (= n s)
			  J_n
			  (new-jacobi J_n-1 J_n (1+ n))))))
	   (new-jacobi J_0 J_1 2))))))

(defun legendre-polynomial (n)
  (jacobi-polynomial 0 0 n))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; roots of separating family
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun zeros-of-separating-family (family n accuracy &optional (float-p t))
  "This function works on the assumption that p_0=const and the zeros of
p_n-1 together with the interval boundaries separate the zeros of p_n.
Then an interval method is performed.  Of course, accuracy problems may
occur for the inexact arithmetic."
  (cond
    ((minusp n) (error "unknown family"))
    ((zerop n) '())
    (t (loop for seps = (append (list (if float-p -1.0 -1))
                                (zeros-of-separating-family
                                 family (1- n) accuracy float-p)
                                (list (if float-p 1.0 1)))
	     then (cdr seps)
	     until (single? seps)
	     collect (interval-method
		      (funcall family n) (first seps) (second seps) accuracy)))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Gauss quadrature
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; s-point Gauss-formula for the weight (1+y)^beta
(defun gauss-points-for-weight (beta s &optional (float-p t))
  (zeros-of-separating-family
   #'(lambda (n) (jacobi-polynomial 0 beta n)) s 1.0e-16 float-p))

(defun weights-for-ips (beta ips)
  "Determines weights for the integration points @arg{ips} such
that they integrate
  @math{int_{-1}^{1} (1+y)^beta ... dy}
 with optimal order."
  (let ((int-weight (poly-expt (make-polynomial '(1 1)) beta)))
    (labels ((weight-for-ip (xi)
	       (let* ((ips-without-xi (remove xi ips))
		      (lagrange (scal
				 (/ 1 (reduce #'* (mapcar #'(lambda (xj) (- xi xj)) ips-without-xi)))
				 (reduce #'poly* (mapcar #'(lambda (xj)
							     (make-polynomial (list (- xj) 1)))
							 ips-without-xi)
					 :initial-value (make-polynomial '(1)))))
		      (lag-int (integrate-simple-polynomial (poly* lagrange int-weight))))
		 (- (evaluate lag-int 1) (evaluate lag-int -1)))))
      (mapcar #'weight-for-ip ips))))


;;; qr for \int_0^1 (1-y)^n f(y) dy
(defun gauss-rule-for-weight (n s)
  (if (minusp n)
      ()
      (loop with coords = (gauss-points-for-weight n s)
	    for coord in coords
	    and weight in (weights-for-ips n coords)
	    collect (list (float (/ weight (expt 2 (+ n 1))) 1.0)
			  (/ (- 1.0 coord) 2)))))

;;; gauss rule for an s^n-point method on the n-simplex
;;; result = ( (weight coordinates) ... )
(defun gauss-rule-for-simplex (n s)
  (labels ((transform (pos factor) ; transformation cube->simplex
	     (if (null pos) '()
		 (cons (* (car pos) factor)
		       (transform (cdr pos) (* factor (- 1.0 (car pos))))))))
    ;; loop through nodes
    (apply #'map-product
	   #'(lambda (&rest args)
	       (make-<ip>
		:weight (reduce #'* (mapcar #'car args)) ; multiply weights
		:coords (coerce (transform (mappend #'cdr args) 1.0) 'double-vec)))
	   (mapcar #'(lambda (k) (gauss-rule-for-weight k s))
		   (loop for k from (1- n) downto 0 collect k)))))

(defun product-rule (&rest quadrature-rules)
  "Computes a product rule for several lower-dimensional quadrature rules."
  (apply #'map-product
	 (lambda (&rest args)
	   (make-<ip>
	    :weight (apply #'* (mapcar #'ip-weight args))
	    :coords (apply #'concatenate 'double-vec (mapcar #'ip-coords args))))
	 quadrature-rules))

(with-memoization (:id 'gauss-rule)
  (defun gauss-rule (factor-dims s)
    "Returns an s-point Gauss integration rule."
    (memoizing-let ((factor-dims factor-dims) (s s))
      (let ((ips (if (null factor-dims)
		     (list (make-<ip> :weight 1.0 :coords (double-vec)))
		     (apply #'product-rule
			    (mapcar #'(lambda (dim) (gauss-rule-for-simplex dim s))
				    factor-dims)))))
	(make-instance
	 '<integration-rule>
	 :order (* 2 s)
	 :points (map 'vector #'ip-coords ips)
	 :weights (map 'vector #'ip-weight ips))))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Gauss-Lobatto quadrature
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun gauss-lobatto-family (s)
  (differentiate (legendre-polynomial (1+ s)) 0))

(defun gauss-lobatto-points (s)
  `(-1.0 ,@(zeros-of-separating-family #'gauss-lobatto-family s 1.0e-16) 1.0))

(defun gauss-lobatto-points-on-unit-interval (s)
  (mapcar #'(lambda (x) (/ (1+ x) 2)) (gauss-lobatto-points s)))

(defun gauss-lobatto-rule-on-unit-interval (s)
  (let* ((points (gauss-lobatto-points s))
         (weights (weights-for-ips 0 points)))
    (loop for point in points
          and weight in weights
          collect (make-<ip>
                   :weight (/ weight 2)
                   :coords (double-vec (/ (+ 1.0 point) 2))))))

(with-memoization (:id 'gauss-lobatto-rule)
  (defun gauss-lobatto-rule (factor-dims s)
    (assert (apply #'= 1 factor-dims)
            () "Should be checked theoretically first, if one can do this analogous to the Gauss rules.")
    (memoizing-let ((factor-dims factor-dims) (s s))
      (let ((ips (if (null factor-dims)
		     (list (make-<ip> :weight 1.0 :coords (double-vec)))
		     (apply #'product-rule
			    (make-list (length factor-dims) :initial-element
                                       (gauss-lobatto-rule-on-unit-interval s))))))
        (make-instance
         '<integration-rule>
         :order (- (* 2 (+ s 2)) 2)
         :points (map 'vector #'ip-coords ips)
         :weights (map 'vector #'ip-weight ips))))))

;;;; calculate order

(defun test-integration-rule (points weights &optional (output t))
  "Tests the given rule against monomials of increasing degree.
As a by-product this determines the order of the rule."
  (let ((dim (length (elt points 0))))
    (loop
      with flag = nil
      for degree from 0
      until flag do
        (when output (format t "Order ~D:~%" degree))
        (loop for mono in (n-variate-monomials-of-degree dim degree) do
          (let* ((exact (fl.discretization::integrate-over-reference-product-cell
                         mono (factor-dimensions (n-cube dim))))
                 (approximate
                   (loop for point across points
                         and weight across weights
                         sum (* weight (evaluate mono point))))
                 (diff (- exact approximate)))
            (unless (mzerop diff 1e-12)
              (setq flag t))
            (format t "~10A ~10A ~10,2G~%"
                    exact approximate diff)))
      finally (return (- degree 2)))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Testing: (test-quadrature)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun test-quadrature ()
  (interval-method (make-polynomial '(0.1 1)) -1.0 1.0 1e-5)
  (gauss-rule () 2)
  (time (let () (gauss-rule-for-weight 2 10) nil)) ; 0.5-0.8 on toba
  (gauss-points-for-weight 2 6)
  
  (gauss-rule-for-simplex 2 3)
  (gauss-lobatto-points 2)
  (gauss-lobatto-points-on-unit-interval 2)
  (gauss-lobatto-points-on-unit-interval 10)

  (length (integration-points (gauss-rule '(1 2) 20)))
  (let* ((n 8)
         (a (time (mapcar (_ (float _ 1.0)) (gauss-points-for-weight 0 n nil))))
         (b (time (gauss-points-for-weight 0 n t))))
    (mapcar #'- a b))
  (let ((rule (gauss-rule '(1) 2)))
    (test-integration-rule (integration-points rule) (integration-weights rule)))

  (let ((rule (gauss-lobatto-rule '(1) 1)))
    (test-integration-rule (integration-points rule)
                           (integration-weights rule)))
  )

(fl.tests:adjoin-test 'test-quadrature)
