;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; polynom.lisp
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

(in-package :algebra)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; abstract ring class, perhaps we'll need it later separately
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; zero and unit for numbers
(defmethod zero ((f number)) 0)
(defmethod unit ((f number)) 1)
(defmethod zero? ((x number)) (zerop x))
(defmethod unit? ((x number)) (= x 1))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; class of multivariate polynomials (over arbitrary rings)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defstruct (<polynomial> (:conc-name poly-) (:print-function print-polynomial))
  (coeffs nil :type list))

(defstruct (<polynomial> (:conc-name poly-) (:print-function print-polynomial))
  (coeffs nil :type list))

(defun shift-polynomial (poly dim)
  "Shifts a polynomial in dimension, e.g. x_1 becomes x_2."
  (labels ((nest (lst k)
	     (if (zerop k)
		 lst
		 (list (nest lst (1- k))))))
    (make-polynomial (nest (poly-coeffs poly) dim))))

(defun polynomial? (x) (typep x '<polynomial>))

(defmethod degree ((poly <polynomial>))
  (- (length (poly-coeffs poly)) 1))

(defmethod multi-degree ((poly <polynomial>))
  (labels ((tensorial-degree (coeff-lst)
	     (apply #'max
		    (mapcar #'(lambda (coeff deg)
				(if (listp coeff)
				    (+ deg (tensorial-degree coeff))
				    deg))
			    coeff-lst
			    (make-set 0 (length coeff-lst))))))
    (tensorial-degree (poly-coeffs poly))))

(defmethod partial-degree ((lst list) (index fixnum))
  (if (zerop index)
      (- (length lst) 1)
      (apply #'max -1
	     (mapcar #'(lambda (term)
		    (partial-degree term (- index 1)))
		  lst))))

(defmethod partial-degree ((poly <polynomial>) (index fixnum))
  (partial-degree (poly-coeffs poly) index))

;;; Returns number of variables on which the polynomial depends.  Does
;;; the name variance fit?
(defmethod variance ((poly <polynomial>))
  (labels ((rec (index coeffs)
	     (cond
	       ((null coeffs) (+ index 1))
	       ((listp coeffs) (rec (+ index 1) (car coeffs)))
	       (t index))))
    (rec 0 (poly-coeffs poly))))


;;; returns the maximal partial degree of a polynomial
(defmethod maximal-partial-degree ((poly <polynomial>))
  (loop for k from 0 below (variance poly)
	maximize (partial-degree poly k)))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; <polynomial> constructor (from coefficient list)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; important, since redundancy may be introduced during arithmetic
;;; operations we have to strip trailing zeros from a coefficient list
(defun simplify (lst)
  (labels ((empty? (lst)
	     (or (null lst)
		 (and (null (cdr lst))
		      (listp (car lst))
		      (empty? (car lst))))))
    (if (empty? lst)
	lst
	(let ((rest (simplify (cdr lst))))
	  (if (empty? rest)
	      (if (listp (car lst))
		  (list (simplify (car lst)))
		  (if (zero? (car lst))
		      rest
		      (cons (car lst) rest)))
	      (if (listp (car lst))
		  (cons (simplify (car lst)) rest)
		  (cons (car lst) rest)))))))

(defun make-polynomial (lst)
  (make-<polynomial> :coeffs (simplify lst)))

(defun eliminate-small-coefficients (poly &optional (threshold 1.0e-12))
  (make-polynomial
   (subst-if 0.0d0 #'(lambda (x) (and (numberp x) (< (abs x) threshold)))
	     (poly-coeffs poly))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; checks
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmethod zero? ((lst list)) (null lst))
(defmethod zero? ((poly <polynomial>))
  (zero? (poly-coeffs poly)))
(defmethod unit? ((lst list)) (and (single? lst) (unit? (car lst))))
(defmethod unit? ((poly <polynomial>))
  (unit? (poly-coeffs poly)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; write method
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; splits the coefficient list of an n-variate polynomial into
;;; monomials given by their exponents (1 1) = x1*x2
(defun split-into-monomials (coeffs)
  (labels ((create-list (lst mono)
	     (if (null lst)
		 '()
		 (append
		  (if (listp (car lst))
		      (create-list (car lst) (cons 0 mono))
		      (list (cons (car lst) mono)))
		  (create-list (cdr lst)
			       (cons (1+ (car mono)) (cdr mono)))))))
    (create-list coeffs '(0))))

(defun write-monomial (coeff&monomial stream)
  (let ((coeff (car coeff&monomial))
	(mono (cdr coeff&monomial)))
    
    ;; princ coefficient
    (if (not (unit? coeff))		; was: (and (number? coeff) (= 1 coeff)))
	(progn
	  (format stream " ~,5E " coeff)
	  (if (not (every #'zerop mono))
	      (princ "*" stream)))
	(if (every #'zerop mono)
	    (princ "1" stream)))
  
    (if (not (every #'zerop mono))
	(do ((i 1 (1+ i))
	     (mono (reverse mono) (cdr mono)))
	    ((null mono))
	  (unless (zerop (car mono))
	    (princ "x" stream)
	    (princ i stream)
	    (when (> (car mono) 1)
	      (princ "^" stream)
	      (princ (car mono) stream))
	    (when (and (not (single? mono)) (not (every #'zerop (cdr mono))))
	      (princ "*" stream)))))))

(defmethod print-polynomial (poly stream depth)
  (declare (ignore depth))
  (let ((*print-circle* nil))
    (print-unreadable-object
     (poly stream :type t :identity t)
     (format stream "{")
     (loop with start-flag = t
	   for coeff&mono in (split-into-monomials (poly-coeffs poly))
	   unless (zerop (car coeff&mono)) do
	   (if start-flag
	       (setq start-flag nil)
	       (princ " + " stream))
	   (write-monomial coeff&mono stream))
     (format stream " }"))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; <function> methods
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; evaluation of a multivariate polynomial
(defun poly-eval (coeffs coords)
  (cond ((not (listp coeffs)) coeffs)
	((null coords) (or (car coeffs) 0.0d0))
	((null coeffs) 0.0d0)
	(t (+ (poly-eval (car coeffs) (cdr coords))
	      (* (car coords)
		 (poly-eval (cdr coeffs) coords))))))

(defmethod evaluate ((poly <polynomial>) (x list))
  (poly-eval (poly-coeffs poly) x))
(defmethod evaluate ((poly <polynomial>) (x array))
  (poly-eval (poly-coeffs poly) (vector->list x)))
(defmethod evaluate ((poly <polynomial>) (x number))
  (poly-eval (poly-coeffs poly) (list x)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; neutral cells for ring operations
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; generates a list that has the same nestedness asthenon lst and contains
;;; cell as single cell
(defun make-inner (lst cell)
  (if (consp lst)
      (list (make-inner (car lst) cell))
      cell))

(defmethod zero ((f list)) '())
(defmethod unit ((f list)) '(1))
(defmethod zero ((f <polynomial>))
  (make-polynomial
   (make-inner (cdr (poly-coeffs f)) '())))

(defmethod unit ((f <polynomial>))
  (make-polynomial
   (make-inner f '(1))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; <vector> methods
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


(defmethod vec+ ((x number) (y number)) (+ x y))
  
(defmethod vec+ ((x list) (y list))
  (labels ((rec (x y)
	     (cond ((null x) y)
		   ((null y) x)
		   (t (cons (vec+ (car x) (car y))
			    (rec (cdr x) (cdr y)))))))
    (rec x y)))

(defmethod vec+ ((x <polynomial>) (y <polynomial>))
  (make-polynomial (vec+ (poly-coeffs x) (poly-coeffs y))))
(defmethod vec+ ((x <polynomial>) (y number))
  (make-polynomial (vec+ (poly-coeffs x) (list y))))
(defmethod vec+ ((x number) (y <polynomial>))
  (make-polynomial (vec+ (list x) (poly-coeffs y))))

;;;;(defmethod vec- (x y) (vec+ x (vec-s* -1 y)))
;(defmethod vec- ((x list) (y list))
;  (labels ((rec (x y)
;             (cond ((null x) (vec-s* -1 y))
;                   ((null y) x)
;                   (t (cons (- (car x) (car y))
;                            (rec (cdr x) (cdr y)))))))
;    (rec x y)))
;
;
;(defmethod vec- ((x <polynomial>) (y <polynomial>))
;  (make-polynomial (vec- (poly-coeffs x) (poly-coeffs y))))
;(defmethod vec- ((x <polynomial>) (y number))
;  (make-polynomial (vec- (poly-coeffs x) (list y))))
;(defmethod vec- ((x number) (y <polynomial>))
;  (make-polynomial (vec- (list x) (poly-coeffs y))))
;
(defmethod vec-s* ((val number) (x number)) (* val x))
(defmethod vec-s* ((val number) (vec list))
  (mapcar #'(lambda (x) (vec-s* val x)) vec))
(defmethod vec-s* ((val number) (poly <polynomial>))
  (make-polynomial
   (vec-s* val (poly-coeffs poly))))

;;; These methods are quite special for our use in polynomial
;;; multiplication.  This is a hint that things should be done in a
;;; better way.
(defmethod vec+ ((x number) (y list))
  (if (null y) x (cons (vec+ (car y) x) (cdr y))))
(defmethod vec+ ((y list) (x number))
  (if (null y) x (cons (vec+ (car y) x) (cdr y))))

(defmethod poly* ((f number) (g number)) (* f g))
(defmethod poly* ((f number) (g list)) (vec-s* f g))
(defmethod poly* ((g list) (f number)) (vec-s* f g))
(defmethod poly* ((f list) (g list))
  (if (null g)
      '()
      (vec+ (mapcar #'(lambda (x) (poly* (car g) x)) f)
	    (cons 0
		  (poly* f (cdr g))))))

(defmethod poly* ((f <polynomial>) (g <polynomial>))
  (make-polynomial (poly* (poly-coeffs f) (poly-coeffs g))))
(defmethod poly* ((f number) (g <polynomial>)) (vec-s* f g))
(defmethod poly* ((g <polynomial>) (f number)) (vec-s* f g))

(defmethod vec-sp ((f <polynomial>) (g <polynomial>)) (poly* f g))
(defmethod vec-sp ((f number) (g <polynomial>)) (poly* f g))
(defmethod vec-sp ((f <polynomial>) (g number)) (poly* f g))

;;; exponentiation
(defmethod poly-expt ((lst list) (n fixnum))
  (cond ((< n 0) (error "not implemented"))
	((= n 0) (list 1))
	((= n 1) lst)
	((evenp n) (poly-expt (poly* lst lst) (/ n 2)))
	(t (poly* lst (poly-expt lst (1- n))))))

(defmethod poly-expt ((poly <polynomial>) (n fixnum))
  (make-polynomial (poly-expt (poly-coeffs poly) n)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; matlisp methods
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmethod m+ ((x <polynomial>) (y <polynomial>)) (vec+ x y))
(defmethod m+ ((x number) (y <polynomial>)) (vec+ x y))
(defmethod m+ ((x <polynomial>) (y number)) (vec+ x y))

(defmethod m- ((x <polynomial>) (y <polynomial>)) (vec- x y))
(defmethod m- ((x number) (y <polynomial>)) (vec- x y))
(defmethod m- ((x <polynomial>) (y number)) (vec- x y))

(defmethod scal ((val number) (poly <polynomial>))
  (make-polynomial
   (vec-s* val (poly-coeffs poly))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; differentiation
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; differentiates an n-variate polynomial wrt the variable given by
;;; index=0,1,...
(defmethod differentiate ((lst list) (index fixnum))
  (cond ((null lst) ())
	((zerop index)
	 (loop for coeff in (cdr lst)
	       and deg from 1
	       collect (vec-s* deg coeff)))
	(t (mapcar #'(lambda (coeff)
		       (differentiate coeff (- index 1)))
		   lst))))

(defmethod differentiate ((poly number) (index fixnum)) 0)

(defmethod differentiate ((poly <polynomial>) (index fixnum))
  (make-polynomial (differentiate (poly-coeffs poly) index)))

;;; warning: does not work for multidimensional polynomials (e.g. p(x1,x2)=x1)
(defmethod gradient ((poly <polynomial>))
  (loop for index from 0 below (variance poly)
	collect (differentiate poly index)))

(defmethod evaluate-gradient ((poly <polynomial>) (x array))
  (let ((result (make-array (length x))))
    (dotimes (index (length result) result)
      (setf (aref result index)
	    (evaluate (differentiate poly index) x)))))

(defmethod k-jet ((poly <polynomial>) (k integer) (dim integer))
  (loop for order upto k
	for Di = poly then
	(let ((Dnext (make-general-tensor (make-fixnum-vec order dim))))
	  (dotimes (index dim Dnext)
	    (if (= order 1)
		(setf (tensor-ref Dnext index)
		      (differentiate poly index))
		(copy! (tensor-map t (rcurry #'differentiate index) Di)
		       (slice Dnext (list (cons (1- order) index)))))))
	collecting Di))

(defmethod evaluate-k-jet ((poly <polynomial>) (k integer) (x array))
  (loop with k-jet = (k-jet poly k (length x))
	for i from 0
	and Di in k-jet collecting
	(if (zerop i)
	    (evaluate poly x)
	    (tensor-map 'double-float (rcurry #'evaluate x) Di))))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Diversities
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun exponents->monomial-list (partition)
  "Converts a monomial given by a list of powers of its components into our
dense polynomial format.

Example: (exponents->monomial-list '(1 2)) -> (0 (0 0 1))
"
  (if (null partition)
      '(1)                              ; special case
      (append (make-list (car partition) :initial-element 0)
              (if (single? partition)
                  (list 1)
                  (list (exponents->monomial-list (cdr partition)))))))

;;; Example: (exponents->monomial-list '(1 2)) -> (0 (0 0 (1)))
;(defun exponents->monomial-list (partition)
;  (if (null partition)
;      (list 1)                          ; special case
;      (append (make-list (car partition) :initial-element 0)
;              (list (exponents->monomial-list (cdr partition))))))

;;; 
(defun n-variate-monomials-of-degree (n degree &optional (type '=))
  "Returns n-variate monomials of degree being equal or being lower or
equal than deg.
Examples:
 (n-variate-monomials-of-degree 2 2) -> (x2^2 x1*x2 x1^2)
 (n-variate-monomials-of-degree 2 2 '<=)  -> (1 x2 x1 x2^2 x1*x2 x1^2)"
  (ecase type
    (= (loop for partition in (n-partitions-of-k n degree)
	     collect (make-polynomial (exponents->monomial-list partition))))
    (<= (loop for deg from 0 upto degree
	      nconcing (n-variate-monomials-of-degree n deg)))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Interpolation polynomials
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun lagrange-polynomials (points)
  (loop for xi in points
	for points-without-xi = (remove xi points)
	collect
	(vec-s*
	 (/ 1 (reduce #'* (mapcar #'(lambda (xj) (- xi xj)) points-without-xi)))
	 (reduce #'poly* (mapcar #'(lambda (xj)
				     (make-polynomial (list (- xj) 1)))
				 points-without-xi)
		 :initial-value (make-polynomial '(1))))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Simple 1d integration
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun integrate-simple-polynomial (poly)
  (make-polynomial
   (cons 0 (if (numberp poly)
	       (list poly)
	       (loop for coeff in (poly-coeffs poly)
		     and n+1 from 1
		     collect (/ coeff n+1))))))


;;; Testing
(defun test-polynom ()
  (integrate-simple-polynomial (make-polynomial '(0 1)))
  (poly-coeffs (make-<polynomial>))
  (shift-polynomial (make-polynomial '(0 1)) 1)
  (let ((p (make-polynomial '((1 2) (1 2))))
	(p-x1 (make-polynomial '(0 1)))
	(p-x2 (make-polynomial '((0 1))))
	(p-1-x3 (make-polynomial '(((1 -1))))))
    (poly* p-x1 p-x2)
    (gradient p)
    (differentiate p-x1 0)
    (gradient p-1-x3)
    (evaluate-gradient p-1-x3 #(1.0 2.0))
    (make-polynomial '(4))
    (k-jet p 2 2)
    (evaluate-k-jet (make-polynomial '(1 2 ((1 2 3)))) 2 #(0.0d0 0.0d0))
    ))

;;; (test-polynom)
(tests::adjoin-femlisp-test 'test-polynom)

