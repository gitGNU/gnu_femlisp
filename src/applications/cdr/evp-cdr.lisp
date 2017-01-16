;;; -*- mode: lisp; fill-column: 64; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; evp-cdr.lisp - computing eigenvalues for cdr problems
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

(in-package :fl.application)

(file-documentation "This file contains routines for computing
eigenvalue/eigenvector pairs for convection-diffusion-reaction problems.")

(defun laplace-eigenvalue-computation
    (domain &key output plot (dirichlet 0.0) (multiplicity 1)
     (shift (if dirichlet 0.0 1.0)) (time 20) (nr-levels 5) (base-level 1))
  "Function performing the eigenvalue demo for the Laplace operator."
  (let ((problem (cdr-model-problem
		  domain :evp t :multiplicity multiplicity
                  :dirichlet dirichlet)))
    (storing
      (solve (blackboard
	      :problem problem :base-level base-level
	      :success-if `(or (>= :time ,time) (>= :nr-levels ,nr-levels))
	      :output output :observe
	      (append *stationary-fe-strategy-observe*
		      (list (list "             lambda" "~A"
				  (_ (map 'vector (rcurry #'- shift)
                                          (slot-value problem 'eigenvalues)))))))))
    (when plot
      (loop for i below multiplicity do
           (plot (getbb *result* :solution) :index i)))))

(defun make-laplace-eigenvalue-demo (domain domain-name)
  (let ((title domain-name)
	(short (format nil "Eigenvalues of Laplace on a ~A." domain-name))
	(long (format nil "Computes eigenvalues for the Laplace
operator on a ~A.  The solution strategy does uniform refinement
and terminates if more than 20 seconds have passed after a
step." domain-name)))
    (let*  ((dim (dimension domain))
            (demo
              (make-demo
               :name title :short short :long long
               :execute (lambda ()
                          (laplace-eigenvalue-computation
                           domain
                           :multiplicity 3 :base-level (if (= dim 1) 1 0)
                           :plot t :output 1)))))
      (adjoin-demo demo *eigenvalue-demo*))))

#+(or)
(let ((FL.DISCRETIZATION:*SUGGESTED-DISCRETIZATION-ORDER* 1))
  (laplace-eigenvalue-computation
   (?2 (cylinder-domain 1.0 1.0)
       (circle-ring-domain 1.0 1.5 :channel-breadth 0.0)
       (reservoir-domain 1.0 1.5 1.0))
   :multiplicity 1 :time 150 :nr-levels 5 :output 1))

;; (laplace-eigenvalue-computation (cylinder-domain 1.0 1.0)
;;                                 :multiplicity 3 :time 1000 :nr-levels 3 :output :all)
;;; First zero of Bessel(0) = lambda_1 = 2.404825557695773
;;; Therefore eigenvalue = lambda_1^2 = 5.78318596297

;;;                                        5.783183821630505
;;;                                        5.783185865288665
;;; cylinder
;;; lambda1^2+pi^2= 15.6527903641

;; 15.338358832149563
;; 15.641985755747394

;;; auf großem Gebiet:
;;; lambda1/4=1.445796490588467
;;; 1.4457860789623667
;;; 1.4457964410120385
;;; 1.4457964903780134
;;;       1:  1.44579649074     scheint auch sehr gut zu stimmen
;;;       2:  1.4457964906571028

;;; Mit Dreieck
;;; 1.4461086803627932
;;; 1.4461203521690598
;;; 1.4461201219288802
;;; 1.4461199921531585
;;; -> sieht auch gut aus

;;; Mit Neumann-Randbedingungen

#+(or)
(laplace-eigenvalue-computation (circle-ring-domain 1.0 1.5)
                                :multiplicity 3 :time 1000 :nr-levels 3 :output 1)

;;; Order=1 auf Kreis(ring) mit r=1.5 (exakt=5.78318596297/2.25=2.57030487243)
   ;; 0       8        27     0.4  #(2.2694504532190947 3.262817890161257
   ;;                                3.83328738108937)  
   ;; 1      32       114     1.3  #(2.5845092366236324 7.08947513627384
   ;;                                7.533151431978448)  
   ;; 2     128       429    11.9  #(2.575413358209864 6.765188827135567
   ;;                                6.766696316748804)  
   ;; 3     512      1632    48.9  #(2.5716794784553527 6.604136808837907
   ;;                                6.60453520088074)  
   ;; 4    2048      6339   199.7  #(2.5706554413092086 6.5496691543420225
   ;;                                6.549770322754479)  

;;;(plot (getbb *result* :solution) :index 0)
(make-laplace-eigenvalue-demo (n-simplex-domain 1) "unit-interval")
(make-laplace-eigenvalue-demo (n-cube-domain 2) "unit-quadrangle")
(make-laplace-eigenvalue-demo (circle-ring-domain 1.0 nil) "unit-circle")
(make-laplace-eigenvalue-demo (cylinder-domain 1.0 1.0) "cylinder")

#+(or)
(let ((problem (cdr-model-problem
		(n-cube-domain 1)
		:multiplicity 2
		:reaction #m(1.0) :source nil :dirichlet nil :evp t)))
  (storing
    (solve (blackboard
	    :problem problem :base-level 2
	    :success-if '(>= :nr-levels 4)
	    :output :all :observe
	    (append *stationary-fe-strategy-observe*
		    (list (list "             lambda" "~A"
				#'(lambda (bb)
				    (declare (ignore bb))
				    (slot-value problem 'eigenvalues)))))))))

;;;; Testing

(defun evp-cdr-test ()
  (plot (laplace-eigenvalue-computation
	 (n-cube-domain 1)
         :multiplicity 3 :base-level 1
         :output t :plot t))

  ;; the following is used in the manual
  (let ((problem (cdr-model-problem 2 :evp t)))
    (storing
      (solve (blackboard :problem problem
			 :success-if '(or (>= :time 5) (>= :nr-levels 5))
			 :output 1))))
  (slot-value (getbb *result* :problem) 'eigenvalues)
  (plot (getbb *result* :solution))
  
  ;; bug shows up when using Lispworks
  (let ((A #m((8.0 -12.0) (-12.0 18.0)))
        (B (eye 2)))
    (hegv A B :V))
  
  )

;;; (evp-cdr-test)
(fl.tests:adjoin-test 'evp-cdr-test)
