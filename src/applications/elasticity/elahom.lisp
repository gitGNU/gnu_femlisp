;;; -*- mode: lisp; fill-column: 64; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; elahom.lisp - homogenized coefficients for linear elasticity
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

(defun elasticity-cell-problem-gamma (dim)
  "Returns a right-hand side for an elasticity cell problem."
  (ensure-coefficient
   'FL.ELLSYS::H
   (let* ((gamma (make-array dim :initial-element nil)))
     (dotimes (j dim)
       (setf (aref gamma j)
	     (let ((entry (zeros dim (* dim dim))))
	       (dotimes (mu dim)
		 (dotimes (k dim)
		   (dotimes (nu dim)
		     (and (= mu nu) (= k j)
			  (setf (mref entry mu (+ nu (* dim k)))
				1.0)))))
	       entry)))
     gamma)))

(defun elasticity-cell-problem-reaction (dim value)
  "We want to have zero at the corners."
  (diagonal-reaction-coefficient (make-double-vec dim value)))

(defun elasticity-inlay-cell-problem (domain &key (inlay-p #'patch-in-inlay-p)
				      (interior 100.0))
  "Generates the inlay cell problem.  The coefficient is of the
form @math{A^{ij}_{\lambda \mu}}$, and the equation to be solved
is in variational form

@math{ 0 = \partial_\lambda \phi^i A^{ij}_{\lambda \mu}
    (\partial_\mu N^j_m + \Gamma^j_{\mu m}) }

With @math{k,\nu} being such that @math{m= k n + \nu} we have
@math{\Gamma^{j_k}{\mu\nu}=\delta_{jk}\delta_{\mu\nu}} leading
to matrix-valued functions @math{N_\nu} in the notation of
\cite{Bakhvalov-Panasenko}.  We should also have
@math{N^{jk}_\nu = N^{j\nu}_k} because of the symmetries of A.
The correction to the arithmetic average has to be computed as

@math{ A_{corr}^{kr}_{iq} := \int_Y A^{kl}_{ij} \grad_{x_j}
N^{lr}_q = \int_Y \div_{x_j} A^{lk}_{ji} N^{lr}_q = F[k*dim+i]
\cdot N[r*dim+q].}"
  (let ((dim (dimension domain)))
    (create-problem '<elasticity-problem>
	(:domain domain :multiplicity (expt dim 2))
      (setup-coefficients (patch)
	(let ((eps (if (funcall inlay-p patch) interior 1.0)))
	  (select-on-patch (patch)
	    (:d-dimensional
	     (list
	      (isotropic-elasticity-tensor-coefficient dim eps)
	      (elasticity-cell-problem-gamma dim)))
            ((and :0-dimensional :boundary (not :inlay))
             (list
              (elasticity-cell-problem-reaction dim 10.0)))))))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; utilities
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun inlay-block-decomposition (blockit asa)
  "Block decomposition for obtaining a coarse-grid solver which is robust
against the size of the coefficient jump."
  (declare (ignore blockit))
  (let ((inlay-block ())
	(rest-blocks ())
	(mesh (mesh asa)))
    (for-each-row-key
     #'(lambda (key)
	 (let ((cell (representative key)))
	   (if (patch-in-inlay-p (patch-of-cell cell mesh))
	       (push key inlay-block)
	       (push (list key) rest-blocks))))
     asa)
    (mapcar (rcurry #'coerce 'vector) (cons inlay-block rest-blocks))))
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Demos
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun elasticity-interior-effective-coeff-demo
    (problem &key order levels plot (output 1) effective-tensor-p)
  "Computes the effective elasticity for a certain periodic
medium.  The approximation is done with finite elements and
blending.  Uniform refinement is used, the linear solver is a
geometric multigrid scheme used in a nested iteration fashion.
A subtle point here is that the coarse-grid smoother has to
treat the inlay as a block to be robust against the large
coefficient jump.

The solution to this cell problem is a tensor field of rank 3
with dim^3 components which are plotted one after the other."
  (let* ((domain (domain problem))
	 (dim (dimension domain))
	 (*output-depth* output)
         saved-effective-tensor)
    (storing
      (solve
       (blackboard
	:problem problem
	:fe-class (lagrange-fe order :nr-comps dim :type :gauss-lobatto)
	:estimator (make-instance '<projection-error-estimator>)
	:indicator (make-instance '<largest-eta-indicator> :fraction 1.0)
	:success-if `(>= :max-level ,(1- levels))
	:solver
        (?1
         (let* ((smoother (make-instance '<jacobi> :damp 1.0))
                (cs (geometric-cs
                     :gamma 1 :smoother smoother :pre-steps 1 :post-steps 0
                     :coarse-grid-iteration
                     (make-instance '<multi-iteration> :base smoother :nr-steps 1)
                     :combination :additive))
                (bpx (make-instance '<cg> :preconditioner cs :restart-cycle 30)))
           (make-instance '<linear-solver>
                          :iteration bpx
                          :success-if `(and (> :step 2) (> :step-reduction 1.0) (< :defnorm 1.0e-8))
                          :failure-if `(and (> :step 200) (> :step-reduction 1.0) (> :defnorm 1.0e-8))))
         (make-instance
          '<linear-solver> :iteration
          (let ((smoother
                  (?2 (make-instance '<cg> :preconditioner (make-instance '<jacobi>))
                      (if (>= dim 3)
                          *gauss-seidel*
                          (geometric-ssc)))))
            (geometric-cs
             :coarse-grid-iteration
             (make-instance '<multi-iteration> :nr-steps 10 :base
                            (?1 *gauss-seidel*
                                #+(or)
                                (make-instance '<custom-ssc>
                                               :block-setup #'inlay-block-decomposition)))
             :smoother smoother :pre-steps 2 :post-steps 2
             :gamma 2 :fmg t))
          :success-if `(and (> :step 2) (> :step-reduction 0.9) (< :defnorm 1.0e-9))
          :failure-if `(and (> :step 100) (> :step-reduction 0.9) (> :defnorm 1.0e-9))
          ))
	:plot-mesh plot
	:observe
	(append *stationary-fe-strategy-observe*
                (list fl.strategy::*mentries-observe*)
		(list
		 (list (format nil "~19@A~19@A~19@A" "A^00_00" "A^00_11" "A^10_10") "~57A"
		       #'(lambda (blackboard)
			   (let ((tensor (effective-tensor blackboard)))
                             (setq saved-effective-tensor tensor)
                             (format nil "~19,10,2E~19,10,2E~19,10,2E"
                                     (and tensor (mref (mref tensor 0 0) 0 0))
                                     (and tensor (> dim 1) (mref (mref tensor 0 0) 1 1))
                                     (and tensor (> dim 1) (mref (mref tensor 1 0) 1 0)))))))))
		   ))
    (when plot
      ;; plot components of cell solution tensor
      (dotimes (j dim)
	(dotimes (k dim)
	  (plot (getbb *result* :solution) :index (+ (* j dim) k))
	  (sleep 1.0))))
    ;; compute the homogenized coefficient
    (awhen (and effective-tensor-p
                (or saved-effective-tensor (effective-tensor *result*)))
      (format t "The effective elasticity tensor is:~%~A~%" it)
      it)))

#+(or)
(time (elasticity-interior-effective-coeff-demo
       (elasticity-inlay-cell-problem (n-cell-with-ball-inlay 1))
       :order 1 :levels 2 :plot nil :output 2))

#|
(prof:with-profiling (:type :time)
  (elasticity-interior-effective-coeff-demo
   (elasticity-inlay-cell-problem (n-cell-with-ball-inlay 2))
   :order 4 :levels 1 :output :all))
(prof:show-call-graph)
|#

#+(or)
(plot (getbb *result* :solution) :index 1 :depth 2)

#+(or)
(let ((dim 2)
      (counter -1))
  (dotimes (i dim)
    (dotimes (j dim)
      (dotimes (k dim)
	(plot (getbb *result* :solution) :component i :index (+ (* j dim) k)
	      :background :white
	      :filename (format nil "ela-cell-sol-x~D" (incf counter)))))))
;; montage -geometry 480x480 -tile 2x4 ela-cell-*.tiff ela-cell-sols.eps

(defparameter *effective-elasticity-demo*
  (make-demo
   :name "effective-elasticity"
   :short "Computes an effective elasticity tensor"
   :long (documentation 'elasticity-interior-effective-coeff-demo 'function)))

(adjoin-demo *effective-elasticity-demo* *interior-coeffs-demo*)
(adjoin-demo *effective-elasticity-demo* *elasticity-demo*)


(defun make-effective-elasticity-porous-domain-demo (dim order levels)
  (let* ((short (format nil "Effective tensor for a ~DD cell with hole" dim))
	 (demo
	  (make-demo
	   :name (format nil "~DD-hole-cell" dim)
	   :short short
	   :long (format nil "We compute the effective
elasticity tensor for a medium with periodically distributed
holes by solving a cell problem on a periodicity cell.~%
Parameters: order=~D, max-levels=~D~%~%"
			      order levels)
		:execute
		(lambda ()
		  (elasticity-interior-effective-coeff-demo
		   (elasticity-inlay-cell-problem (n-cell-with-ball-hole dim))
		   :order order :levels levels :plot t)))))
    (adjoin-demo demo *effective-elasticity-demo*)))

(make-effective-elasticity-porous-domain-demo 2 5 2)
(make-effective-elasticity-porous-domain-demo 3 4 1)

(defun make-effective-elasticity-inlay-domain-demo (dim order levels)
  (let* ((short (format nil "Effective tensor for a ~DD cell with inlay." dim))
	 (demo
	  (make-demo
	   :name (format nil "~DD-inlay-cell" dim) :short short
	   :long
	   (format nil "We compute the effective elasticity
tensor for a mixture of two isotropic media forming a ball
inlay.  The Lame constants are 1/1 outside the inlay and 100/100
inside.
Parameters: order=~D, levels=~D~%~%"
		   order levels)
		:execute
		(lambda ()
		  (elasticity-interior-effective-coeff-demo
		   (elasticity-inlay-cell-problem (n-cell-with-ball-inlay dim))
		   :order order :levels levels)))))
    (adjoin-demo demo *effective-elasticity-demo*)))

(make-effective-elasticity-inlay-domain-demo 2 4 2)
(make-effective-elasticity-inlay-domain-demo 3 4 2)

#+(or)(demo *effective-elasticity-demo*)

(defun elahom-performance-calculation (&key (levels 2) (output 2) (threshold 2e-10))
  (elasticity-interior-effective-coeff-demo
   (elasticity-inlay-cell-problem (n-cell-with-ball-hole 3))
   :order 5 :levels levels :plot nil :output output)
  (with-items (&key defnorm initial-defnorm)
      (getbb *result* :solver-blackboard)
    (format t "Initial defect norm: ~A  --> Defect norm at the end: ~G~%"
            initial-defnorm defnorm)
    (when (member levels '(1 2))
      (let ((serial-defnorm (ecase levels
                              (1 1.56e-09)
                              (2 3.13e-09))))
        (assert (< (abs (- defnorm serial-defnorm)) threshold))))))

(defun elahom-performance-test (nrs-of-kernels &rest args &key initialize &allow-other-keys)
  (when initialize
    (format t "We perform one computation on two levels for initializing FE and interpolation:~%")
    (elahom-performance-calculation :levels 2))
  ;; then we test performance for varying number of kernels
  (loop for n in nrs-of-kernels do
    (sb-ext:gc :full t)
    (new-kernel n)
    (format t "Testing for ~D kernels:~%" n)
    (apply #'elahom-performance-calculation (sans args :initialize))))

;; (elahom-performance-test '(1 2) :levels 3)
(defun elahom-longtime-stability-test (nr-of-trials &rest args
                                       &key (levels 1) (threshold 1e-10)
                                       &allow-other-keys)
  (loop for i below nr-of-trials do
    (format t "~&~%*** Rechnung ~D~%~%" i)
       (apply #'elahom-performance-calculation
              :levels levels :threshold threshold
              (sans args :levels :threshold))))

;;; (elahom-performance-calculation :levels 2)
;;; (elahom-longtime-stability-test 2 :levels 1)

(defun mconvert (mat)
  (flet ((key-converter (key)
           (if (typep key 'fl.mesh::identification)
               (mapcar #'midpoint (fl.mesh::cells key))
               (midpoint key))))
    (map-matrix (make-instance 'fl.matlisp::<ht-sparse-matrix> :test 'equalp)
                (lambda (e i j)
                  (values e (key-converter i) (key-converter j)))
                mat)))

;; (setq *mymat* (mconvert *mymat*))
;; (row-keys *mymat*)
;; (x<-0 *mymat*)
;; (loop for i below 100 do
;;   (x<-0 *mymat*)
;;   (assemble-interior (ansatz-space *mymat*) :surface :matrix *mymat*)
;;   (assert (< (norm (m-! (mconvert *mymat0*) (mconvert *mymat*))) 1e-10)))



(defun test-homogenization-elasticity ()

  (time (elasticity-interior-effective-coeff-demo
	 (elasticity-inlay-cell-problem (n-cell-with-ball-inlay 2)) :order 3 :levels 2
         :output :all))
  ;;(profile:report-time)
  ;;(profile:profile :methods 'fl.matlisp:sparse-matrix->matlisp)
  ;;(profile:unprofile)
  ;; we should have N^{lr}_q = N^{lq}_r
  (with-items (&key mesh solution ansatz-space rhs) *result*
    (let ((average-tensor
	   (average-coefficient ansatz-space :coefficient 'FL.ELLSYS::A))
	  (correction-tensor
	   (convert-elasticity-correction
	    (correction-tensor solution rhs)))
	  (dim (dimension mesh)))
      (check-elasticity-tensor average-tensor dim 0.0)
      (check-elasticity-tensor correction-tensor dim 1.0e-4)
      (check-elasticity-tensor (m- average-tensor correction-tensor) dim)))
  (effective-tensor *result*)

  (isotropic-elasticity-tensor :dim 2 :lambda 100 :mu 100)

  (plot (getbb *result* :solution) :component 0 :index 1)

  #+(or)
  (sb-sprof:with-profiling
      (:max-samples 10000)
    (time
     (with-items (&key rhs matrix) *result*
       (linsolve matrix rhs
                 :output t
                 :threshold 1e-10
                 :iteration
                 (make-instance '<multi-iteration> :nr-steps 5 :base
                                (make-instance '<cg> :preconditioner (make-instance '<jacobi>)))))))
  
  ;; Testing for bug
  #+(or)
  (progn
    (sb-ext:gc :full t)
    (sb-sprof:reset)
    (sb-sprof:start-profiling)
    (list sb-sprof:*sample-interval* sb-sprof:*alloc-interval* sb-sprof:*max-samples*)
    (lret ((fl.matlisp::*blas3-operation-count* 0))
      (time
       (elasticity-interior-effective-coeff-demo
        (elasticity-inlay-cell-problem (n-cell-with-ball-hole 3))
        :order 5 :levels 2 :plot nil :output 2)))
    (sb-sprof:stop-profiling))
  #+(or)
  (sb-sprof:report :type :flat :sort-by :cumulative-samples :max 200)

  ;; testing for bug
  (loop repeat 3
        for old = nil then tensor
        and tensor = (elasticity-interior-effective-coeff-demo
                      (elasticity-inlay-cell-problem (n-cell-with-ball-hole 3))
                      :order 5 :levels 1 :plot nil :output 2)
        do
           (when old (assert (mequalp old tensor))))

  (time
   (flet ((computation ()
            (elasticity-interior-effective-coeff-demo
             (elasticity-inlay-cell-problem (n-cell-with-ball-hole 3))
             :order 5 :levels 2 :plot nil :output 1)))
     (with-workers (#'computation)
       (loop repeat 2 do (work-on)))))
  (time
   (elasticity-interior-effective-coeff-demo
    (elasticity-inlay-cell-problem (n-cell-with-ball-hole 3))
    :order 5 :levels 2 :plot nil :output 2))
  
  (assemble-interior (getbb *result* :ansatz-space) :surface)
  
  (loop for i below 1 do
    (format t "~&~%*** Rechnung ~D~%~%" i)
    (elasticity-interior-effective-coeff-demo
     (elasticity-inlay-cell-problem (n-cell-with-ball-hole 3))
     :order 5 :levels 1 :plot nil :output 2)
        do
           (let ((defnorm (getbb (getbb *result* :solver-blackboard) :defnorm)))
             (assert (< (abs (- defnorm 1.56164e-09)) 1e-5))
             (format t "~G~%" defnorm)
             defnorm))
  
;; >       0   3.40515e-01     UNDEFINED     
;; >       1   4.68740e-01      1.37656e+00  
;; >       2   3.52337e-01      7.51667e-01  
;; >       3   3.48071e-01      9.87893e-01  
;; >       4   2.84252e-01      8.16649e-01  
;; >       5   2.37511e-01      8.35566e-01  
;; >       6   1.85954e-01      7.82927e-01  
;; >       7   1.76617e-01      9.49789e-01  
;; >       8   1.40476e-01      7.95370e-01  
;; >       9   1.33209e-01      9.48273e-01  
;; >      10   1.00517e-01      7.54576e-01  
;; >      11   8.32486e-02      8.28207e-01  
;; >      12   6.90603e-02      8.29567e-01  
;; >      13   6.70766e-02      9.71276e-01  
;; >      14   5.59697e-02      8.34415e-01  
;; >      15   4.41744e-02      7.89255e-01  
;; >      16   3.32768e-02      7.53305e-01  
;; >      17   2.51508e-02      7.55806e-01  
;; >      18   2.08140e-02      8.27569e-01  
;; >      19   1.59011e-02      7.63960e-01  
;; >      20   1.31436e-02      8.26587e-01  
;; >      21   9.56384e-03      7.27642e-01  
;; >      22   6.48499e-03      6.78074e-01  
;; >      23   5.64712e-03      8.70799e-01  
;; >      24   4.30671e-03      7.62638e-01  
;; >      25   2.66636e-03      6.19119e-01  
;; >      26   2.29514e-03      8.60774e-01  
;; >      27   1.69797e-03      7.39814e-01  
;; >      28   1.10722e-03      6.52084e-01  
;; >      29   7.31214e-04      6.60404e-01  
;; >      30   4.86785e-04      6.65722e-01  
;; >      31   3.04622e-04      6.25784e-01  
;; >      32   2.00025e-04      6.56633e-01  
;; >      33   1.39029e-04      6.95059e-01  
;; >      34   1.02067e-04      7.34142e-01  
;; >      35   7.15768e-05      7.01270e-01  
;; >      36   5.85013e-05      8.17323e-01  
;; >      37   4.51134e-05      7.71152e-01  
;; >      38   2.76396e-05      6.12670e-01  
;; >      39   1.89543e-05      6.85766e-01  
;; >      40   1.39972e-05      7.38468e-01  
;; >      41   8.84626e-06      6.32004e-01  
;; >      42   6.74961e-06      7.62991e-01  
;; >      43   5.01588e-06      7.43137e-01  
;; >      44   3.78946e-06      7.55492e-01  
;; >      45   2.72662e-06      7.19527e-01  
;; >      46   2.05114e-06      7.52264e-01  
;; >      47   1.61976e-06      7.89689e-01  
;; >      48   1.27186e-06      7.85212e-01  
;; >      49   1.02466e-06      8.05645e-01  
;; >      50   9.33554e-07      9.11082e-01  
;; >      51   7.28848e-07      7.80724e-01  
;; >      52   6.36419e-07      8.73185e-01  
;; >      53   4.85990e-07      7.63632e-01  
;; >      54   4.26887e-07      8.78386e-01  
;; >      55   3.59166e-07      8.41360e-01  
;; >      56   3.03773e-07      8.45774e-01  
;; >      57   2.86535e-07      9.43252e-01  
;; >      58   2.17011e-07      7.57364e-01  
;; >      59   1.60925e-07      7.41554e-01  
;; >      60   1.18693e-07      7.37569e-01  
;; >      61   9.81600e-08      8.27004e-01  
;; >      62   1.04974e-07      1.06942e+00  
;; >      63   9.50038e-08      9.05022e-01  
;; >      64   7.03516e-08      7.40513e-01  
;; >      65   6.18515e-08      8.79177e-01  
;; >      66   5.18942e-08      8.39013e-01  
;; >      67   3.56091e-08      6.86187e-01  
;; >      68   2.52713e-08      7.09686e-01  
;; >      69   1.95321e-08      7.72898e-01  
;; >      70   1.52263e-08      7.79551e-01  
;; >      71   1.32333e-08      8.69106e-01  
;; >      72   1.16680e-08      8.81719e-01  
;; >      73   1.11058e-08      9.51819e-01  
;; >      74   9.34913e-09      8.41821e-01  
;; >      75   8.35822e-09      8.94010e-01  
;; >      76   6.47601e-09      7.74807e-01  
;; >      77   5.06761e-09      7.82520e-01  
;; >      78   3.86510e-09      7.62706e-01  
;; >      79   3.06596e-09      7.93243e-01  
;; >      80   2.63208e-09      8.58484e-01  
;; >      81   2.29936e-09      8.73591e-01  
;; >      82   2.11525e-09      9.19929e-01  
;; >      83   1.69522e-09      8.01430e-01  
;; >      84   1.55839e-09      9.19283e-01  
  ;; >      85   1.56164e-09      1.00209e+00


  ;; Bug: mit 2 Workern, Diskretisierung seriell
  ;; >      58   2.17011e-07      7.57364e-01  
  ;; >      59   3.26301e-07      1.50362e+00  
  ;; >      60   2.36785e-07      7.25665e-01  

  )

;;; (test-homogenization-elasticity)
(fl.tests:adjoin-test 'test-homogenization-elasticity)
