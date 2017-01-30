;;; -*- mode: lisp; fill-column: 64; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; heisig-neuss-2017.lisp - Calculations for [MHeisig-NNeuss-2017]-Part I
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;
;;; Copyright (C) 2017
;;; Nicolas Neuss, Friedrich-Alexander-Universitaet Erlangen-Nuernberg
;;; All rights reserved.
;;; 
;;; Redistribution and use in source and binary forms, with or without
;;; modification, are permitted provided that the following conditions
;;; are met:
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
;;; MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
;;; IN NO EVENT SHALL THE AUTHOR, THE FAU ERLANGEN-NUERNBERG, OR OTHER
;;; CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
;;; EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
;;; PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
;;; PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
;;; LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
;;; NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
;;; SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(in-package :cl-user)

(defpackage "FL.KONWIHR"
  (:use "COMMON-LISP"
	"FL.MACROS" "FL.UTILITIES" "FL.MATLISP"
	"FL.DEBUG" "FL.TESTS" "FL.DEMO" "FL.PARALLEL"
        "FL.DICTIONARY"        
	"FL.FUNCTION" "FL.MESH"
	"FL.PROBLEM" "FL.CDR" "FL.ELLSYS" "FL.ELASTICITY" "FL.NAVIER-STOKES"
	"FL.DISCRETIZATION" "FL.ELLSYS-FE" "FL.ELASTICITY-FE" "FL.NAVIER-STOKES-FE"
	"FL.ITERATION" "FL.MULTIGRID" "FL.GEOMG"
	"FL.STRATEGY" "FL.PLOT"
	"FL.DOMAINS" "FL.APPLICATION")
  (:import-from "FL.APPLICATION"
                "ELASTICITY-INTERIOR-EFFECTIVE-COEFF-DEMO"
                "ELASTICITY-INLAY-CELL-PROBLEM"
                "*ARTICLES-DEMO*"
                )
  (:export "KONWIHR-PAPER-MAX-LEVELS" "KONWIHR-SPEED"
           "*KONWIHR-INITIALIZED*" "INITIALIZE-KONWIHR-PAPER-CALCULATION"
           "ELASTICITY-INTERIOR-EFFECTIVE-COEFF-DEMO"
           "ELASTICITY-INLAY-CELL-PROBLEM"
           "HEISIG-NEUSS-2017-DEMO" "*HEISIG-NEUSS-2017-DEMO*"
           "*ARTICLES-DEMO*"
           )
  (:documentation "Special package for the KONWIHR demo problem."))

(in-package :fl.konwihr)

(file-documentation
 "Demos for the calculations reported in the KONWIHR paper by M. Heisig and N. Neuss.")

(defun konwihr-paper-max-levels ()
  (whereas ((mem (fl.port:dynamic-space-size)))
    (floor (+ 1.9 (log (/ mem 1e9) 8)))))

(defvar *konwihr-speed* nil
  "For memoizing the measured CL speed (for having a nicer demo experience:-)")

(defun konwihr-speed ()
  "Returns estimated speed of the current Femlisp in MLOPS."
  (unless *konwihr-speed*
    (format t "~&Estimating CPU speed. This may take a little time, so please wait...")
    (force-output)
    (setq *konwihr-speed* (common-lisp-speed))
    (format t " Measured ~D MFLOPS.~%" (* 10 (round *konwihr-speed* 10)))
    (force-output))
  *konwihr-speed*)

(defun initialize-elahom-calculation (dim order &optional (levels 2))
  "Initializes all local finite element data and local
interpolation matrices for dimension @arg{dim} and finite
elements of order @arg{order} on cubic cells.  This file should
be loaded before starting elahom performance measurements."
  (let* ((domain (n-cube-domain dim))
         (problem (elasticity-model-problem
                   domain :lambda 1.0 :mu 1.0
                   :force (ellsys-one-force-coefficient dim 1)
                   :dirichlet nil)))
    (storing
      (solve
       (blackboard
	:problem problem
	:fe-class (lagrange-fe order :nr-comps dim :type :gauss-lobatto)
	:estimator (make-instance '<projection-error-estimator>)
	:indicator (make-instance '<largest-eta-indicator> :fraction 1.0)
	:success-if `(>= :max-level ,(1- levels))
	:solver
        (let* ((smoother (make-instance '<jacobi> :damp 1.0))
               (cs (geometric-cs
                    :gamma 1 :smoother smoother :pre-steps 1 :post-steps 0
                    :coarse-grid-iteration
                    (make-instance '<multi-iteration> :base smoother :nr-steps 1)
                    :combination :additive))
               (bpx (make-instance '<cg> :preconditioner cs :restart-cycle 30)))
          (make-instance '<linear-solver>
                         :iteration bpx
                         :success-if `(> :step 2)
                         :failure-if `(> :step 100)))
	:plot-mesh nil
	:observe *stationary-fe-strategy-observe*
        :output 0)))))

(defvar *konwihr-initialized* nil
  "A flag keeping track, if the setup for these tests has already been performed.")

(defun initialize-konwihr-paper-calculation ()
  (unless *konwihr-initialized*
    (format t "We initialize FEs and interpolation matrices for
the KONWIHR calculations.

This calculation may take about ~D seconds, but it is necessary only once
for all these demos, so please wait..."
            (round (/ 20000 (konwihr-speed))))
    (force-output)
    (initialize-elahom-calculation 3 5 2)
    (setq *konwihr-initialized* t)
    (format t "  Done.~%~%")
    (force-output)))

(defun heisig-neuss-2017-demo
    (problem &key order levels (output 1) distributed-p
               (initial-mesh-refinements 0))
  "This is a slightly trimmed copy of the function
@function{elasticity-interior-effective-coeff-demo}."
  (let* ((domain (domain problem))
	 (dim (dimension domain))
	 (*output-depth* output)
         saved-effective-tensor)
    (storing
      (solve
       (blackboard
	:problem problem
        :initial-mesh-refinements initial-mesh-refinements
	:fe-class (lagrange-fe order :nr-comps dim :type :gauss-lobatto)
	:estimator nil
	:indicator (make-instance '<uniform-refinement-indicator>)
	:success-if `(>= :max-level ,(1- levels))
	:solver
        (let* ((smoother (make-instance
                          (if distributed-p
                              (intern "<DISTRIBUTED-JACOBI>" (find-package :femlisp-ddo))
                              '<JACOBI>)
                          :damp 1.0))
               (cs (geometric-cs
                    :gamma 1 :smoother smoother :pre-steps 1 :post-steps 0
                    :coarse-grid-iteration
                    (make-instance '<multi-iteration> :base smoother :nr-steps 1)
                    :combination :additive))
               (bpx (make-instance '<cg> :preconditioner cs :restart-cycle 30)))
          (make-instance '<linear-solver>
                         :iteration bpx
                         :success-if `(or (zerop :defnorm) (and (> :step 2) (> :step-reduction 1.0) (< :defnorm 1.0e-8)))
                         :failure-if `(and (> :step 200) (> :step-reduction 1.0) (> :defnorm 1.0e-8))))
	:plot-mesh nil
	:observe
	(append *stationary-fe-strategy-observe*
                (list fl.strategy::*mentries-observe*)
		(list
		 (list (format nil "~19@A~19@A~19@A" "A^00_00" "A^00_11" "A^10_10") "~57A"
		       #'(lambda (blackboard)
			   (let ((tensor (fl.application::effective-tensor blackboard)))
                             (setq saved-effective-tensor tensor)
                             (format nil "~19,10,2E~19,10,2E~19,10,2E"
                                     (and tensor (mref (mref tensor 0 0) 0 0))
                                     (and tensor (> dim 1) (mref (mref tensor 0 0) 1 1))
                                     (and tensor (> dim 1) (mref (mref tensor 1 0) 1 0)))))))))
		   ))))

(defvar *heisig-neuss-2017-demo*
  (make-demo
   :name "Heisig-Neuss-2017"
   :short "Calculations from the paper by M. Heisig and N. Neuss 2017"
   :long
   (format
    nil
    "These are the benchmark calculations from a paper by
M. Heisig and N. Neuss.  This paper reports work done in a
KONWIHR project and was submitted to the ACM JAE in January
2017.

The benchmark problem considered here is calculating a 3D linear
elasticity problem on a representative cell with a hole.  Finite
elements of order 5 are used.

These calculations need a lot of memory, namely about
8^(nr_of_levels-2) GB.  For the serial and shared-memory
parallel calculations this must be provided at startup,
e.g. using the --dynamic-space-size switch for SBCL/Femlisp.

The currently running Femlisp has ~3,2F GB available, thus the
allowed choice for the number of levels will be restricted to
the range 1-~D."
    (* (fl.port:dynamic-space-size) 1e-9)
    (konwihr-paper-max-levels))))

(adjoin-demo *heisig-neuss-2017-demo* *articles-demo*)

(defun konwihr-demo (&optional parallel-p)
  (let* ((max-levels (konwihr-paper-max-levels))
         (query (format nil "Levels (1-~D): " (or max-levels 4)))
         (levels (user-input query #'parse-integer (_ (<= 1 _ max-levels))))
         (thread-query
           (format nil "Threads~@[ [recommended=~D]~]: "
                   (aand (fl.parallel::get-workers) (length it))))
         (threads (and parallel-p 
                       (user-input thread-query #'parse-integer #'plusp))))
    (fl.parallel::end-kernel)
    (when threads (new-kernel threads))
    (initialize-konwihr-paper-calculation)
    (heisig-neuss-2017-demo
     (elasticity-inlay-cell-problem (n-cell-with-ball-hole 3))
     :order 5 :levels levels :output 2)))

(let ((comment
        (_ (format nil "~
Note that performance measurements may vary depending on your
clock frequency.

According to the performance of your computer, you should expect
a time for serial execution of 2^(3*levels) * ~3,1F seconds."
                   (/ 5e3 (konwihr-speed))))))
  (let ((initialization-demo
          (make-demo
           :name "Initialize"
           :short "Initialization of fe/interpolation data"
           :long "
This performs a (mostly silent) FE/multigrid calculation on a
mesh consisting of one cube which is once refined.  This is
sufficient for initializing all data for the performance
measurements of the KONWIHR paper.

If this setup demo has not been performed, the calculation of
this data will be performed before the first test, in order not
to disturb the succeeding tests."
           :execute (_ (initialize-konwihr-paper-calculation))))
        (serial-demo
          (make-demo
           :name "Serial"
           :short "Performance of the serial code"
           :long comment
           :execute (_ (konwihr-demo nil))
           :test-input "1~%"))
        (parallel-demo
          (make-demo
           :name "Shared-memory-parallel"
           :short "Benefits of shared memory parallelization"
           :long comment
           :execute (_ (konwihr-demo t))
           :test-input "1~%2~%")))
    (adjoin-demo initialization-demo *heisig-neuss-2017-demo*)
    (adjoin-demo serial-demo *heisig-neuss-2017-demo*)
    (adjoin-demo parallel-demo *heisig-neuss-2017-demo*)
    ))

