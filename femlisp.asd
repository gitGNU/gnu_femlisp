;;; -*- Lisp -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; femlisp.asd - ASDF definition of the Femlisp system
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

(in-package :cl-user)

(defpackage #:femlisp.system
  (:use #:cl #:asdf))

(asdf::defsystem
    "femlisp"
    :pathname (translate-logical-pathname "femlisp:src;")
    :depends-on ()
    :components
    (
     ;;
     ;; Basic functionality
     ;;
     (:module
      "basic"
      :depends-on ()
      :components
      ((:file "debug")
       (:file "tests")
       (:file "patches")
       (:file "port")
       (:file "multi-processing")
       (:file "macros-defp")
       (:file "macros" :depends-on ("macros-defp"))
       (:file "utilities-defp" :depends-on ("macros-defp"))
       (:file "utilities" :depends-on ("utilities-defp" "macros" "tests"))
       (:file "general" :depends-on ("utilities"))
       (:file "mflop" :depends-on ("utilities"))
       (:file "demo" :depends-on ("tests" "mflop" "macros" "utilities"))
       ))
     (:module
      "matlisp"
      :depends-on ("basic")
      :components
      ((:file "matlisp-defp")
       (:file "vector" :depends-on ("matlisp-defp"))
       (:file "matrix" :depends-on ("vector"))
       (:file "number-blas" :depends-on ("matrix"))
       (:file "standard-matrix" :depends-on ("matrix"))
       (:file "standard-matrix-blas" :depends-on ("standard-matrix"))
       (:file "standard-matrix-lr" :depends-on ("standard-matrix-blas"))
       (:file "array-blas" :depends-on ("matrix"))
       (:file "compat" :depends-on ("matrix"))
       ))
     (:module
      "algebra"
      :depends-on ("basic" "matlisp")
      :components
      ((:file "algebra-defp")
       (:file "tensor" :depends-on ("algebra-defp"))
       (:file "sparse-tensor" :depends-on ("tensor"))
       (:file "crs" :depends-on ("algebra-defp"))
       (:file "sparse" :depends-on ("crs"))
       (:file "sparselu" :depends-on ("sparse"))))
     (:module
      "function"
      :depends-on ("basic" "matlisp" "algebra")
      :components
      ((:file "function-defp")
       (:file "function" :depends-on ("function-defp"))
       (:file "polynom" :depends-on ("function"))
       (:file "spline" :depends-on ("polynom"))))
     ;;
     ;; Mesh
     ;;
     (:module
      "mesh"
      :depends-on ("basic" "matlisp" "algebra" "function")
      :components
      ((:file "mesh-defp")
       (:file "cell" :depends-on ("mesh-defp"))
       (:file "skeleton" :depends-on ("cell"))
       (:file "identify" :depends-on ("skeleton"))
       (:file "refine" :depends-on ("skeleton"))
       (:file "vertex" :depends-on ("refine"))
       (:file "simplex" :depends-on ("vertex"))
       (:file "tensorial" :depends-on ("simplex"))
       (:file "skeleton-build" :depends-on ("tensorial"))
       (:file "domain" :depends-on ("skeleton-build"))
       (:file "mesh" :depends-on ("domain" "identify"))
       (:file "meshgen" :depends-on ("mesh"))
       (:file "extend" :depends-on ("meshgen"))
       ))
     ;;
     ;; Problems
     ;;
     (:module
      "problem"
      :depends-on ("mesh")
      :components
      ((:file "problem-defp")
       (:file "problem" :depends-on ("problem-defp"))
       (:file "pde-problem" :depends-on ("problem"))
       (:file "cdr" :depends-on ("pde-problem"))
       (:file "elasticity" :depends-on ("pde-problem"))
       (:file "navier-stokes" :depends-on ("pde-problem"))
       (:file "cdrsys" :depends-on ("cdr"))
       (:file "time" :depends-on ("pde-problem" "cdr"))
       ))
     ;;
     ;; Discretization
     ;;
     (:module
      "discretization"
      :depends-on ("algebra" "function" "mesh" "problem")
      :components
      ((:file "discretization-defp")
       (:file "discretization" :depends-on ("discretization-defp"))
       (:file "quadrature" :depends-on ("discretization-defp"))
       (:file "fe" :depends-on ("discretization" "quadrature"))
       (:file "lagrange" :depends-on ("fe"))
       (:file "fetransfer" :depends-on ("fe" "lagrange"))
       (:file "sparseas" :depends-on ("fe" "fetransfer"))
       (:file "sparseif" :depends-on ("fe" "sparseas"))
       (:file "feeval" :depends-on ("fe" "sparseif"))
       (:file "constraints" :depends-on ("fe" "sparseif"))
       (:file "fedisc" :depends-on ("constraints"))
       ;;
       (:file "cdr-fe" :depends-on ("fe"))
       (:file "system-fe" :depends-on ("fe" "sparseas"))
       (:file "elasticity-fe" :depends-on ("system-fe"))
       (:file "navier-stokes-fe" :depends-on ("system-fe"))
       ))
     ;;
     ;; Iteration
     ;;
     (:module
      "iteration"
      :depends-on ("basic" "algebra" "function" "problem")
      :components
      ((:file "iteration-defp")
       (:file "iterate" :depends-on ("iteration-defp"))
       (:file "linit" :depends-on ("iteration-defp"))
       (:file "blockit" :depends-on ("linit"))
       (:file "krylow" :depends-on ("linit"))
       (:file "solve" :depends-on ("iteration-defp" "iterate"))
       (:file "linsolve" :depends-on ("linit" "solve"))
       (:file "nlsolve" :depends-on ("solve"))
       (:file "multigrid-defp" :depends-on ("iteration-defp"))
       (:file "multigrid" :depends-on ("linit" "multigrid-defp" "linit"))
       (:file "amg" :depends-on ("multigrid"))
       (:file "selection-amg" :depends-on ("amg"))
       (:file "aggregation-amg" :depends-on ("amg"))
       (:file "stueben" :depends-on ("selection-amg"))))
     (:module
      "special-iteration"
      :pathname #.(translate-logical-pathname #p"femlisp:src;iteration;")
      :depends-on ("mesh" "problem" "discretization" "iteration")
      :components
      ((:file "geomg-defp")
       (:file "geoblock" :depends-on ("geomg-defp"))
       (:file "vanka" :depends-on ("geoblock"))
       (:file "geomg" :depends-on ("geomg-defp"))
       ))
     ;;
     ;; Post-processing
     ;;
     (:module
      "graphic"
      :depends-on ("basic")
      :components
      ((:file "graphics-defp" :depends-on ())
       (:file "graphics" :depends-on ("graphics-defp"))
       (:file "dx" :depends-on ("graphics"))
       (:file "gnuplot" :depends-on ("graphics"))))
     (:module
      "plot"
      :depends-on ("graphic" "mesh" "problem" "discretization")
      :components
      ((:file "plot-defp")
       (:file "plot" :depends-on ("plot-defp"))
       (:file "plot-dx" :depends-on ("plot"))
       (:file "plot-gnuplot" :depends-on ("plot"))
       (:file "meshplot" :depends-on ("plot"))
       (:file "feplot" :depends-on ("plot" "meshplot"))
       (:file "coeffplot" :depends-on ("plot"))
       (:file "asaplot" :depends-on ("plot"))))
     ;;
     ;; Strategy
     ;;
     (:module
      "strategy"
      :depends-on ("mesh" "problem" "discretization" "iteration"
			  "special-iteration" "plot")
      :components
      ((:file "strategy-defp")
       (:file "strategy" :depends-on ("strategy-defp"))
       (:file "strategy-utilities" :depends-on ("strategy"))
       (:file "error-estimator" :depends-on ("strategy-utilities"))
       (:file "error-indicator" :depends-on ("strategy-utilities"))
       (:file "fe-approximation" :depends-on
	      ("strategy" "strategy-utilities" "error-estimator" "error-indicator"))
       (:file "fe-interpolation" :depends-on ("fe-approximation"))
       (:file "fe-stationary" :depends-on ("fe-approximation"))
       (:file "gps" :depends-on ("fe-stationary"))
       ))
     ;;
     ;; Applications
     ;;
     (:module
      "applications"
      :depends-on ("basic" "algebra" "function" "iteration" "mesh"
			   "problem" "discretization" "strategy" "plot")
      :components
      ((:file "application-defp")
       (:file "app-utils" :depends-on ("application-defp"))
       (:module
	"domains"
	:depends-on ("app-utils")
	:components
	((:file "hole-domain")
	 (:file "inlay-domain")
	 (:file "bl-cell")))
       (:module
	"demos"
	:depends-on ("app-utils")
	:components
	((:file "application-demos")
	 (:file "discretization-demos" :depends-on ("application-demos"))
	 (:file "multigrid-demos" :depends-on ("application-demos"))
	 (:file "refinement-demos" :depends-on ("application-demos"))))
       (:module
	"cdr"
	:depends-on ("demos" "domains")
	:components
	((:file "demo-cdr")
	 (:file "tools" :depends-on ("demo-cdr"))
	 (:file "model-problem" :depends-on ("tools"))
	 (:file "unstructured" :depends-on ("tools"))
	 (:file "locref" :depends-on ("tools"))
	 (:file "hom-cdr" :depends-on ("tools"))
	 (:file "bl-cdr" :depends-on ("demo-cdr"))
	 (:file "amg-cdr" :depends-on ("demo-cdr"))
	 (:file "bratu" :depends-on ("demo-cdr"))
	 ))
       (:module
	"elasticity"
	:depends-on ("demos" "domains")
	:components
	((:file "demo-elasticity")
	 (:file "model-problem")
	 (:file "elahom" :depends-on ("demo-elasticity"))
	 ))
       (:module
	"navier-stokes"
	:depends-on ("demos" "domains")
	:components
	((:file "demo-ns")
	 (:file "driven-cavity")
	 (:file "hom-ns")
	 ))
       (:module
	"articles"
	:depends-on ("demos" "domains" "cdr" "elasticity")
	:components
	((:file "heuveline-rannacher-2003")))
       (:module
	"books"
	:depends-on ("demos" "domains" "cdr" "elasticity")
	:components
	())
       (:module
	"courses"
	:depends-on ("demos" "domains" "cdr" "elasticity")
	:components
	())
       (:module
	"talks"
	:depends-on ("demos" "domains" "cdr" "elasticity")
	:components
	((:file "effcoeff")))
       )				; application components
      )					; applications module
     )					; femlisp modules
    )

;;;(asdf:operate 'asdf::load-op 'femlisp)

