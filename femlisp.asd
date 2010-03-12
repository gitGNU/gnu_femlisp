;;; -*- Lisp -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; femlisp.asd - Femlisp system definition (source files)
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

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Setup the logical host "FEMLISP"
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(in-package :cl-user)

(eval-when (:load-toplevel :execute)
  (let ((directory (pathname-directory *load-truename*)))
    (setf (logical-pathname-translations "FEMLISP")
	  `((#-gcl "**;*.*.*" #+gcl "**;*.*" 
	     ,(make-pathname :directory `(,@directory :wild-inferiors)
			     :name :wild :type :wild :version :wild))))))

(in-package :asdf)

(defsystem :femlisp-basic
  :depends-on (#+(or clisp ccl) :cffi #+(or clisp ccl) :closer-mop)
  :components
  ((:file "setup")
   (:module
    "basic"
    :depends-on ("setup")
    :pathname #.(translate-logical-pathname #p"femlisp:src;basic;")
    :components
    ((:file "debug" :depends-on ())
     (:file "tests" :depends-on ())
     (:file "patches" :depends-on ())
     (:file "macros" :depends-on ())
     (:file "port" :depends-on ("debug"))
     (:file "amop" :depends-on ("debug" "port"))
     (:file "utilities-defp" :depends-on ("macros" "debug"))
     (:file "utilities" :depends-on ("utilities-defp" "macros" "tests"))
     (:file "mflop" :depends-on ("utilities"))
     (:file "multiprocessing" :depends-on ("utilities" "mflop" "port" "debug"))
     (:file "parcells" :depends-on ("multiprocessing"))
     (:file "general" :depends-on ("amop" "utilities-defp" "multiprocessing"))
     (:file "demo" :depends-on ("tests" "mflop" "macros" "utilities"))
     ))
   ))

(defsystem :femlisp-matlisp
  :depends-on (:femlisp-basic)
  :components
  ((:module
    "alien"
    :depends-on ()
    :pathname #.(translate-logical-pathname "femlisp:src;alien;")
    :components
    (;; alien should be recompiled if superlu.so or umfpack.so has changed
     (:file "alien" :depends-on ())
     (:file "alienc" :depends-on ("alien"))
     (:file "lapack" :depends-on ("alien"))
     (:file "superlu" :depends-on ("alien"))
     (:file "umfpack" :depends-on ("alien"))
     ))
   (:module
    "matlisp"
    :depends-on ("alien")
    :pathname #.(translate-logical-pathname "femlisp:src;matlisp;")
    :components
    ((:file "matlisp-defp")
     (:file "ctypes" :depends-on ("matlisp-defp"))
     (:file "vector" :depends-on ("matlisp-defp"))
     (:file "blas-basic" :depends-on ("vector"))
     (:file "matrix" :depends-on ("vector"))
     (:file "number-blas" :depends-on ("matrix" "blas-basic"))
     (:file "array-blas" :depends-on ("matrix" "ctypes" "number-blas"))
     (:file "store-vector" :depends-on ("matrix" "ctypes" "blas-basic"))
     (:file "standard-matrix" :depends-on ("matrix" "store-vector"))
     (:file "standard-matrix-blas" :depends-on ("standard-matrix"))
     (:file "standard-matrix-lr" :depends-on ("standard-matrix-blas"))
     (:file "compat" :depends-on ("standard-matrix"))
     (:file "call-matlisp" :depends-on ("standard-matrix"))
     (:file "tensor" :depends-on ("store-vector" "standard-matrix"))
     (:file "sparse-tensor" :depends-on ("tensor"))
     (:file "compressed" :depends-on ("store-vector" "standard-matrix"))
     (:file "ggev" :depends-on ("standard-matrix"))
     (:file "hegv" :depends-on ("standard-matrix"))
     (:file "sparse-vector" :depends-on ("compat" "sparse-tensor"))
     (:file "sparse-matrix" :depends-on ("sparse-vector" "standard-matrix-blas" "compressed"))
     (:file "sparselu" :depends-on ("sparse-matrix"))
     ))
   ))

(defsystem :femlisp
  :depends-on (:femlisp-basic :femlisp-matlisp)
  :components
  ((:module
    "function"
    :depends-on ()
    :pathname #.(translate-logical-pathname "femlisp:src;function;")
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
    :pathname #.(translate-logical-pathname "femlisp:src;mesh;")
    :depends-on ( "function")
    :components
    ((:file "mesh-defp")
     (:file "cell" :depends-on ("mesh-defp"))
     (:file "skeleton" :depends-on ("cell"))
     (:file "identify" :depends-on ("skeleton"))
     (:file "refine" :depends-on ("identify"))
     (:file "vertex" :depends-on ("refine"))
     (:file "simplex" :depends-on ("vertex"))
     (:file "product-cell" :depends-on ("simplex"))
     (:file "skeleton-build" :depends-on ("product-cell" "identify"))
     (:file "domain" :depends-on ("skeleton-build"))
     (:file "mesh" :depends-on ("domain"))
     (:file "meshgen" :depends-on ("mesh"))
     (:file "triangulate" :depends-on ("meshgen"))
     (:file "triangle" :depends-on ("triangulate"))
     (:file "extend" :depends-on ("meshgen"))
     ))
   ;;
   ;; Problems
   ;;
   (:module
    "problem"
    :depends-on ("mesh")
    :pathname #.(translate-logical-pathname "femlisp:src;problem;")
    :components
    ((:file "problem-defp")
     (:file "problem" :depends-on ("problem-defp"))
     (:file "pde-problem" :depends-on ("problem"))
     (:file "pdef" :depends-on ("pde-problem"))
     (:file "evp" :depends-on ("pdef"))
     (:file "time" :depends-on ("pdef"))
     (:file "ellsys" :depends-on ("evp" "time" "pdef"))
     (:file "cdr" :depends-on ("ellsys"))
     ;;(:file "cdrsys" :depends-on ("cdr"))
     (:file "elasticity" :depends-on ("ellsys"))
     (:file "navier-stokes" :depends-on ("ellsys"))
     ;;(:file "navier-stokes-ellsys" :depends-on ("ellsys"))
     ))
   ;;
   ;; Iteration
   ;;
   (:module
    "iteration"
    :depends-on ("function" "problem")
    :pathname #.(translate-logical-pathname "femlisp:src;iteration;")
    :components
    ((:file "iteration-defp")
     (:file "iterate" :depends-on ("iteration-defp"))
     (:file "linit" :depends-on ("iteration-defp"))
     (:file "blockit" :depends-on ("linit"))
     (:file "krylow" :depends-on ("linit"))
     (:file "solve" :depends-on ("iteration-defp" "iterate"))
     (:file "linsolve" :depends-on ("linit" "solve"))
     (:file "nlsolve" :depends-on ("linsolve"))
     (:file "evpsolve" :depends-on ("nlsolve"))
     (:file "multigrid-defp" :depends-on ("iteration-defp"))
     (:file "multigrid" :depends-on ("linit" "linsolve" "multigrid-defp"))
     (:file "amg" :depends-on ("multigrid"))
     (:file "selection-amg" :depends-on ("amg"))
     (:file "aggregation-amg" :depends-on ("amg"))
     (:file "stueben" :depends-on ("selection-amg"))))
   ;;
   ;; Discretization
   ;;
   (:module
    "discretization"
    :depends-on ("function" "mesh" "problem" "iteration")
    :pathname #.(translate-logical-pathname "femlisp:src;discretization;")
    :components
    ((:file "discretization-defp")
     (:file "discretization" :depends-on ("discretization-defp"))
     (:file "quadrature" :depends-on ("discretization-defp"))
     (:file "fe" :depends-on ("discretization" "quadrature"))
     (:file "ansatz-space" :depends-on ("fe"))
     (:file "lagrange" :depends-on ("ansatz-space"))
     (:file "fetransfer" :depends-on ("ansatz-space" "lagrange"))
     (:file "sparseas" :depends-on ("ansatz-space" "fetransfer"))
     (:file "sparseif" :depends-on ("ansatz-space" "sparseas"))
     (:file "feeval" :depends-on ("ansatz-space" "sparseif"))
     (:file "constraints" :depends-on ("ansatz-space" "sparseif"))
     (:file "fedisc" :depends-on ("constraints"))
     ;;
     ;;(:file "cdr-fe" :depends-on ("fedisc"))
     (:file "system-fe" :depends-on ("fedisc"))
     (:file "ellsys-fe" :depends-on ("fedisc"))
     (:file "elasticity-fe" :depends-on ("system-fe"))
     (:file "navier-stokes-fe" :depends-on ("system-fe"))
     ;;(:file "cdrsys-fe" :depends-on ("system-fe"))
     ))
   (:module
    "special-iteration"
    :depends-on ("mesh" "problem" "discretization" "iteration")
    :pathname #.(translate-logical-pathname #p"femlisp:src;iteration;")
    :components
    ((:file "geomg-defp")
     (:file "geoblock" :depends-on ("geomg-defp"))
     (:file "vanka" :depends-on ("geoblock"))
     (:file "geomg" :depends-on ("geoblock"))
     ))
   ;;
   ;; Post-processing
   ;;
   (:module
    "graphic"
    :depends-on ()
    :pathname #.(translate-logical-pathname #p"femlisp:src;graphic;")
    :components
    ((:file "graphics-defp" :depends-on ())
     (:file "graphics" :depends-on ("graphics-defp"))
     (:file "dx" :depends-on ("graphics"))
     (:file "gnuplot" :depends-on ("graphics"))))
   (:module
    "plot"
    :depends-on ("graphic" "mesh" "problem" "discretization")
    :pathname #.(translate-logical-pathname #p"femlisp:src;plot;")
    :components
    ((:file "plot-defp")
     (:file "plot" :depends-on ("plot-defp"))
     (:file "plot-dx" :depends-on ("plot"))
     (:file "plot-gnuplot" :depends-on ("plot"))
     (:file "meshplot" :depends-on ("plot-dx"))
     (:file "feplot" :depends-on ("meshplot"))
     (:file "coeffplot" :depends-on ("feplot"))
     (:file "function-plot" :depends-on ("plot-dx"))
     (:file "asaplot" :depends-on ("plot-dx"))))
   ;;
   ;; Strategy
   ;;
   (:module
    "strategy"
    :depends-on ("mesh" "problem" "discretization" "iteration"
			"special-iteration" "plot")
    :pathname #.(translate-logical-pathname #p"femlisp:src;strategy;")
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
     (:file "fe-evp" :depends-on ("fe-stationary"))
     (:file "rothe" :depends-on ("fe-stationary"))
     (:file "rothe-ellsys" :depends-on ("rothe"))
     (:file "gps" :depends-on ("fe-stationary"))
     ))
   ;;
   ;; Applications
   ;;
   (:module
    "applications"
    :depends-on ("function" "iteration" "mesh"
		 "problem" "discretization" "strategy" "plot")
    :pathname #.(translate-logical-pathname #p"femlisp:src;applications;")
    :components
    ((:module
      "domains"
      :depends-on ()
      :components
      ((:file "domains-defp")
       (:file "hole-domain" :depends-on ("domains-defp"))
       (:file "inlay-domain" :depends-on ("hole-domain"))
       (:file "bl-cell" :depends-on ("domains-defp"))))
     (:file "application-defp" :depends-on ("domains"))
     (:file "app-utils" :depends-on ("application-defp"))
     (:module
      "demos"
      :depends-on ("app-utils" "domains")
      :components
      ((:file "application-demos")
       (:file "discretization-demos" :depends-on ("application-demos"))
       (:file "multigrid-demos" :depends-on ("application-demos"))
       (:file "refinement-demos" :depends-on ("application-demos"))
       (:file "problem-demos" :depends-on ("application-demos"))))
     (:module
      "cdr"
      :depends-on ("demos" "domains")
      :components
      ((:file "tools")
       (:file "model-problem" :depends-on ("tools"))
       (:file "unstructured" :depends-on ("tools"))
       (:file "locref" :depends-on ("tools"))
       (:file "hom-cdr" :depends-on ("tools"))
       (:file "bl-cdr")
       (:file "mg-cdr")
       (:file "amg-cdr" :depends-on ("tools"))
       (:file "bratu")
       (:file "evp-cdr")
       (:file "heat")
       ))
     (:module
      "elasticity"
      :depends-on ("demos" "domains")
      :components
      ((:file "model-problem")
       (:file "elahom")
       ))
     (:module
      "navier-stokes"
      :depends-on ("demos" "domains")
      :components
      ((:file "driven-cavity")
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
     )					; application components
    )					; applications module
   )					; femlisp modules
  )
