;;; -*- mode: lisp; fill-column: 64; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; bl-cdr.lisp
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

;;; In this file we compute effective boundary-layer constants
;;; for oscillating boundaries.  This was treated theoretically
;;; by [Jaeger&Mikelic_1995] and by [Allaire&Amar_1999].  A
;;; further presentation was given in [Neuss_2002].

(defun boundary-layer-cell-problem (bl-cell-domain &key &allow-other-keys)
  "Returns a boundary layer problem for the given boundary layer
cell domain.  We assume that the oscillating boundary on which
Dirichlet nodes are prescribed is identified by $x_n<1$ for its
coordinates, and the patch for which the source is prescribed is
$x_n=0$."
  (let ((dim (dimension bl-cell-domain)))
    (make-instance
     '<cdr-problem> :domain bl-cell-domain
     :patch->coefficients
     #'(lambda (patch)
	 (cond ((and (= (dimension patch) (1- dim))
		     (bl-patch-on-artificial-boundary bl-cell-domain patch))
		(list (scalar-source 1.0)))
	       ((bl-patch-on-lower-boundary bl-cell-domain patch)
		(list (constraint-coefficient 1 1)))
	       ((= (dimension patch) dim)
		(list (scalar-diffusion dim 1.0)))
	       (t nil))))))

(defun sinusoidal-boundary-layer-cell-problem (dim &rest key-args)
  "Returns a boundary layer problem for a sinusoidally
oscillating domain."
  (boundary-layer-cell-problem
   (apply #'sinusoidal-bl-cell dim key-args)))

(defun compute-cbl (blackboard)
  (with-items (&key solution rhs) blackboard
    (and solution rhs (dot solution rhs))))

(defparameter *cbl-observe*
  (list "                Cbl" "~19,10,2E" #'compute-cbl)
  "Observe list for Cbl.")

(defun cdr-bl-computation (dim/domain order max-levels &rest rest
			   &key output plot (plot-mesh t) solver
			   &allow-other-keys)
  "Performs the bl-diffusion demo."
  (let ((*output-depth* output)
	(domain (if (numberp dim/domain)
		    (apply #'sinusoidal-bl-cell dim/domain rest)
		    dim/domain)))
    (defparameter *result*
      (solve (blackboard
	      :problem (boundary-layer-cell-problem domain)
	      :fe-class (lagrange-fe order)
	      :estimator
	      #+(or)(make-instance '<projection-error-estimator>)
	      #-(or)(make-instance '<duality-error-estimator> :functional :load-functional)
	      :indicator
	      #+(or)(make-instance '<uniform-refinement-indicator>)
	      (make-instance '<largest-eta-indicator> :pivot-factor 0.01
			     :from-level 1 :block-p t)
	      :solver solver
	      :plot-mesh plot-mesh
	      :observe (append *stationary-fe-strategy-observe*
			       (list *cbl-observe* *eta-observe*))
	      :success-if `(= :max-level ,(1- max-levels)))))
    (when plot
      (plot (getbb *result* :solution))))
  *result*)

#+(or) (time (cdr-bl-computation
	2 4 3 :plot t :amplitude 0.15 :shift 1.0 :extensible nil
	:solver (lu-solver) :output :all))


#|
(require :sb-sprof)

(sb-sprof:reset)
(sb-sprof:with-profiling
    (:max-samples 10000 :mode :alloc :report :flat :loop t)
  (cdr-bl-computation
    2 4 5 :plot nil :plot-mesh nil :amplitude 0.15 :shift 1.0 :extensible nil
    :solver (lu-solver) :output :all))

 (prof:with-profiling (:type :time)
	 (cdr-bl-computation
	  2 4 3 :plot t :amplitude 0.15 :shift 1.0 :extensible nil
	  :solver (lu-solver) :output :all))
(prof:show-call-graph)
(profile:unprofile)
(profile:report-time)
(profile:reset-time)
(profile:profile sparse-matrix->ccs)
(profile:profile fl.discretization::increment-global-by-local-vec)
(profile:profile fl.discretization::increment-global-by-local-mat)
(profile:profile fl.discretization::get-local-from-global-vec)
(profile:profile fl.discretization::fe-cell-geometry)
(profile:profile fl.discretization::do-fe-dofs-mblocks)
(profile:profile fl.discretization::do-fe-dofs-vblocks)
(profile:profile fl.discretization::fe-cell-geometry)
(profile:profile fl.discretization::assemble-interior)
(profile:profile :methods 'fl.discretization::discretize-locally)
(profile:profile fl.mesh::local->Dglobal)
(profile:profile fl.mesh::local->global)
(profile:profile fl.mesh::l2g)
(profile:profile fl.mesh::l2Dg)
(profile:profile fl.matlisp::gesv!)
(profile:profile fl.matlisp::gemm!)
(profile:profile fl.mesh::euclidean->barycentric)
(profile:profile fl.mesh::weight-vector-product-cell)
(profile:profile fl.mesh::corners)
(profile:profile fl.mesh::vertices)
(profile:profile fl.application::bottom-mapping)
(profile:profile fl.mesh::weight-lists-grad-product-cell)
|#

;; before change (+ 1 sec on another run)
;;    2        36       784    1.8   9.3593619483d-01   6.6484919470d-04
;;    8       172      4462   10.9   9.4073588298d-01   1.6845432818d-04
;;   20       444     12426   34.6   9.4064954875d-01   2.4130130651d-06

(defun make-cdr-bl-demo (dim order levels)
  (let ((title (format nil "bl-diffusion-~Dd" dim))
	(short "Computes a boundary law coefficient")
	(long "Computes a boundary law coefficient constant
for an oscillating function which has to be 1-periodic and
negative.  This involves solving a cell problem on a
semi-infinite domain with a smooth and exponentially decaying
solution.  The constant we search for is equal to the energy of
the solution.  Our adaptive strategy with high order
approximations and local mesh refinement is well suited for
computing this constant with high accuracy.  Error estimation is
achieved by using a duality error estimator for the load
functional."))
    (let ((demo
	   (make-demo
	    :name (format nil title dim) :short short
	    :long (format nil "~A~%~%Parameters: dim=~D, p=~D, ~D levels~%"
			  long dim order levels)
	    :execute (lambda () (cdr-bl-computation dim order levels :plot t :output 1)))))
      (adjoin-demo demo *laplace-demo*)
      (adjoin-demo demo *adaptivity-demo*)
      (adjoin-demo demo *boundary-coeffs-demo*))))

;;; 2D and 3D demo
(make-cdr-bl-demo 2 4 3)
(make-cdr-bl-demo 3 3 2)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Testing
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun test-bl-cdr ()
  
  (loop
   repeat 2 do
   (cdr-bl-computation
    2 4 3 :plot t :amplitude 0.15 :extensible nil :output :all)
   finally
   (let ((work (* (common-lisp-speed) (getbb *result* :time))))
     ;; ensure that we have no dramatic performance drop.
     ;; @var{work} was about 1000 for Femlisp-0.9.4 on ortler
     (format t "~F ~F  ->  ~F~%"
	     (common-lisp-speed :cache 0.0 :memory 1.0)
	     (common-lisp-speed :cache 1.0 :memory 0.0)
	     work)
     ;;; 0.5/0.5 -> 2500 on ortler
     (when (> work 2500)
       (error "Performance problem: Work=~A, expected 1000+/-10%" work))))
  
  (multiple-value-bind (f Df)
      (cubic-spline #(1.2 1.2 1.2))
    (describe 
     (oscillating-boundary-domain
     2 f :grad-f Df)))
  
  (cdr-bl-computation (spline-interpolated-bl-cell #(1.2 1.2 1.2))
		      2 2 :plot t :output :all)
  
  ;; testing if identification with only one cell width works
  (let* ((problem (sinusoidal-boundary-layer-cell-problem 2))
	 (mesh (make-hierarchical-mesh-from-domain (domain problem))))
    (refine mesh :indicator #'(lambda (cell) (<= (aref (midpoint cell) 1) 0.0)))
    (multiple-value-bind (mat rhs)
	(discretize-globally problem mesh (lagrange-fe 1))
      #+(or)(show mat)
      #-(or)
      (plot (linsolve mat rhs :iteration (make-instance '<lu>))
	    :depth 3)
      ))

  (plot (getbb *result* :mesh))
  (plot (fl.strategy::eta->p2-vec (getbb *result* :eta) (getbb *result* :problem)
			       (getbb *result* :mesh)))
  (display-ht (getbb *result* :eta))
  (apply #'fl.strategy::compute-local-estimate
	 (slot-value (getbb *result* :strategy) 'fl.strategy::estimator) *result*)

  (dohash (key (getbb *result* :eta))
    (let ((mp (midpoint key)))
      (when (< (norm (m- mp #(0.9375 -0.3125))) 1.0e-10)
	(display-ht (matrix-row (getbb *result* :matrix) key)))))
 
  (plot (getbb *result* :solution))

  ;; reiteration
  (solve (s1-reduction-amg-solver 4 :reduction 1.0e-3) *result*)
  (plot (getbb *result* :solution))
  (plot (getbb *result* :mesh))
  (getbb *result* :global-eta)
  (show (getbb *result* :solution))
  (plot (getbb *result* :solution) :depth 2 :plot :file :filename "bl-cell-sol" :format "eps")
  (plot (mesh (getbb *result* :solution)) :plot :file :filename "bl-cell-mesh.ps" :format "eps")
  (compute-cbl *result*)
  (nr-of-levels (mesh (getbb *result* :solution)))
  (nr-of-cells (mesh (getbb *result* :solution)))
  (destructuring-bind (&key matrix solution rhs &allow-other-keys)
      *result*
    (plot (m- rhs (m* matrix solution))))
  
  (describe (sinusoidal-boundary-layer-cell-problem 2))
  )

;;; (test-bl-cdr)
;;; (fl.tests:adjoin-test 'test-bl-cdr)