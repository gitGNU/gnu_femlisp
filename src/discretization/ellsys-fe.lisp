;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; ellsys-fe.lisp - Discretization of general elliptic systems
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;
;;; Copyright (C) 2006 Nicolas Neuss, University of Karlsruhe.
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
;;; NO EVENT SHALL THE AUTHOR, THE UNIVERSITY OF KARLSRUHE OR OTHER
;;; CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
;;; EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
;;; PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
;;; PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
;;; LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
;;; NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
;;; SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(in-package "COMMON-LISP-USER")
(defpackage "FL.ELLSYS-FE"
  (:use "COMMON-LISP" "FL.MACROS" "FL.UTILITIES" "FL.DEBUG"
	"FL.MATLISP" "FL.FUNCTION"
	"FL.MESH" "FL.PROBLEM" "FL.ELLSYS"
	"FL.DISCRETIZATION")
  (:export)
  (:documentation "Finite element discretization of a general second order
elliptic system, see the description in the ELLSYS package.  The result is
a local matrix A and local rhs b.  They will usually depend on u_old which
is stored in the solution vector."))

(in-package :fl.ellsys-fe)

(defun ensure-m-matrix-property (mat)
  "Problematic, because, when applied to the convection contribution,
it reduces the discretization to first order accuracy."
  (for-each-row-key
   (lambda (i)
     (for-each-key-in-row
      (lambda (j)
        (unless (= j i)
          (let ((entry (mref mat i j)))
            (when (plusp entry)
              (incf (mref mat i i) entry)
              (setf (mref mat i j) 0.0)))))
      mat i))
   mat)
  mat)

(defparameter *upwinding* nil
  "Discretizes convection with upwinding when set to T.")

(defmethod discretize-locally
    ((problem <ellsys-problem>) coeffs vecfe qrule fe-geometry
     &key matrix rhs mass-matrix local-u local-v residual-p fe-parameters
     &allow-other-keys)
  "Local discretization for a pde system of the form described in the
documentation of the package @package{ELLSYS}."
  (declare (optimize debug))
  (assert (and (null local-u) (null local-v)) () "NYI")

  (let ((cell (getf fe-geometry :cell))
        (ip-values (ip-values vecfe qrule))
        (ip-gradients
          (and (intersection
                (mapcar #'coefficient-name coeffs)
                '(FL.ELLSYS::A FL.ELLSYS::AI FL.ELLSYS::B FL.ELLSYS::C  FL.ELLSYS::H))
               (ip-gradients vecfe qrule)))
        (gradients (getf fe-geometry :gradients))
        (gradient-inverses (getf fe-geometry :gradient-inverses)))
    (dbg :disc "Number of quadrature points = ~A" (length ip-values))

    ;; loop over quadrature points
    (loop
       for i from 0
       and global across (getf fe-geometry :global-coords)
       and weight across (getf fe-geometry :weights)
       do
          (let* (
                 (Dphi (and gradients (aref gradients i)))
                 (Dphi^-1 (and gradients (aref gradient-inverses i)))
                 ;; nr-comps x (n-basis x 1)
                 (shape-vals (aref ip-values i))
                 ;; nr-comps x (n-basis x dim)
                 (shape-grads (and ip-gradients (aref ip-gradients i)))
                 (gradients (and Dphi^-1 shape-grads (map 'vector (rcurry #'m* Dphi^-1) shape-grads)))
                 (coeff-input (construct-coeff-input
                               cell global Dphi shape-vals gradients fe-parameters))
                 (left-vals shape-vals)
                 (right-vals (if residual-p
                                 (map 'vector #'m* shape-vals local-u)
                                 shape-vals))
                 (right-gradients gradients)
                 (left-gradients (if residual-p
                                     (map 'vector #'m* gradients local-u)
                                     gradients)))
            (assert (vectorp left-vals))
            (loop
              for coeff in coeffs
              for name = (coefficient-name coeff) do
                ;; matrix part
                (when matrix
                  (case name
                    ;; diffusion
                    (FL.ELLSYS::A
                     (let ((diff-tensor (evaluate coeff coeff-input)))
                       (dbg :disc "Discretizing diffusion")
                       (for-each-entry-and-key
                        #'(lambda (diffusion i j)
                            (gemm! weight (m* (aref left-gradients j) diffusion) (aref right-gradients i)
                                   1.0 (mref matrix i j) :NT))
                        diff-tensor)))
                    ;; resistance
                    (FL.ELLSYS::AI
                     (let ((resistance-tensor (evaluate coeff coeff-input)))
                       (dbg :disc "Discretizing resistance")
                       (for-each-entry-and-key
                        #'(lambda (resistance i j)
                            (let ((factor (m/ (m*-tn Dphi (m* resistance Dphi)))))
                              (gemm! weight (m* (aref shape-grads j) factor) (aref shape-grads i)
                                     1.0 (mref matrix i j) :NT)))
                        resistance-tensor)))
                    (FL.ELLSYS::B
                     (dbg :disc "Discretizing convection")
                     (let ((velocity-tensor (evaluate coeff coeff-input)))
                       (for-each-entry-and-key
                        #'(lambda (velocity i j)
                            (let ((convection-contribution
                                    (m*-nt (m* (aref left-gradients i) (scal -1.0 velocity))
                                           (aref right-vals j))))
                              (when *upwinding*
                                (scal! 2.0 (ensure-m-matrix-property convection-contribution)))
                              (axpy! weight convection-contribution
                                     (mref matrix i j))))
                        velocity-tensor)))
                    (FL.ELLSYS::C
                     (dbg :disc "Discretizing non-conservative convection")
                     (let ((velocity-tensor (evaluate coeff coeff-input)))
                       (for-each-entry-and-key
                        #'(lambda (velocity i j)
                            (if *upwinding*
                                (let ((contribution
                                        (m* (aref left-vals i)
                                            (transpose (m* (aref right-gradients j) velocity)))))
                                  ;; enhancement factor for consistency
                                  (axpy! (* 1.0
                                            weight)
                                         (ensure-m-matrix-property contribution)
                                         (mref matrix i j)))
                                (gemm! weight (aref left-vals i)
                                       (m* (aref right-gradients j) velocity)
                                       1.0 (mref matrix i j) :NT)))
                        velocity-tensor)))
                    (FL.ELLSYS::R
                     (dbg :disc "Discretizing reaction")
                     (let ((reaction-tensor (evaluate coeff coeff-input)))
                       (for-each-entry-and-key
                        #'(lambda (reaction i j)
                            (unless (mzerop reaction)
                              (when (standard-matrix-p reaction)
                                (assert (= 1 (nrows reaction) (ncols reaction)))
                                (setq reaction (vref reaction 0)))
                              (gemm! (* weight reaction) (aref left-vals i) (aref right-vals j)
                                     1.0 (mref matrix i j) :NT)))
                        reaction-tensor)))))
                (when rhs
                  (case name
                    (FL.ELLSYS::F
                     (dbg :disc "Discretizing standard source")
                     (let ((source-vector (evaluate coeff coeff-input)))
                       (for-each-entry-and-key
                        #'(lambda (source k)
                            (gemm! weight (aref left-vals k) (ensure-matlisp source)
                                   1.0 (aref rhs k)))
                        source-vector)))
                    (FL.ELLSYS::G
                     (dbg :disc "Discretizing g-source")
                     (let ((g-source-vector (evaluate coeff coeff-input)))
                       (for-each-entry-and-key
                        #'(lambda (g-source k)
                            (gemm! weight (aref left-gradients k) g-source
                                   1.0 (aref rhs k)))
                        g-source-vector)))
                    (FL.ELLSYS::H
                     (dbg :disc "Discretizing h-source")
                     (let ((diff-tensor (evaluate (get-coefficient coeffs 'FL.ELLSYS::A) coeff-input))
                           (gamma (evaluate coeff coeff-input)))
                       (for-each-entry-and-key
                        #'(lambda (diffusion i j)
                            (gemm! weight (m* (aref left-gradients j) diffusion) (aref gamma j)
                                   1.0 (aref rhs i)))
                        diff-tensor)))))
              )
            (when mass-matrix
              (dbg :disc "Discretizing mass matrix")
              (loop for i from 0 and vals across shape-vals do
                (gemm! weight vals vals 1.0 (mref mass-matrix i i) :NT)))
            ))
    ;; nonstandard rhs (outside of ip loop)
    (whereas ((fe-rhs-function (and rhs (get-coefficient coeffs 'FL.ELLSYS::FE-RHS))))
             (m+! (evaluate fe-rhs-function (list :cell cell :fe vecfe :geometry fe-geometry))
                  rhs))
    ))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; GPS interface
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmethod select-discretization ((problem <ellsys-problem>) blackboard)
  (declare (ignore blackboard))
  (let ((dim (dimension (domain problem))))
    (lagrange-fe (or *suggested-discretization-order*
		     (cond ((<= dim 2) 4)
                           ((<= dim 3) 3)
                           (t 1)))
		 :nr-comps (nr-of-components problem))))

;;; Testing

(defun discretize-ellsys (dim nr-comps order level)
  "Discretizes an elliptic system with @arg{nr-comps} components on the
@arg{dim}-dimensional unit cube with Lagrange fe of @arg{order}."
  (let* ((problem (ellsys-model-problem
		   dim `((u ,nr-comps))
		   :a (isotropic-diffusion dim #(1.0))
		   :f (coerce (loop repeat nr-comps collect #m(1.0)) 'vector)
		   :dirichlet (constraint-coefficient nr-comps 1)))
	 (h-mesh (uniformly-refined-hierarchical-mesh (domain problem) level))
	 (fedisc (lagrange-fe order :nr-comps nr-comps)))
    (time (discretize-globally problem h-mesh fedisc))))

;;;; Testing

(defun ellsys-fe-tests ()
  (dbg-on :disc)
  (multiple-value-bind (matrix rhs)
      (discretize-ellsys 1 1 1 2)
    (fl.matlisp:show matrix)
    (fl.matlisp:show rhs)
    (assert (zerop (norm (m- (fe-value (gesv matrix rhs) #(0.5))
			     (vector #m(0.125)))))))
  (dbg-off)

  (let* ((dim 1) (order 1) (level 2)
	 (problem (fl.cdr::cdr-model-problem dim :diffusion nil))
	 (h-mesh (uniformly-refined-hierarchical-mesh (domain problem) level))
	 (fedisc (lagrange-fe order)))
    (discretize-globally problem h-mesh fedisc))
  
  (time (multiple-value-bind (matrix rhs)
	    (discretize-ellsys 2 1 1 5)
	  (fe-value (gesv matrix rhs) #d(0.5 0.5))))

  #+(or)  ; plot is not yet available when loading
  (multiple-value-bind (matrix rhs)
      (discretize-ellsys 2 1 1 5)
    (fl.plot:plot (gesv matrix rhs)))
  
  #+(or)
  (multiple-value-bind (matrix rhs)
      (time (discretize-elasticity 2 1 5))
    (fl.plot:plot (gesv matrix rhs) :component 0))

  (let* ((dim 1) (levels 1) (order 1)
         (problem (ellsys-model-problem
                   dim '((u 1))
                   :r (diagonal-sparse-tensor (vector #m(1.0)))
                   :f (lambda (x)
                        (vector (ensure-matlisp (float #I"sin(2*pi*x[0])" 1.0))))))
         (h-mesh (uniformly-refined-hierarchical-mesh (domain problem) levels))
         (fe-class (lagrange-fe order :nr-comps 1))
         (ansatz-space (make-fe-ansatz-space fe-class problem h-mesh)))
    (loop for i below 1000 do
         (let ((rhs (make-ansatz-space-vector ansatz-space)))
           (assemble-interior ansatz-space :surface :level 1 :rhs rhs)
           (let ((val (fe-value rhs #d(0.5))))
             (unless (mzerop val 1.0e-15)
               (format t "~A ~A~%" (mzerop val 1.0e-15) val)
               (error "Should not happen"))))))
  )

;;; (ellsys-fe-tests)
(fl.tests:adjoin-test 'ellsys-fe-tests)
