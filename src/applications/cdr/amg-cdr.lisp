;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; amg-cdr.lisp - Solving CDR problems with AMG
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

(in-package :application)

;;; Testing the S1-reduction AMG
(time
 (let* ((dim 1) (level 2) (order 2)
	(problem (laplace-test-problem-on-domain (n-cube-domain dim))))
   (multiple-value-bind (A b)
       (problem-discretization problem :level level :order order)
     (let ((amg (s1-reduction-amg-solver order :reduction 1.0d-10 :output t)))
       (solve amg :matrix A :rhs b :output t)))))

;;; k-fold jump in refinement depth
(time
 (let* ((dim 2) (order 2) (k 2)
	(problem (laplace-test-problem-on-domain (n-cube-domain dim)))
	(fedisc (lagrange-fe order))
	(h-mesh (make-hierarchical-mesh-from-domain (domain problem))))
   #-(or) (loop repeat 1 do (refine h-mesh))
   #-(or) (loop repeat 2 do (refine h-mesh :test (rcurry #'inside-cell? (make-double-vec dim 0.25))))
   #-(or) (plot h-mesh)
   #-(or)
   (multiple-value-bind (mat rhs)
       (discretize-globally problem h-mesh fedisc)
     (let ((amg (s1-reduction-amg-solver order :output t :maxsteps 10 :reduction 1.0e-10)))
       (solve amg :matrix mat :rhs rhs)))))

(let* ((mat multigrid::mat)
       (rks (row-keys mat))
       (cks (col-keys mat))
       (rks (append cks (set-difference rks cks))))
  (display mat :col-order cks :row-order rks)
  rks)
  
(time
 (let* ((dim 2) (order 2)
	(problem (laplace-test-problem-on-domain (n-cube-domain dim)))
	(fedisc (lagrange-fe order))
	(h-mesh (make-hierarchical-mesh-from-domain (domain problem))))
   #-(or) (loop repeat 1 do (refine h-mesh))
   #-(or) (loop repeat 2 do (refine h-mesh :test (rcurry #'inside-cell? #(0.25 0.25))))
   #+(or) (plot h-mesh)
   #-(or)
   (multiple-value-bind (mat rhs)
       (discretize-globally problem h-mesh fedisc)
     (let ((amg (s1-reduction-amg-solver order :output t :reduction 1.0e-5 :maxsteps 5)))
       (solve amg :matrix mat :rhs rhs)))))

    

;;; AMG cycle for higher order discretizations
(let* ((dim 2) (level 2) (order 3)
       ;;(amg (make-instance '<s1-reduction> :max-depth 2))
       (problem (laplace-test-problem-on-domain (n-cube-domain dim))))
  (multiple-value-bind (A b)
      (problem-discretization problem :level level :order order)
    #+(or)
    (let ((mg-data (multilevel-decomposition amg A)))
      (show (aref (get-al mg-data :a-vec) 0)))
    #+(or)
    (let ((mg-data (multilevel-decomposition amg A)))
      (show (aref (get-al mg-data :i-vec) 0)))
    #+(or) (linsolve A b :output t :iteration amg :maxsteps 10)
    #-(or) (solve (amg-solver order) :matrix A :rhs b)
    ))

