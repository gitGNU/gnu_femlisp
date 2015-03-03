;;; -*- mode: lisp; fill-column: 64; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; sturm.lisp - An electromagnetic potential calculation
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;
;;; Copyright (C) 2010 Nicolas Neuss, Karlsruhe Institute of Technology
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
;;; IN NO EVENT SHALL THE AUTHOR, THE KARLSRUHE INSTITUTE OF TECHNOLOGY,
;;; OR OTHER CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
;;; SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
;;; LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
;;; DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
;;; THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
;;; (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
;;; OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(in-package :fl.application)

(defun periodic-polygonal (positions)
  "Constructs a periodic polygon from @arg{positions}.  The result is a
@class{<skeleton>}."
  (let* ((vertices (mapcar #'make-vertex positions))
         (lines (loop for (a b) on vertices collect
                     (make-line a (or b (first vertices))))))
  (make-instance '<skeleton> :dimension 1 :cells lines)))

(defun sturm-domain ()
  (let* ((outer-boundary
          (periodic-polygonal
           (list #d(0.0 0.0) #d(1.0 0.0) #d(1.0 1.0) #d(0.0 1.0))))
         (inner-boundary
          (periodic-polygonal
           (list #d(0.2 0.2) #d(0.8 0.2) #d(0.8 0.8) #d(0.2 0.8))))
         (holes '((:A 0.1 #(0.65 0.65)) (:B 0.1 #(0.35 0.35))))
         (hole-boundaries
           (loop for (name radius center) in holes collecting
                (make-instance
                 '<skeleton> :cells
                 (list (let* ((ball-mapping (circle-function radius center (* 2 pi)))
                              (vtx (make-vertex (evaluate ball-mapping #d(0.0)))))
                         (make-line vtx vtx :mapping ball-mapping))))))
         (outer-region-patch
          (make-instance '<boundary-cell>
                         :dimension 2
                         :boundary (concatenate 'vector
                                                (cells-of-dim outer-boundary 1)
                                                (cells-of-dim inner-boundary 1))
                         :midpoint #d(0.9 0.9)
                         :holes '(#(0.5 0.5))))
         (inner-region-patch
          (make-instance '<boundary-cell>
                         :dimension 2
                         :boundary (coerce (mappend #'cells-of-highest-dim
                                                    (list* inner-boundary hole-boundaries))
                                           'vector)
                         :midpoint #d(0.3 0.7)
                         :holes (mapcar #'third holes)))
         (result (make-instance '<skeleton> :dimension 2)))
    (skel-add! result (mark-skeleton outer-boundary :part :outer-boundary))
    (skel-add! result (mark-skeleton inner-boundary :part :inner-boundary))
    (loop for hole in holes and bdry in hole-boundaries do
         (skel-add! result (mark-skeleton bdry :part (first hole))))
    (setf (skel-ref result outer-region-patch) '(:part :outer))
    (setf (skel-ref result inner-region-patch) '(:part :inner))
    (change-class
     result '<domain>
     :classifiers
     (list (lambda (patch classifiers)
             (awhen (get-cell-property patch result :name)
               (pushnew it classifiers))
             classifiers)))))

(defun sturm-mesh (domain)
  (let ((mesh (triangulate domain :threshold 0.2 :meshsize .1)))
    (check mesh)
    (change-class mesh '<hierarchical-mesh>)))

(defun sturm-problem (domain)
  (create-problem
      '<cdr-problem>
      (:domain domain :components '(u) :multiplicity 1)
    (setup-coefficients (patch)
      (select-on-patch (patch)
        (:A                             ; hole A
         (coeff FL.CDR::SCALAR-CONSTRAINT () 1.0))
        (:B                             ; hole B
         (coeff FL.CDR::SCALAR-CONSTRAINT () -1.0))
        ((and :d-dimensional :outer)
         (coeff FL.CDR::ISOTROPIC-DIFFUSION () 10.0))
        ((and :d-dimensional :inner)
         (coeff FL.CDR::ISOTROPIC-DIFFUSION () 1.0))
        ))))

(defun solve-sturm-problem (nr-levels order)
  (let* ((domain (sturm-domain))
         (problem (sturm-problem domain))
         (mesh (sturm-mesh domain))
         (fe-class (lagrange-fe order))
         (as (make-fe-ansatz-space fe-class problem mesh)))
    (storing
      (solve (blackboard
              :ansatz-space as
              :success-if `(> :nr-levels ,nr-levels)
              :output t)))))

(defun install-electromagnetic-potential-demo ()
  (adjoin-demo
   (make-demo :name "Electromagnetic potential"
              :short "Electromagnetic potential"
              :long
              "Solves the Laplace equation for computing the
electromagnetic potential between two electrodes in a domain
with a piecewise constant permeability.  Triangle is used for
mesh generation.

This problem was suggested by Sebastian Sturm."
              :execute (_ (solve-sturm-problem 2 1)
                          (plot (getbb *result* :solution))))
   *laplace-demo*))

(install-electromagnetic-potential-demo)

;;;; Testing

(defun test-sturm ()

  (plot (periodic-polygonal 
         (list #d(0.0 0.0) #d(1.0 0.0) #d(1.0 1.0) #d(0.0 1.0))))

  (plot (triangulate (sturm-domain) :threshold 0.2 :meshsize 0.1))
  (describe (sturm-domain))
  (plot (sturm-mesh (sturm-domain)))
  
  (solve-sturm-problem 1 4)
  (plot (getbb *result* :solution))
  (plot (getbb *result* :solution) :program :vtk)
  (plot (getbb *result* :mesh) :program :vtk)
  
  )

;;; (test-sturm)
(fl.tests:adjoin-test 'test-sturm)
