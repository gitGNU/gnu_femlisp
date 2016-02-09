;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; blending.lisp - Introduces blending functions
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;
;;; Copyright (C) 2010 Nicolas Neuss, KIT Karlsruhe
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
;;; MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
;;; IN NO EVENT SHALL THE AUTHOR, THE KARLSRUHE INSTITUTE OF TECHNOLOGY OR
;;; OTHER CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
;;; SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
;;; LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
;;; DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
;;; THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
;;; (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
;;; OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(in-package :fl.mesh)

(file-documentation
 "Here we introduce blending maps for simplices, which interpolate between
 given maps on the simplex boundary.  The central idea is to construct for
 each vertex of the simplex a function interpolating linearly between the
 faces of those boundaries not opposite to that vertex.  These n functions
 are then linearly combined with the weight of the standard P^1 basis
 functions.  The resulting map should be globally C^0 and in the simplex as
 good as the boundary functions permit.")

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Basic
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; We use from cell.lisp, simplex.lisp the following functions:
;;;   euclidean->barycentric (pos)
;;;   barycentric-coordinates ((refcell <simplex>) local-pos)
;;;   barycentric-gradients ((cell <simplex>) local-pos)

(defun barycentric->euclidean (s)
  (funcall (l2g-evaluator (n-simplex (1- (length s)))) s))

(defun without-component (s k)
  (concatenate (if (listp s) 'list 'double-vec)
               (subseq s 0 k) (subseq s (1+ k))))

(defun sub-side-mapping (k sk j dim)
  "Returns the side mapping j for the intersection of the simplex with a
level surface sk=const.  @arg{dim} is the dimension of the original
simplex."
  (lambda (s1)
    (assert (> dim 1))
    (assert (= (1- dim) (length s1)))
    (assert (and (<= 0 k dim) (<= 0 j dim) (/= k j)))
    (lret ((result
            (loop for i from 0 upto dim
               do
                 (dbg-show :blending k j i s1)
               unless (= i j) collect
               (if (= i k)
                   sk
                   (* (pop s1) (- 1.0 sk))))))
      (assert (= (length result) dim)))))

(defun ip (mappings s)
  "mappings contains n+1 mappings from barycentric side coordinates of
length n, s contains n+1 coordinates 0.0<s[i]<1.0"
  (declare (optimize debug))
  (assert (= (length mappings) (length s)))
  ;; handle the trivial case that one coordinate is zero
  (loop+ ((si s) (map mappings) k)
       when (zerop si) do
       (return-from ip (funcall map (without-component s k))))
  ;; now for the interpolating case
  (let ((dim (1- (length mappings))))
    (if (zerop dim)
        ;; we must evaluate
        (funcall (first mappings) s)
        ;; we recurse on the lower-dimensional sides
        (lret (result)
          (loop for k from 0 upto dim do
               (let* ((sk (elt s k))
                      (summand
                       (scal sk         ; with weight sk
                             (if (= dim 1)
                                 (funcall (elt mappings (if (= k 0) 1 0)) '(1.0))
                                 (ip (mapcar (lambda (j)
                                               (compose (elt mappings j)
                                                        (sub-side-mapping k sk j dim)))
                                             (remove k (range<= 0 dim)))
                                     (scal (/ (- 1.0 sk)) (without-component s k)))))))
                 (if result
                     (m+! summand result)
                     (setf result summand))))))))

(defun interpolate (mappings)
  "Interpolates a vector of functions f of length n+1 on the reference
n-simplex.  The functions are understood as providing functional values for
the n+1 n-1-dimensional boundary simplices."
  (let ((mappings
         (map 'list
              (lambda (f)
                (_ (evaluate f (barycentric->euclidean (coerce _ 'double-vec)))))
              mappings)))
    (lambda (x)
      (ip mappings (euclidean->barycentric x)))))

(defun blending-map (cell)
  "Returns a mapping for the cell by interpolating the mappings of its
boundary."
  (assert (simplex-p cell))
  (let ((eval (interpolate (map 'list #'cell-mapping (boundary cell)))))
    (make-instance
     '<special-function>
     :domain-dimension (dimension cell)
     :image-dimension (embedded-dimension cell)
     :evaluator eval
     :gradient (numerical-gradient eval))))

(defun change-to-blending (cell)
  (change-class cell (mapped-cell-class (class-of cell))
                :mapping (blending-map cell)))

(defun curved-triangle-domain ()
  "A quarter of a unit disk."
  (lret ((domain
          (make-domain (copy-skeleton (skeleton (n-simplex 2))))))
    (let* ((triangle (first (cells-of-highest-dim domain)))
           (line (elt (boundary triangle) 0)))
      (change-class line 'fl.mesh::<mapped-1-simplex>
                    :mapping (circle-function 1.0 #d(0.0 0.0) (* 0.5 pi)))
      (change-to-blending triangle))
    (check domain)
    ))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Tests
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(in-suite mesh-suite)

(defun test-f (cell f)
  "Returns a mapping for cell transformed with f"
  (lambda (x)
    (evaluate f (l2g cell x))))

(test blending

  (is (equalp #d(0.6 0.4) (euclidean->barycentric #d(0.4))))
  (let* ((dim 2)
         (cell (n-simplex dim))
         (A #m((2.0 -1.0) (-0.5 3.0)))
         (f (curry #'m* A))
         (bdry-f (vector-map (_ (test-f _ f))
                             (boundary cell)))
         (ip (interpolate bdry-f))
         (x #d(0.4 0.5)))
    (is-true (mzerop (m- (m* A x) (funcall ip x)) 1.0e-15)))
  
  (finishes
    ;; (dbg-on :blending)
    ;; (dbg-off :blending)
    (loop for (dim pos) in 
          '((1 (0.6 0.4)) (2 (0.3 0.4 0.3)) (3 (0.2 0.2 0.3 0.3)))
          do
             (let* ((boundary
                      (map 'list
                           (lambda (side)
                             (_ (evaluate (cell-mapping side)
                                          (barycentric->euclidean (coerce _ 'double-vec)))))
                           (boundary (n-simplex dim)))))
               (loop repeat 10 do (ip boundary pos))))
    (let* ((domain (curved-triangle-domain))
           (mesh (uniformly-refined-mesh domain 4 :parametric :from-domain)))
      ;; (fl.plot:plot mesh)
      (check mesh)))

  )
