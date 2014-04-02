;;; -*- mode: lisp -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; circle-ring-domain.lisp - Circle and cylinder
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

(in-package :fl.domains)

(file-documentation
 "Special domains including circle and cylinder.")

#+(or)
(defun get-angle (point)
  (cond
    ((plusp (aref point 0)) 0.0)
    ((minusp (aref point 0)) pi)
    ((plusp (aref point 1)) (* 0.5 pi))
    ((minusp (aref point 1)) (* 1.5 pi))
    ))

(defun linear-interpolate (s x y)
  (let ((s (vref s 0)))
    (+ (* (- 1 s) x)
       (* s y))))

(defun zap-to-unit-vectors (v)
  (cond ((mzerop (aref v 0) 1.0e-14) (double-vec 0.0 (signum (aref v 1))))
        ((mzerop (aref v 1) 1.0e-14) (double-vec (signum (aref v 0)) 0.0))
        (t v)))

(defun triangle-mapping (r phi1 phi2 &optional flipped-p)
  (let ((eval
         (lambda (v)
           (let* ((s (norm v 1))
                  (vn (scal (if (plusp s) (/ s) 1.0) v))
                  (phi (linear-interpolate (vref vn 1)
                                           (if (and (zerop phi1) (> phi2 pi))
                                               (* 2 pi)
                                               phi1)
                                           phi2)))
             (scal (* s r)
                   (zap-to-unit-vectors
                    (double-vec (if flipped-p (- (cos phi)) (cos phi))
                                (sin phi))))))
         #+(or)
         (lambda (v)
           (let ((phi (linear-interpolate
                       (* (/ 2 pi) (phase (complex (vref v 0) (vref v 1))))
                       (if (and (zerop phi1) (> phi2 pi)) (* 2 pi) phi1)
                       phi2))
                 (r (norm v 1)))
             (scal (* r r2)
                   (zap-to-unit-vectors
                    (double-vec (if flipped-p (- (cos phi)) (cos phi))
                                (sin phi))))))))
  (make-instance
   '<special-function>
   :domain-dimension 2 :image-dimension 2
   :evaluator eval :gradient (numerical-gradient eval))))

(defun polar->rect (r1 r2 phi1 phi2 &optional flipped-p delta)
  (let ((rk -1)
        (phik -1))
    (when (/= r1 r2) (incf rk))
    (when (/= phi1 phi2) (setf phik (1+ rk)))
    (let ((eval
           (lambda (v)
             (let* ((sr (if (minusp rk) 0.0 (vref v rk)))
                    (sphi (if (minusp phik)
                              (if (plusp r1) 0.0 1.0)
                              (vref v phik)))
                    (r (linear-interpolate sr r1 r2))
                    (special-p (and (zerop phi1) (> phi2 pi)))
                    (phi (if (minusp phik)
                             (if (plusp r1) phi1 phi2)
                             (linear-interpolate
                              sphi (if special-p (* 2 pi) phi1) phi2))))
               (when delta
                 (dbg :domains "sr=~A, sphi=~A" sr sphi)
                 (dbg :domains "r=~A, phi=~A" r phi)
                 (multiple-value-setq (r phi)
                   (let ((z (linear-interpolate
                             sr
                             (complex r1 0.0)
                             (complex #I"sqrt(r2^^2-delta^^2)" (if special-p (- delta) delta)))))
                     (dbg :domains "z=~A, |z|=~A, arg(z)=~A" z (abs z) (phase z))
                     (values
                      ;; should be small modifications of r and phi
                      (linear-interpolate sphi (abs z) r)
                      (+ phi
                         (linear-interpolate sphi (phase z) 0.0)))))
                 (dbg :domains "r=~A, phi=~A" r phi)
                 )
               (scal r (zap-to-unit-vectors
                        (double-vec (if flipped-p (- (cos phi)) (cos phi))
                                    (sin phi))))))))
      (make-instance
       '<special-function>
       :domain-dimension (- 2 (count-if #'minusp (list rk phik)))
       :image-dimension 2
       :evaluator eval :gradient (numerical-gradient eval)))))

(defun circle-ring-domain (r1 r2 &key (interior-p t) (channel-breadth 0.0) flipped-p)
  (declare (optimize debug))
  (let* ((pi/2 (* 0.5 pi))
         (thick-channel-p (plusp channel-breadth))
         (ring-p r2)
         (dphi (if ring-p
                   (asin (/ channel-breadth 2 r2))
                   0.0))
         (pc-nodes
          `((0.0 0.0)                                                 ; 0
            (,r1 0.0) (,r1 ,pi/2) (,r1 ,(* 2 pi/2)) (,r1 ,(* 3 pi/2)) ; 1-4
            ,(if thick-channel-p                                      ; 5
                 (list r2 dphi)
                 (list r2 0.0))
            (,r2 ,pi/2) (,r2 ,(* 2 pi/2)) (,r2 ,(* 3 pi/2))           ; 6-8
            ,@(when thick-channel-p                                   ;9
                    `((,r2 ,(- (* 2 pi) dphi))))
            ))
         (table (make-hash-table :test 'equalp))
         )
    (labels
        ((insert-cell (indices &optional new-p)
           (or (if new-p
                   (setf
                    (gethash indices table)
                    (case (length indices)
                      (1 (destructuring-bind (r phi)
                             (elt pc-nodes (first indices))
                           (make-vertex (evaluate (polar->rect r r phi phi flipped-p) #()))))
                      (2 (destructuring-bind (k1 k2) indices
                           (make-line (insert-cell (list k1))
                                      (insert-cell (list k2))
                                      :mapping
                                      (destructuring-bind ((r1 phi1) (r2 phi2))
                                          (list (elt pc-nodes k1) (elt pc-nodes k2))
                                      (unless (or (/= r1 r2)
                                                  (and thick-channel-p
                                                       (= k1 5) (= k2 9)))
                                        (polar->rect r1 r2 (if (zerop r1) phi2 phi1) phi2 flipped-p))))))
                      (3 (make-simplex
                          (map 'vector (_ (insert-cell (remove _ indices))) indices)
                          :mapping
                          ;; we use the special shape of the triangle
                          (destructuring-bind (k1 k2 k3) indices
                            (unless (plusp k1)
                              (destructuring-bind ((r2 phi2) (r3 phi3))
                                  (list (elt pc-nodes k2) (elt pc-nodes k3))
                                (declare (ignore r3))
                                (triangle-mapping r2 phi2 phi3 flipped-p))))))
                      (4 (destructuring-bind (k1 k2 k3 k4) indices
                           (make-instance
                            (fl.mesh::product-cell-class '(1 1) t)
                            :boundary 
                            (vector (insert-cell (list k2 k3))
                                    (insert-cell (list k1 k4))
                                    (insert-cell (list k4 k3))
                                    (insert-cell (list k1 k2)))
                            :mapping
                            (destructuring-bind ((r1 phi1) (r2 phi2))
                                (list (elt pc-nodes k1) (elt pc-nodes k3))
                              (polar->rect r1 r2 phi1 phi2 flipped-p
                                           (and thick-channel-p (= k1 1)
                                                (* 0.5 channel-breadth)))))))
                      ))
                   (gethash indices table))
               (error "should not happen"))))
      (loop for indices in
           (append
            (when (or interior-p ring-p)
              '((1) (2) (3) (4)
                (1 2) (2 3) (3 4) (1 4)))
            (when interior-p
              '((0) (0 1) (0 2) (0 3) (0 4)
                (0 1 2) (0 2 3) (0 3 4) (0 1 4)))
            (when ring-p
              `((5) (6) (7) (8)
                ;; index 9 should be considered more like 5.5
                ,@(when thick-channel-p '((9)))
                (1 5) (2 6) (3 7) (4 8)
                ,@(when thick-channel-p '((1 9)))
                (5 6) (6 7) (7 8)
                ,@(if thick-channel-p
                      '((9 8) (5 9))
                      '((5 8)))
                (1 5 6 2) (2 6 7 3) (3 7 8 4)
                ,@(if thick-channel-p
                      '((1 9 8 4) (1 5 9))
                      '((1 5 8 4))))))
         do (insert-cell indices t))
      (make-domain (make-instance '<skeleton>
                                    :dimension 2
                                    :cells (hash-table-values table))))))

(defun cylinder-domain (r h &key flipped-p)
  (make-domain
   (fl.mesh::cartesian-product (circle-ring-domain r nil :flipped-p flipped-p)
                               (box-domain `((0.0 ,h))))))

(defun flat-disk-domain (r1 r2 h &rest paras)
  (make-domain
   (fl.mesh::cartesian-product (apply #'circle-ring-domain r1 r2 paras)
                               (box-domain `((,(- h) 0.0))))))

;;;; Testing

(defun test-circle-ring-domain ()
  (fl.plot:plot (uniformly-refined-mesh (flat-disk-domain 1.0 1.5 .3) 3
                                        :parametric :from-domain))
  (cylinder-domain 1.0 1.0)

  (vref 0.0 0)
  
  (evaluate (polar->rect 0.0 1.0 0.0 0.0) #d(1.0))

  (let ((domain
         (?2 (circle-ring-domain 1.0 1.5 :interior-p t :channel-breadth 0.0 :flipped-p nil)
             (circle-ring-domain 1.0 1.5 :interior-p t :channel-breadth 0.4 :flipped-p nil)
             (circle-ring-domain 1.0 2.0 :interior-p t)
             (cylinder-domain 1.0 1.0))))
    (check domain)
    (fl.plot:plot
     (uniformly-refined-hierarchical-mesh domain 2 :parametric :from-domain)
     :program :dx))
  (check (circle-ring-domain 1 1.5 :channel-breadth 0.1 :flipped-p t))

  (let* ((domain (circle-ring-domain 1 1.5 :channel-breadth 0.1))
         (mesh (uniformly-refined-mesh domain 0 :parametric :from-domain)))
    (check mesh))
  )

(fl.tests:adjoin-test 'test-circle-ring-domain)
             