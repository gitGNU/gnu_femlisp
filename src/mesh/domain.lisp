;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;  domain.lisp
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

(in-package :mesh)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Domain class
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass <domain> (<skeleton>)
  ((boundary :accessor domain-boundary)
   (extensible-p :reader extensible-p :initform nil))
  (:documentation "A <domain> is a special <skeleton>.  Its cells are
called patches, and the values are property lists carrying geometric
information, e.g. metric, volume-form, embedding or identification."))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Patches
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defparameter *known-patch-properties* '(METRIC VOLUME IMBEDDING EXTENSION)
  "This list contains some pre-defined keywords for use as patch
properties.")

(defun patch-identification (patch domain)
  "Returns a list of identified cells."
  (getf (skel-ref domain patch) 'IDENTIFIED))
(defun (setf patch-identification) (patch-list patch domain)
  "Sets the identification of patch to a list of identified patches."
  (setf (getf (skel-ref domain patch) 'IDENTIFIED) patch-list))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Domain generation
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun ensure-secondary-information (domain)
  "Preliminary.  Maybe it would be better to update the boundary
automatically when inserting cells.  Unfortunately, it is not clear
how this could be done in an easy way."
  (setf (domain-boundary domain) (skeleton-boundary domain))
  (doskel (patch domain)
    (when (get-cell-property patch domain 'EXTENSION)
      (setf (slot-value domain 'extensible-p) t))))
   
(defmethod initialize-instance :after ((domain <domain>) &key cells &allow-other-keys)
  "When a domain is constructed from a list of cells, we assume that its
definition is finished and setup the boundary slot."
  (when cells (ensure-secondary-information domain)))

(defmethod update-instance-for-different-class :after
    ((skel <skeleton>) (domain <domain>) &rest initargs)
  "Sometimes a skeleton is transformed into a domain by change-class.
Usually, this means that the definition is finished such that we can
compute the boundary afterwards."
  (declare (ignore initargs))
  (ensure-secondary-information domain))

(defun domain-characteristics (domain)
  "Returns a property list of characteristics.  The property curved means
that curved patches exist.  The property exact is set to t if the domain
mappings are exact.  Otherwise, only the boundary of the domain should be
assumed to be provided in an exact form."
  (let ((properties (list :curved nil :exact t)))
    (doskel (cell domain)
      (when (and (not (vertex? cell)) (mapping cell))
	(setf (getf properties :curved) t)))
    (doskel (cell domain)
      (unless (or (vertex? cell) (aand (mapping cell) (differentiable-p it)))
	(loop for side across (boundary cell)
	      when (and (mapping side) (not (vertex? side))) do
	      (setf (getf properties :exact) nil))))
    properties))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Frequently used domains
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; n-simplex
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun n-simplex-domain (dim)
  (change-class (skeleton (n-simplex dim)) '<domain>))

(defparameter *unit-interval-domain* (n-simplex-domain 1))
(defparameter *unit-triangle-domain* (n-simplex-domain 2))
(defparameter *unit-tetrahedron-domain* (n-simplex-domain 3))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; n-cube
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun n-cube-domain (dim)
  (change-class (skeleton (n-cube dim)) '<domain>))
(defparameter *unit-quadrangle-domain* (n-cube-domain 2))
(defparameter *unit-cube-domain* (n-cube-domain 3))

(defun box-domain (dimensions)
  "Generates a box domain for the given dimensions.  Here,
dimensions is expected to be a list of 2-element lists denoting
the interval along the respective axis.  The algorithm works by
copying the unit cube and modifying the coordinates of the
vertices of the copy."
  (let ((new-skel (copy-skeleton (skeleton (n-cube (length dimensions))))))
    (doskel (vertex new-skel :dimension 0)
      (loop with pos = (vertex-position vertex)
	    for dims in dimensions
	    and i from 0 do
	    (setf (aref pos i)
		  (coerce (elt dims (if (zerop (aref pos i)) 0 1))
			  'double-float))))
    (change-class new-skel '<domain>)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; n-cell
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun n-cell-domain (dim)
  "Generates an n-dimensional cell domain which is a n-dimensional unit
cube with its opposite sides identified."
  (let* ((domain (n-cube-domain dim))
	 (cube (car (cells-of-highest-dim domain))))
    (mapc (rcurry #'identify domain)
	  (iterate-identifications
	   (loop with bdry = (boundary cube)
		 for i from 0 below (length bdry) by 2 collect
		 (list (aref bdry i) (aref bdry (1+ i))))))
    (ensure-secondary-information domain)
    domain))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; unit ball
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun n-ball-phi (x)
  (let ((norm2 (norm x 2)))
    (if (zerop norm2)
	x
	(scal (/ (norm x 1) norm2) x))))

(defun n-ball-Dphi (x)
  (let ((dim (length x)))
    (if (zerop (norm x 2))
	(eye dim)
	(loop with 1/norm2 = (/ (norm x 2))
	      with alpha = (* (norm x 1) 1/norm2)
	      with result = (make-real-matrix dim)
	      for i from 0 below dim do
	      (loop with factor = (- (signum (aref x i)) (* alpha (aref x i) 1/norm2))
		    for k from 0 below dim do
		    (setf (mref result k i)
			  (* (aref x k) 1/norm2 factor))
		    (when (= i k)
		      (incf (mref result k i) alpha)))
	      finally (return result)))))

(defun n-ball-domain (dim)
  "Generates an n-dimensional ball domain with 2^n simplex patches."
  (let* ((cell-list (make-hash-table :test #'equalp))
	 (domain (make-instance '<domain> :dimension dim)))
    ;; generate and intern the cells into domain
    (flet ((insert! (corners)
	     (let ((cell
		    (cond ((single? corners)
			   (let ((id (car corners)))
			     (cond ((zerop id) (make-vertex (make-double-vec dim)))
				   ((plusp id) (make-vertex (unit-vector dim (1- id))))
				   (t (make-vertex (scal -1.0 (unit-vector dim (- -1 id))))))))
			  (t (make-simplex
			      (map 'cell-vec #'(lambda (side-ids) (gethash side-ids cell-list))
				   (mapcar #'(lambda (corner) (remove corner corners)) corners)))))))
	       (unless (vertex? cell)
		 (change-class
		  cell (mapped-cell-class (class-of cell))
		  :mapping (make-instance
			    '<special-function>
			    :domain-dimension (1- (length corners))
			    :image-dimension dim
			    :evaluator #'(lambda (x) (n-ball-phi (l2g cell x)))
			    :gradient #'(lambda (x) (m* (n-ball-Dphi (l2g cell x))
							(l2Dg cell x))))))
	       (setf (gethash corners cell-list) cell)
	       (insert-cell! domain cell))))
    
      (loop initially (insert! '(0))
	    for k from 1 upto (1+ dim)
	    for sign-lists = (apply #'map-product #'list (make-list k :initial-element '(-1 1))) do
	    (loop for k-subset in (k-subsets (range 1 dim) k) do
		  (loop for signs in sign-lists
			for corners = (mapcar #'* k-subset signs) do
			(insert! corners)
			(insert! (cons 0 corners))))))
    ;; finally ensure the boundary
    (ensure-secondary-information domain)
    domain))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Special domains (mainly for testing purposes)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defparameter *circle-domain*
  (let* ((pi/2 (* 0.5 pi))
	 (pi/2-scaling (make-real-matrix `((,pi/2))))
	 (circle-boundary
	  (make-instance
	   '<special-function> :domain-dimension 1 :image-dimension 2
	   :evaluator #'(lambda (phi)
			  (double-vec (cos (vref phi 0)) (sin (vref phi 0)))))))
    ;; corners
    (let ((center (make-vertex #d(0.0 0.0)))
	  (east-vtx (make-vertex #d(1.0 0.0)))
	  (north-vtx (make-vertex #d(0.0 1.0)))
	  (west-vtx (make-vertex #d(-1.0 0.0)))
	  (south-vtx (make-vertex #d(0.0 -1.0))))
      ;; inner line segments
      (let ((seg-ce (make-line center east-vtx))
	    (seg-cn (make-line center north-vtx))
	    (seg-cw (make-line center west-vtx))
	    (seg-cs (make-line center south-vtx)))
	;; curved boundaries
	(let ((seg-en
	       (make-line
		east-vtx north-vtx
		:mapping (transform-function
			  circle-boundary :domain-transform
			  (list pi/2-scaling #d(0.0)))))
	      (seg-nw
	       (make-line
		north-vtx west-vtx
		:mapping (transform-function
			  circle-boundary :domain-transform
			  (list pi/2-scaling (double-vec pi/2)))))
	      (seg-ws
	       (make-line
		west-vtx south-vtx
		:mapping (transform-function
			  circle-boundary :domain-transform
			  (list pi/2-scaling (double-vec pi)))))
	      (seg-es
	       (make-line
		east-vtx south-vtx
		:mapping (transform-function
			  circle-boundary :domain-transform
			  (list pi/2-scaling (double-vec (* 2 pi)))))))
	  ;; Now the four triangle cells.  Note, that we don't bother about
	  ;; the precise mappings for now.  Nevertheless, these would be needed,
	  ;; if we wanted to work with completely nonlinear (and not only
	  ;; isoparametric) cell mappings.
	  (let ((tri-1 (make-simplex (vector seg-en seg-cn seg-ce)))
		(tri-2 (make-simplex (vector seg-nw seg-cw seg-cn)))
		(tri-3 (make-simplex (vector seg-ws seg-cs seg-cw)))
		(tri-4 (make-simplex (vector seg-es seg-cs seg-ce))))
	    ;; Finally, construct the domain
	    (make-instance '<domain> :cells (list tri-1 tri-2 tri-3 tri-4)))))))
  "This definition of a circle domain gives somewhat better results than
the general n-ball-domain, probably because the boundary parametrization is
better.")

(defparameter *rotated-square-domain*
  ;; corners
  (let ((center (make-vertex #d(0.0 0.0)))
	(east-vtx (make-vertex #d(1.0 0.0)))
	(north-vtx (make-vertex #d(0.0 1.0)))
	(west-vtx (make-vertex #d(-1.0 0.0)))
	(south-vtx (make-vertex #d(0.0 -1.0))))
    ;; inner line segments
    (let ((seg-ce (make-line center east-vtx))
	  (seg-cn (make-line center north-vtx))
	  (seg-cw (make-line center west-vtx))
	  (seg-cs (make-line center south-vtx)))
      ;; boundary segments
      (let ((seg-en (make-line east-vtx north-vtx))
	    (seg-nw (make-line north-vtx west-vtx))
	    (seg-ws (make-line west-vtx south-vtx))
	    (seg-es (make-line east-vtx south-vtx)))
	;; finally the four triangle cells
	(let ((tri-1 (make-simplex (vector seg-en seg-cn seg-ce)))
	      (tri-2 (make-simplex (vector seg-nw seg-cw seg-cn)))
	      (tri-3 (make-simplex (vector seg-ws seg-cs seg-cw)))
	      (tri-4 (make-simplex (vector seg-es seg-cs seg-ce))))
	  ;; then construct the domain
	  (make-instance '<domain> :cells (list tri-1 tri-2 tri-3 tri-4)))))))

(defun triangle-domain (corner1 corner2 corner3)
  (let ((vtx1 (make-vertex corner1))
	(vtx2 (make-vertex corner2))
	(vtx3 (make-vertex corner3)))
    (let ((seg1 (make-line vtx2 vtx3))
	  (seg2 (make-line vtx1 vtx3))
	  (seg3 (make-line vtx1 vtx2)))
      (make-instance '<domain> :cells (list (make-simplex (vector seg1 seg2 seg3)))))))

(defun L-domain (dim)
  "Creates an L-domain by cutting out a small cube of the uniform refinement of
the unit cube."
  (let* ((skel (refine-globally (skeleton (n-cube dim))))
	 (upper-right-cell
	  (find-cell-from-position skel (make-double-vec dim 0.75))))
    (change-class (skeleton-without-cell skel upper-right-cell)
		  '<domain>)))


;;;; Testing
(defun test-domain ()
  (display-ht (etable (n-ball-domain 2) 2))
  (corners *unit-triangle*)
  (let ((*print-skeleton-values* t))
    (describe (triangle-domain #d(0.0 0.0) #d(1.0 0.0) #d(0.0 1.0))))
  (check *unit-quadrangle-domain*)
  (check (n-cube-domain 2))
  (assert (= -1 (dimension (skeleton-boundary (n-cell-domain 2)))))
  )

;;; (mesh::test-domain)
(fl.tests:adjoin-test 'test-domain)

