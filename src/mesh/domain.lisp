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

(in-package :fl.mesh)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Domain class
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass <domain> (<skeleton> property-mixin)
  ((classifiers :initform () :initarg :classifiers :documentation
                "A list of functions of two arguments -patch and classifications so
far- which are called from the right to classify the patch."))
  (:documentation "A @class{<domain>} is a special @class{<skeleton>}.  We
call its cells @emph{patches}, and the properties of a patch carries
geometric information.  Properties supported up to now are:

@itemize
@item @code{IDENTIFIED}: @emph{list of identified patches}
@item @code{EXTENSION}: @emph{extender}
@item @code{METRIC}: @emph{metric tensor function}
@item @code{VOLUME}: @emph{volume function}
@end itemize

Metric and volume should be functions depending on keyword arguments like
@code{:LOCAL} and @code{:GLOBAL} and allowing arbitrary other keys."))

(defun domain-boundary (domain)
  "The boundary of the domain as a skeleton."
  (get-property domain 'boundary))

(defun domain-substance (domain)
  "The substance of the domain as a hash-table of cells."
  (get-property domain 'substance))

(defun domain-substance-boundaries (domain)
  "The boundaries of the substance cells as a hash-table of cells."
  (get-property domain 'substance-boundaries))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Domain generation
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun standard-classifier (domain)
  (lambda (patch classifications)
    (declare (ignore classifications))
    (let ((result (if (member-of-skeleton? patch (domain-boundary domain))
                      (list :skeleton-boundary :boundary :external-boundary)
                      (list :skeleton-interior :interior)))
          (dim (dimension patch)))
      (when (gethash patch (domain-substance domain))
        (_f union result '(:substance :d-dimensional)))
      (when (gethash patch (domain-substance-boundaries domain))
        (_f union result '(:d-1-dimensional)))
      (when (= dim 0)
        (when (mzerop (vertex-position patch))
          (push :origin result)))
      (when (= dim 1) (push :1-dimensional result))
      (when (= dim 0) (push :0-dimensional result))
      (awhen (get-cell-property patch domain :part)
        (pushnew it result))
      (awhen (get-cell-property patch domain 'classification)
        (setf result (append result it)))
      result)))

(defun ensure-secondary-information (domain)
  "Preliminary.  Maybe it would be better to update the boundary
automatically when inserting cells.  Unfortunately, it is not clear
how this could be done in an easy way."
  (with-slots (classifiers boundary substance substance-boundaries extensible-p) domain
    (setf boundary (skeleton-boundary domain))
    (setf substance (skeleton-substance domain))
    (setf substance-boundaries
          (lret ((result (make-hash-table)))
            (loop for cell being each hash-key of substance do
                 (loop for side across (boundary cell) do
                      (setf (gethash side result) t)))))
    (doskel (patch domain)
      (when (get-cell-property patch domain 'EXTENSION)
	(setf extensible-p t)))
    (_f append classifiers (list (standard-classifier domain)))))

(defmethod initialize-instance :after ((domain <domain>) &key &allow-other-keys)
  "When a domain is constructed from a list of cells, we assume that its
definition is finished and setup the boundary slot."
  (ensure-secondary-information domain))

(defmethod update-instance-for-different-class :after
    ((skel <skeleton>) (domain <domain>) &rest initargs)
  "Sometimes a skeleton is transformed into a domain by change-class.
Usually, this means that the definition is finished such that we can
compute the boundary afterwards."
  (declare (ignore initargs))
  (ensure-secondary-information domain))

(defun patch-classification (patch domain)
  "Returns a list of classifications for @arg{patch} in @arg{domain}."
  (labels ((classify (classifiers)
	     (and classifiers
		  (funcall (car classifiers) patch (classify (cdr classifiers))))))
    (classify (slot-value domain 'classifiers))))

(defun make-domain (skel &rest args)
  "Changes a skeleton into a domain."
  (assert (typep skel '<skeleton>))
  (apply #'change-class skel '<domain> args))

#+(or)
(defun test-condition (condition classifications)
  "Test if @arg{condition} which at the moment should be either an AND or
an OR combination of symbols in the list @arg{classification} holds."
  (cond
    ((eq t condition))
    ((symbolp condition) (member condition classifications))
    ((consp condition)
     (funcall (ecase (car condition)
		(and #'subsetp)
		(or #'intersection))
	      (cdr condition)
	      classifications))
    (t (error "Unknown patch classification scheme."))))

(defun test-condition (condition classifications)
  "Test if @arg{condition} holds for @arg{classifications}.
@arg{condition} should be a logical combination of the keyword symbols in
the list @arg{classification}."
   (funcall (eval `(lambda (classifications)
		     (declare (ignorable classifications))
		     ,(map-tree (lambda (leaf)
				  (if (keywordp leaf)
				      (find leaf classifications)
				      leaf))
				condition)))
	    classifications))

(defun find-patch (domain classifiers &key midpoint)
  "Searches for a patch of @arg{domain} having the given list
@arg{classifiers}."
  (find-cell (lambda (patch)
               (and (subsetp classifiers (patch-classification patch domain))
                    (or (null midpoint)
                        (equalp midpoint (midpoint patch)))))
             domain))

(defun make-classifier (test classification)
  "Returns a classifier for patches."
   (lambda (patch classifiers)
     (if (funcall test patch)
         (adjoin classification classifiers)
         classifiers)))

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

(defmethod describe-object :after ((domain <domain>) stream)
  (format stream "~&Domain characteristics: ~S~%" (domain-characteristics domain))
  (format stream "~&Patch characterizations:~%")
  (doskel (patch domain)
    (format t "~A ->~% ~S~%" patch (patch-classification patch domain))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Frequently used domains
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; n-simplex
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun n-simplex-domain (dim)
  (make-domain (skeleton (n-simplex dim))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; product-cells, n-cubes
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun simplex-product-domain (dims)
  "Generates a product-cell domain for the given factor dimensions."
  (make-domain (skeleton (ensure-simplex-product dims))))

(defun classify-top-bottom-lateral (&optional (bottom 0.0) (top 1.0) (left 0.0) (right 1.0))
  "Returns a classifier for identifying top, bottom, and lateral parts of a
cube domain."
  (lambda (cell classifications)
    (let* ((midpoint (midpoint cell))
	   (dim (length midpoint))
           (cell-dim (dimension cell)))
      (when (and (= cell-dim 0) (mzerop midpoint))
        (pushnew :origin classifications))
      (when (< cell-dim dim)
	(let ((coord (aref midpoint (1- dim))))
	  (pushnew (cond ((= coord top) :top)
			 ((= coord bottom) :bottom)
			 (t :lateral))
		   classifications))
	(when (< cell-dim (1- dim))
	  (pushnew :lateral classifications))
	;; the first coordinate classifies left and right
        (cond ((= (elt midpoint 0) left) (pushnew :left classifications))
              ((= (elt midpoint 0) right) (pushnew :right classifications)))))
    classifications))

(defun n-cube-domain (dim)
  (make-domain (skeleton (n-cube dim))
               :classifiers (list (classify-top-bottom-lateral)
                                  (lambda (cell classifications)
                                    (when (mzerop (midpoint cell))
                                      (pushnew :origin classifications))
                                    classifications))))

(defun ensure-domain (domain)
  "If @arg{domain} is an integer, return the corresponding
@arg{n-cube-domain}, if @arg{domain} is a domain return it unchanged,
otherwise signal an error."
  (if (integerp domain)
      (n-cube-domain domain)
      (if (typep domain '<domain>)
	  domain
	  (error "Expected either a dimension or a domain."))))

(defun box-domain (dimensions)
  "Generates a box domain for the given dimensions.  Here,
dimensions is expected to be a list of 2-element lists denoting
the interval along the respective axis.  The algorithm works by
copying the unit cube and modifying the coordinates of the
vertices of the copy."
  (let* ((dim (length dimensions))
	 (new-skel (copy-skeleton (skeleton (n-cube dim)))))
    (doskel (vertex new-skel :dimension 0)
      (loop with pos = (vertex-position vertex)
	    for dims in dimensions
	    and i from 0 do
	    (setf (aref pos i)
		  (coerce (elt dims (if (zerop (aref pos i)) 0 1))
			  'double-float))))
    (make-domain new-skel)))


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
			      (map-cell-vec #'(lambda (side-ids)
                                                (the <cell> (gethash side-ids cell-list)))
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
	    (loop for k-subset in (k-subsets (range<= 1 dim) k) do
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

(defun cube-without-corner (dim)
  "Creates a domain by cutting out one cube from the uniform refinement of
the unit cube."
  (let* ((skel (refine (skeleton (n-cube dim))))
	 (upper-right-cell
	  (find-cell-from-position skel (make-double-vec dim 0.75))))
    (make-domain (skeleton-without-cell skel upper-right-cell))))

(defun L-domain (dim)
  "Creates an L-domain by cutting out cubes from the uniform refinement of
the unit cube."
  (make-domain
   (skeleton (remove-if
	      (lambda (cell)
		(let ((midpoint (midpoint cell)))
		  (and (> (aref midpoint 0) 0.5)
		       (> (aref midpoint 1) 0.5))))
	      (cells-of-highest-dim
	       (refine (skeleton (n-cube dim))))))))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Tests
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(in-suite mesh-suite)

(test domain
  (is (= -1 (dimension (skeleton-boundary (n-cell-domain 2)))))
  (let* ((domain (n-cube-domain 3))
	 (cell (find-cell (lambda (c) (mequalp (midpoint c) #d(0.5 0.5 0.0))) domain)))
    (is-true (member :bottom (patch-classification cell domain))))
  (is-false (test-condition
             '(and (or :left :right) :d-1-dimensional)
             '(:top :d-1-dimensional :bottom)))
  (is-false (test-condition '(and :d-dimensional (not :inlay))
                            '(:d-dimensional :inlay)))
  (finishes
    (doskel (cell *rotated-square-domain* :direction :down)
      (print cell))
    (describe (n-cube-domain 2))
    (corners (n-simplex 2))
    (let ((*print-skeleton-values* t))
      (describe (triangle-domain #d(0.0 0.0) #d(1.0 0.0) #d(0.0 1.0))))
    (check (n-cube-domain 2))
    (let ((domain (n-cube-domain 1)))
      (doskel (patch domain)
        (format t "~A : ~S~%" patch (patch-classification patch domain))))
    (let ((domain (n-cube-domain 3)))
      (doskel (cell domain)
        (print (patch-classification cell domain))))
    (let ((domain (n-cube-domain 1)))
      (doskel (patch domain)
        (print (funcall (standard-classifier domain) patch ()))))
    )
  )

