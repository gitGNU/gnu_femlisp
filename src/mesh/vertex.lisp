;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; vertex.lisp
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
;;;; <vertex>
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; The class definition is contained in cell.lisp

(defun print-vertex (vtx stream depth)
  (declare (ignore depth))
  (print-unreadable-object
   (vtx stream :type t :identity t)
   (princ (vertex-position vtx) stream)))

(defun reference-vertex ()
  (find-reference-cell-with #'(lambda (refcell) (zerop (dimension refcell)))))

(defun make-vertex (position)
  "The vertex constructor."
  (make-<vertex>
   :cell-class *vertex-class*
   :mapping (if (typep position 'double-vec)
		position
		(map 'double-vec #'identity position))))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Routines
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun vertex? (cell) (typep cell '<vertex>))

(defmethod vertices ((vtx <vertex>)) (list vtx))

(defmethod l2g ((vtx <vertex>) local-pos)
  "Same as local->global for vertices."
  (assert (equalp local-pos (double-vec)))
  (vertex-position vtx))

(defmethod manifold-dimension ((vtx <vertex>))
  "Anchor for recursive definition."
  (let ((mapping (mapping vtx)))
    (if (vectorp mapping)
	(length mapping)
	(image-dimension mapping))))

(defmethod l2Dg ((vtx <vertex>) local-pos)
  (assert (equalp local-pos (double-vec)))
  (make-float-matrix (manifold-dimension vtx) 0))

(defmethod local->Dglobal ((vtx <vertex>) local-pos)
  "Not perfect, should take mapping into account."
  (l2Dg vtx local-pos))

(defmethod l2jet ((vtx <vertex>) (local-pos array) (k integer))
  "The jet of a simplex with linear cell mapping is 0 above the first
derivative."
  (loop for i from 0 upto k
	collect
	(case i
	  ((0) (l2g vtx local-pos))
	  ((1) (l2Dg vtx local-pos))
	  (t (make-real-tensor (cons (manifold-dimension vtx)
				     (make-list (1- i) :initial-element 0)))))))

(defmethod g2l ((vtx <vertex>) global-pos)
  (when (equalp global-pos (vertex-position vtx))
    (double-vec)))

(defmethod global->local ((vtx <vertex>) global-pos)
  (when (equalp global-pos (evaluate (mapping vtx) (double-vec)))
    (double-vec)))

(defmethod coordinates-inside? ((vtx <vertex>) local-pos)
  (equalp local-pos (double-vec)))

(defmethod local-coordinates-of-midpoint ((cell <vertex>))
  (double-vec))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; refinement
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmethod primary-refine-info ((refcell <vertex>))
  "Construction of the refinement information." 
  (vector
   (make-<child-info>
    :class *vertex-class*
    :barycentric-corners
    (list (list (make-double-vec 1 1.0d0)))  ; not perfect, because no real factor
    :boundary-paths ())))

(defmethod update-refine-info! ((refcell <vertex>)) nil)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Activation of the vertex class
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun ensure-vertex ()
  (or (reference-vertex)
      (initialize-cell-class *reference-vertex*)
      *reference-vertex*))

(ensure-vertex)


;;;; Testing: (test-vertex)
(defun test-vertex ()
  (assert (= 0 (manifold-dimension *reference-vertex*)))
  (assert (eq *reference-vertex* (reference-cell *reference-vertex*)))
  (reference-cell? *reference-vertex*)
    (describe (etable (skeleton *reference-vertex*) 0))
  (describe *reference-vertex*)
  (skeleton-boundary (skeleton *reference-vertex*))
  ;;
  (refine-info *reference-vertex*)
  (describe (refcell-skeleton *reference-vertex*))
  (describe (refcell-refinement-skeleton *reference-vertex* 1))
  (describe (refcell-skeleton *reference-vertex*))
  (refcell-refinement-vertices *reference-vertex* 2)
  )

(tests::adjoin-femlisp-test 'test-vertex)

