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
;;;; Vertex class and reference vertex generation
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; The class <vertex> has already been defined in cell.lisp.  For
;;; completeness, we define also mapped vertices.

(defclass <mapped-vertex> (<mapped-cell> <vertex>)
  ()
  (:documentation "At the moment, this class is not needed because vertices
are not mapped."))

(defmethod initialize-instance :after ((cell <mapped-vertex>) &key &allow-other-keys)
  (error "Up to now, the <mapped-vertex> should not be used."))

(let ((empty-cell-vec #()))
  (defmethod boundary ((vtx <vertex>))
    "Returns the empty boundary for vertices."
    empty-cell-vec))

(defun make-vertex (position &optional mapping)
  "General vertex constructor."
  (let ((pos (coerce position 'double-vec)))
    (if mapping
	(make-instance '<mapped-vertex> :position pos :mapping mapping)
	(make-instance '<vertex> :position pos))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Routines
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmethod copy-cell ((vtx <vertex>))
  "Copy constructor for vertex.  The vertex position is freshly allocated."
  (make-instance (class-of vtx) :position (copy-seq (vertex-position vtx))))

(defmethod vertices ((vtx <vertex>)) (list vtx))

(defmethod cell-mapping ((vtx <vertex>))
  "For vertices, this returns a <special-function> evaluating to the vertex
position."
  (with-slots (position) vtx
    (make-instance '<constant-function> :domain-dimension 0
		   :image-dimension (length position) :value position)))

(defmethod l2g ((vtx <vertex>) local-pos)
  "Same as local->global for vertices."
  (assert (equalp local-pos (double-vec)))
  (vertex-position vtx))

(defmethod manifold-dimension ((vtx <vertex>))
  "Anchor for recursive definition."
  (length (vertex-position vtx)))

(defmethod l2Dg ((vtx <vertex>) local-pos)
  (assert (equalp local-pos (double-vec)))
  (make-real-matrix (manifold-dimension vtx) 0))

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
  (when (equalp global-pos (vertex-position vtx))
    (double-vec)))

(defmethod global->local ((vtx <mapped-vertex>) global-pos)
  (assert "Should not be called under the current system.")
  (when (equalp global-pos (evaluate (mapping vtx) (double-vec)))
    (double-vec)))

(defmethod coordinates-inside? ((vtx <vertex>) local-pos)
  (equalp local-pos (double-vec)))

(defmethod local-coordinates-of-midpoint ((cell <vertex>))
  (double-vec))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; refinement
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmethod generate-refine-info ((refcell <vertex>))
  "Construction of the refinement information which is to copy the vertex."
  (with-cell-information (refine-info)
    refcell
    (setq refine-info
	  (vector
	   (make-<child-info>
	    :reference-cell refcell
	    :barycentric-corners
	    (list (list (make-double-vec 1 1.0))) ; not perfect, because no real factor
	    :boundary-paths ())))))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Initialization of the vertex class
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defvar *reference-vertex*
  (let ((refcell (make-vertex (double-vec))))
    (initialize-cell-class refcell ())
    refcell)
  "The reference vertex.")


;;;; Testing: (test-vertex)
(defun test-vertex ()
  (assert (= 0 (manifold-dimension *reference-vertex*)))
  (assert (reference-cell-p *reference-vertex*))
  (reference-cell-p *reference-vertex*)
  (reference-cell *reference-vertex*)
  (dimension (skeleton *reference-vertex*))
  (describe (etable (skeleton *reference-vertex*) 0))
  (describe *reference-vertex*)
  (skeleton-boundary (skeleton *reference-vertex*))
  ;;
  (refine-info *reference-vertex*)
  (cell-class-information *reference-vertex*)
  *
  (describe (refcell-skeleton *reference-vertex*))
  (describe (refcell-refinement-skeleton *reference-vertex* 1))
  (describe (refcell-skeleton *reference-vertex*))
  (describe (refine-globally (skeleton *reference-vertex*)))
  (refcell-refinement-vertices *reference-vertex* 2)
  )

(fl.tests:adjoin-test 'test-vertex)

