;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; amg.lisp
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

(in-package :fl.multigrid)

;;; This file provides the algebraic multigrid iteration.

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; <algebraic-mg>
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defvar *amg-cg-max-size* 1000
  "Maximum size of the coarse grid in the AMG algorithm.")

(defclass <algebraic-mg> (<correction-scheme> <mg-iteration>)
  ((max-depth :reader max-depth :initform most-positive-fixnum
	      :initarg :max-depth)
   (cg-max-size :reader cg-max-size :initform *amg-cg-max-size*
		:initarg :cg-max-size)
   (output :initform nil :initarg :output))
  (:documentation "The algebraic multigrid iteration is a multigrid
iteration where the hierarchy of problems is derived from the fine-grid
matrix.  Usually, an algebraic multigrid will use the same iterator as its
geometric counterpart."))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; multilevel-decomposition
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmethod multilevel-decomposition ((amg <algebraic-mg>) mat)
  "This is the standard method for building up an AMG hierarchy.  Thus, it
should be usable by a large variety of AMG algorithms.  Its purpose is to
repeat coarsening steps until either max-depth is larger than the given
value or the size of the problem is smaller than min-size.  The generic
function 'coarsen' is called to do a single coarsening step."
  (with-slots (output max-depth cg-max-size) amg
    (let ((mats (list mat))
	  (pmats (list nil))
	  (rmats (list nil)))
      (when output
	(format t "~&~A-Coarsening~%" (class-name (class-of amg)))
	(format t "Level   n-rows  n-entries~%")
	(format t "~3D~6T~8D~16T~9D~%" 0 (nrows mat) (nr-of-entries mat)))
      (loop for i from 1 below max-depth
	    until (<= (nr-nonempty-rows mat) cg-max-size) do
	    (let ((coarsening (coarsen amg mat)))
	      (unless coarsening (return nil))
	      (destructuring-bind (&key coarse-grid-matrix prolongation restriction)
		  coarsening
		(push coarse-grid-matrix mats)
		(push prolongation pmats)
		(push restriction rmats)
		(setq mat coarse-grid-matrix)))
	    (when output
	      (format t "~3D~6T~8D~16T~9D~%" i (nrows mat) (nr-of-entries mat))))
      (when output (terpri))
      (blackboard :a-vec (coerce mats 'simple-vector)
		  :i-vec (coerce pmats 'simple-vector)
		  :r-vec (coerce rmats 'simple-vector)))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; coarsen - do one coarsening step
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defgeneric coarsen (amg mat)
  (:documentation "Given AMG and matrix, this generic function returns
coarse-grid matrix, interpolation and restriction matrices for one coarsening
step."))

(defgeneric preprocess-matrix (amg mat)
  (:documentation "Eliminates Dirichlet or slave degrees of freedom which
can be handled well by smoothing."))

(defgeneric prolongation (amg mat)
  (:documentation "Computes a prolongation matrix from amg and mat.  This
is often the essence of an AMG method."))

(defgeneric restriction (amg mat prol)
  (:documentation "Compute a restriction matrix from amg, mat and
prolongation."))

(defgeneric coarse-grid-matrix (amg mat prolongation restriction)
  (:documentation "Computes a coarse-grid matrix from amg, mat,
prolongation and restriction."))
  

(defmethod coarsen ((amg <algebraic-mg>) mat)
  "This is the default method for AMG coarsening.  It performs the following steps:

1. (prolongation) Build up a prolongation P.  Its domain is the coarse-grid.

2. (restriction) Compute a restriction R.  Usually this will be R=I^t.

3. (compute-coarse-grid-matrix) Build the coarse-grid matrix A_C by computing
the Galerkin product A_C = I^t A I."
  (let ((new-mat (preprocess-matrix amg mat)))
    (unless (mzerop new-mat)
      (let ((pmat (prolongation amg new-mat)))
	(unless (mzerop pmat)
	  (let ((rmat (restriction amg new-mat pmat)))
	    (let ((cgm (coarse-grid-matrix amg new-mat pmat rmat)))
	      #+(or)(assert (every #'plusp
			     (matlisp::store
			      (eig (sparse-matrix->matlisp
				    cgm :keys (row-keys cgm))))))
	      (list :coarse-grid-matrix cgm :prolongation pmat :restriction rmat))))))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Preprocessing
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun slave-dof-p (key mat)
  "Checks if nothing depends on key and if key depends on other keys."
  (for-each-key-in-col
   #'(lambda (row-key)
       (or (eq key row-key)
	   (mzerop (mref mat row-key key))
	   (return-from slave-dof-p nil)))
   mat key)
  (for-each-key-in-row
   #'(lambda (col-key)
       (or (eq key col-key)
	   (mzerop (mref mat key col-key))
	   (return-from slave-dof-p t)))
   mat key)
  nil)

(defun dirichlet-dof-p (key mat)
  "Checks if key does not depend on other keys, so that it can be kept on
the fine grid."
  (for-each-key-in-row
   #'(lambda (col-key)
       (or (eq key col-key)
	   (mzerop (mref mat key col-key))
	   (return-from dirichlet-dof-p nil)))
   mat key)
  t)

(defun slave-or-dirichlet-dof-p (key mat)
  "Checks if key is a hanging node or a Dirichlet node."
  (or (slave-dof-p key mat) (dirichlet-dof-p key mat)))

(defmethod preprocess-matrix ((amg <algebraic-mg>) mat)
  "Default method which eliminates slave nodes and Dirichlet nodes for
sparse matrices."
  (let ((active (make-hash-table)))
    (for-each-row-key
     #'(lambda (key)
	 (unless (slave-or-dirichlet-dof-p key mat)
	   (setf (gethash key active) t)))
     mat)
    (extended-extract mat active)
    ))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; prolongation - compute a prolongation
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; Interestingly, the construction of the prolongation follows an abstract
;;; pattern which is identical for AMG of selection or aggregation type.
;;; Yet, to keep maximum flexibility, we implement this pattern as pipe
;;; operations working on blackboards.

(defgeneric filtered-matrix (amg blackboard)
  (:documentation "Precondition: items :matrix.  Postcondition: :matrix and
:filtered-matrix."))

(defgeneric choose-coarse-grid (amg blackboard)
  (:documentation "Pre-condition: an item :matrix has to be supplied on
@arg{blackboard}.  Post-condition: Items :matrix and :coarse-grid in the
result.  In the case of selection-amg this will be a set of selected
indices of the fine-grid matrix graph.  In the case of amg of aggregation
type, this is a set of mutually disjoint sets of fine-grid indices."))

(defgeneric tentative-prolongation (amg blackboard)
  (:documentation "Precondition: the items :matrix and
:coarse-grid-nodes have to be supplied on the blackboard.  Post-condition:
Items :matrix and :tentative-prolongation on the blackboard.  In the case
of selection-amg this will usually be injection, in the case of
aggregation-amg this will be the piecewise constant prolongation."))

(defgeneric improved-prolongation (amg blackboard)
  (:documentation "Precondition: items :mat and :prol in the rest
parameters.  Postcondition: same."))

(defmethod prolongation (amg mat)
  "General definition for construction of a prolongation.  For
incorporating other AMG algorithms, you should first try to keep this
routinge as it is and define modifications for the called generic
functions."
  (let ((bb (blackboard :matrix mat)))
    (filtered-matrix amg bb)
    (choose-coarse-grid amg bb)
    (tentative-prolongation amg bb)
    (improved-prolongation amg bb)
    (getbb bb :prolongation)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; restriction
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmethod restriction ((amg <algebraic-mg>) mat prol)
  "This is the default method for defining the restriction in AMG algorithms as
the transpose of the prolongation."
  (declare (ignore mat))
  (transpose prol))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; compute-coarse-grid-matrix
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun galerkin-product (R A P)
  "Builds the Galerkin product R A P.  This function works for every type
of matrices for which the row- and column-loop macros are defined.  This
procedure should be inlined into an environment where types are known for
avoiding generic arithmetic."
  (let* ((row-key->size (row-key->size R))
	 (col-key->size (col-key->size P))
	 (result (make-sparse-matrix
		  :print-row-key (print-row-key R)
		  :print-col-key (print-col-key P)
		  :row-key->size row-key->size
		  :col-key->size col-key->size
		  :keys->pattern
		  #'(lambda (row-key col-key)
		      (full-crs-pattern (funcall row-key->size row-key)
					(funcall col-key->size col-key))))))
    (for-each-row-key
     #'(lambda (l)
	 (for-each-key-and-entry-in-row
	  #'(lambda (m A_lm)		; loop through all entries A_lm of Amat
	      (for-each-key-and-entry-in-col
	       #'(lambda (k R_kl)	; loop through column l of R
		   (assert (= (ncols R_kl) (nrows A_lm)))
		   (let ((R_kl*A_lm (m* R_kl A_lm)))
		     (for-each-key-and-entry-in-row
		      #'(lambda (n P_mn) ; loop through row m of P
			  (let ((result-entry (mref result k n))) 
			    (assert (= (nrows R_kl*A_lm) (nrows result-entry)))
			    (assert (= (ncols R_kl*A_lm) (nrows P_mn)))
			    (assert (= (ncols P_mn) (ncols result-entry)))
			    (gemm! 1.0 R_kl*A_lm P_mn 1.0 result-entry)))
		      P m)))
	       R l))
	  A l))
     A)
    result))

(defmethod coarse-grid-matrix ((amg <algebraic-mg>) (A <sparse-matrix>)
			       (P <sparse-matrix>) (R <sparse-matrix>))
  "This is the standard method of generating a coarse-grid matrix."
  (galerkin-product R A P))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; AMG variants
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; see selection-amg.lisp
;;; see aggregation-amg.lisp

;;; Testing


