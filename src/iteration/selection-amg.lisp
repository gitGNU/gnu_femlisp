;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; selection-amg.lisp
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

;;;; This file contains coarsening and interpolation routines for an
;;;; algebraic multigrid of selection type.  Such AMG algorithms were first
;;;; proposed in 1983 by Ruge and Stueben, for reference on an often-used
;;;; variant see their paper on AMG which appeared in 'Multigrid methods'
;;;; edited by S. F. McCormick in 1987.

(in-package multigrid)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; <selection-amg> class
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass <selection-amg> (<algebraic-mg>)
  ()
  (:documentation "This variant of algebraic multigrid coarsens in a
special way by selecting coarse-grid nodes from the fine-grid nodes.  This
selection is kept in a table, which is then used by the method build-ip to
build the actual interpolation matrix."))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; prolongation
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmethod tentative-prolongation ((amg <selection-amg>) &rest parameters)
  "Simply injection, i.e. P_{CC}=Id_{CC}."
  (let ((mat (getf parameters :matrix))
	(coarse (getf parameters :coarse-nodes)))
    (let ((prol (make-sparse-matrix
		 :print-row-key (print-row-key mat)
		 :print-col-key (print-col-key mat)
		 :row-key->size (row-key->size mat)
		 :col-key->size (col-key->size mat)
		 :keys->pattern (keys->pattern mat))))
      ;; set prolongation entries
      (dolist (node coarse)
	(setf (mref prol node node)
	      (eye (funcall (row-key->size mat) node))))
      ;; augment parameters list with result
      (list* :prolongation prol parameters))))

(defmethod improved-prolongation ((amg <selection-amg>) &rest parameters
				  &key matrix filtered-keys filtered-matrix
				  prolongation fine-nodes
				  &allow-other-keys)
  "Define P_{FC} in a Stueben like manner.  For the fine grid point i we
have to solve for e_i in

$$ a_{ii} e_i = - \sum_{j \in N_i} a_{ij} e_j $$

with N_i being the neighborhood of $i$.  Now $\sum_{j \in N_i} a_{ij} e_j$ is replaced by
$\alpha \sum_{j \in P_i} a_{ij} e_j$ where $\alpha$ is chosen as

$$ \alpha = \frac{\sum_{j \in N_i} a_{ij}}}{\sum_{j \in P_i} a_{ij}} $$

to make the prolongation exact for constant functions."

  (dolist (i fine-nodes)
    (let ((sum-neighboring 0.0)
	  (sum-prolongating 0.0))
      (for-each-key-and-entry-in-row
       #'(lambda (j entry)
	   (when (and (not (eql j i)) (gethash j filtered-keys))
	     (incf sum-neighboring (mref entry 0 0))))
       matrix i)
      (for-each-key-and-entry-in-row
       #'(lambda (j entry)
	   (when (matrix-column prolongation j)
	     (incf sum-prolongating (mref entry 0 0))))
       matrix i)
      (let ((scaled-diagonal-inverse (scal! (- (/ sum-neighboring sum-prolongating))
					    (m/ (mref matrix i i)))))
	(for-each-key-and-entry-in-row
	 #'(lambda (j entry)
	     (when (matrix-column prolongation j)
	       (setf (mref prolongation i j)
		     (m* scaled-diagonal-inverse entry))))
	 filtered-matrix i))))
    
  ;; pass on parameters (including modified prolongation)
  parameters)


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Variants
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; The most interesting variant due to Stueben can be found in
;;; stueben.lisp.

(defclass <custom-coarsening-selection-amg> (<selection-amg>)
  ((select :accessor select :initarg :select))
  (:documentation "If you should already have a suitable selection of
coarse grid points, you may use it in this type of selection-amg.  Your
routine should accept a blackboard of arguments and return two lists: the
coarse nodes and the fine nodes."))

(defmethod choose-coarse-grid ((amg <custom-coarsening-selection-amg>) &rest parameters)
  (multiple-value-bind (coarse fine)
      (apply (select amg) parameters)
    (list* :coarse-nodes coarse :fine-nodes fine parameters)))



