;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; blockit.lisp
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

(in-package :iterations)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; block iterations
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass <block-iteration> (<linear-iteration>)
  ((inner-iteration
    :reader inner-iteration :initform *lu-iteration*
    :initarg :inner-iteration
    :documentation "Iteration which is used to solve for each block.")
   (ordering :accessor ordering :initform nil :initarg :ordering)))

(defgeneric setup-blocks (blockit matrix)
  (:documentation "Setup routine for determining the blocking of unknowns.
Returns a list of blocks where each block is a vector of keys.  May return
a second value which is a list of pair.  Each pair is of the form
start-index/end-index and can be used to filter out different fe
components."))

(defmethod setup-blocks ((bgs <block-iteration>) (smat <sparse-matrix>))
  "If a setup function is provided, it is called.  The default is to use
the standard blocking introduced by the block sparse matrix."
  (values (mapcar #'vector (row-keys smat)) nil))

#+(or)
(defmethod setup-blocks :around ((blockit <block-iteration>) (smat <sparse-matrix>))
  (multiple-value-bind (blocks ranges)
      (call-next-method)
    (if (functionp (ordering blockit))
	(funcall (ordering blockit) blocks)
	blocks)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; block Gauss-Seidel
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass <block-sor> (<block-iteration>)
  ((omega :initform 1.0d0 :initarg :omega)
   (store-p :reader store-p :initform t :initarg :store-p)))

(defclass <block-gauss-seidel> (<block-sor>))

(defmethod make-iterator ((bgs <block-sor>) (smat <sparse-matrix>))
  (multiple-value-bind (blocks ranges)
      (setup-blocks bgs smat)
    (dbg :iter "blocks=~A~%ranges=~A" blocks ranges)
    (let ((diagonal-inverses
	   (loop initially (dbg :iter "computing block-inverses for bgs")
		 for keys in blocks
		 and ranges-tail = ranges then (cdr ranges-tail) ; may be NIL
		 collecting
		 (and (store-p bgs)
		      (m/ (sparse-matrix->matlisp
			   smat :keys keys :ranges (car ranges-tail)))))))
      (make-instance
       '<iterator>
       :matrix smat
       :residual-before nil
       :initialize nil
       :iterate
       #'(lambda (x b r)
	   (dbg :iter "iterate bgs")
	   (loop for block in blocks
		 and ranges-tail = ranges then (cdr ranges-tail) ; may be NIL
		 for block-ranges = (car ranges-tail)
		 and diag-inverse in diagonal-inverses do
		 ;; compute residual on the block (disregarding the ranges)
		 (loop for key across block do
		       (copy! (vec-ref b key) (vec-ref r key))
		       (for-each-key-and-entry-in-row
			#'(lambda (col-key mblock)
			    #+(or)(gemm! -1.0d0 mblock (vec-ref x col-key) 1.0d0 (vec-ref r key))
			    #-(or)(x-=Ay (vec-ref r key) mblock (vec-ref x col-key)))
			smat key))
		 ;; now invert local system
		 (let ((local-r (sparse-vector->matlisp r :keys block :ranges block-ranges))
		       (local-x (sparse-vector->matlisp x :keys block :ranges block-ranges)))
		   (unless diag-inverse
		     (setq diag-inverse 
			   (m/ (sparse-matrix->matlisp
				smat :keys block :ranges (car ranges-tail)))))
		   (gemm! (slot-value bgs 'omega) diag-inverse local-r 1.0d0 local-x)
		   (set-svec-to-local-block x local-x :keys block :ranges block-ranges)))
	   x)
       :residual-after nil))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; custom blocking
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass <setup-blocks-mixin> ()
  ((block-setup :initform nil :initarg :block-setup))
  (:documentation "Executes the given function for determining the block
decomposition."))

(defmethod setup-blocks ((blockit <setup-blocks-mixin>) (smat <sparse-matrix>))
  "Use the setup function."
  (funcall (slot-value blockit 'block-setup) smat))

(defclass <custom-block-gauss-seidel> (<setup-blocks-mixin> <block-gauss-seidel>)
  ())

