;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; multigrid.lisp
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

(in-package :multigrid)

;;; This file provides an abstract interface for a multigrid iteration.
;;; The standard methods work for multilevel hierarchies obtained by
;;; algebraically coarsening or uniform geometrical refinement.  It
;;; incorporates the correction scheme, FAS, as well as full multigrid
;;; (F-cycle type iterations).

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; <mg-iteration>
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defparameter *default-smoother* *gauss-seidel*
  "The default smoothing iteration for the multigrid scheme.  This can be
Gauss-Seidel, but for many applications something special will give better
results.")

(defparameter *default-coarse-grid-iteration* *lu-iteration*
  "The default coarse grid solver is a sparse LU decompositions.")

(defclass <mg-iteration> (<linear-iteration>)
  ((pre-smooth :reader pre-smooth :initform *default-smoother*
	       :initarg :pre-smooth)
   (pre-steps :reader pre-steps :initform 1 :initarg :pre-steps)
   (post-smooth :reader post-smooth :initform *default-smoother*
		:initarg :post-smooth)
   (post-steps :reader post-steps :initform 1 :initarg :post-steps)
   (gamma :reader gamma :initform 1 :initarg :gamma)
   (base-level :reader base-level :initform 0 :initarg :base-level)
   (coarse-grid-iteration :reader coarse-grid-iteration
			  :initform *default-coarse-grid-iteration*
			  :initarg :coarse-grid-iteration
			  :initarg :coarse-grid-solver)
   (fmg :reader fmg :initform nil :initarg :fmg))
  (:documentation "The multigrid iteration is a linear iteration specially
suited for the solution of systems of equations with elliptic terms.  In
ideal situations, it solves such systems with optimal complexity.  It is a
complicated linear iteration, which consists of applying simple linear
iterators as smoothers on a hierarchy of grids.  This grid hierarchy is
obtained either by discretizing on successively refined meshes (geometric
multigrid) or it is constructed from matrix information alone (algebraic
multigrid).

The <mg-iteration> is not intended to be used directly.  Only incorporating
mixins like <correction-scheme> or <fas> results in concrete classes like
<algebraic-mg>."))

(defclass <correction-scheme> ()
  ()
  (:documentation "This is a mixin-class which yields the correction
scheme.  It should be merged !before! <mg-iteration> for standard CLOS
class precedence."))

(defclass <fas> ()
  ()
  (:documentation "This is a mixin-class for <mg-iteration> which yields
the behaviour of Brandt's FAS scheme.  It should be merged !before!
<mg-iteration> or the derived class <geometric-mg> when using standard CLOS
class precedence."))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Multigrid data
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; This is a blackboard with the following items: a-vec, i-vec, r-vec,
;;; coarse-grid-it, pre-smooth-vec, post-smooth-vec, sol-vec, rhs-vec,
;;; res-vec, nr-levels, base-level, top-level, current-level, residual-p.

;;; The following macros are anaphoric macros making access shorter.  They
;;; assume that a variable mg-data is available.  They are not exported.

(defmacro A_ (level) `(aref (getbb mg-data :a-vec) ,level))
(defmacro I_ (level) `(aref (getbb mg-data :i-vec) ,level))
(defmacro R_ (level) `(aref (getbb mg-data :r-vec) ,level))
(defmacro sol_ (level) `(aref (getbb mg-data :sol-vec) ,level))
(defmacro rhs_ (level) `(aref (getbb mg-data :rhs-vec) ,level))
(defmacro res_ (level) `(aref (getbb mg-data :res-vec) ,level))
(defmacro pre-smoother_ (level) `(aref (getbb mg-data :pre-smooth-vec) ,level))
(defmacro post-smoother_ (level) `(aref (getbb mg-data :post-smooth-vec) ,level))

(define-symbol-macro current-level (getbb mg-data :current-level))
(define-symbol-macro A_l (A_ current-level))
(define-symbol-macro I_l (I_ current-level))
(define-symbol-macro R_l (R_ current-level))
(define-symbol-macro sol_l (sol_ current-level))
(define-symbol-macro rhs_l (rhs_ current-level))
(define-symbol-macro res_l (res_ current-level))
(define-symbol-macro pre-smoother_l (pre-smoother_ current-level))
(define-symbol-macro post-smoother_l (post-smoother_ current-level))

(define-symbol-macro residual-p (getbb mg-data :residual-p))

;;; The FAS scheme needs an additional restriction or projection of the
;;; solution vector.
(defmacro FAS-R_ (level) `(aref (getbb mg-data :fas-r-vec) ,level))
(define-symbol-macro FAS-R_l (FAS-R_ current-level))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Generic function interface for customization
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defgeneric multilevel-decomposition (mg-it A)
  (:documentation "The central generic function constructing the multilevel
hierarchy."))

(defgeneric ensure-residual (mg-it mg-data)
  (:documentation "Ensure residual on current-level."))

(defgeneric ensure-sol-rhs-res (mg-it mg-data level)
  (:documentation "Ensure allocation of sol/rhs/res on level."))

(defgeneric prolongate (mg-it mg-data)
  (:documentation "Prolongation from current level.  The current level is
incremented after the operation."))

(defgeneric restrict (mg-it mg-data)
  (:documentation "Restriction from current level.  The current level is
decremented after the operation."))

(defgeneric smooth (mg-it mg-data which)
  (:documentation "Smooth on current level.  The argument which should be
one of the keywords :pre or :post and determines which smoothing is
to be performed."))

(defgeneric lmgc (mg-it mg-data)
  (:documentation "Performs an lmgc cycle on current level."))

(defgeneric f-cycle (mg-it mg-data)
  (:documentation "Performs an F-cycle on current level."))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Default methods for <mg-iteration> and specializations for
;;;; <correction-scheme>
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmethod ensure-residual ((mg-it <mg-iteration>) mg-data)
  (unless residual-p
    (compute-residual A_l sol_l rhs_l res_l)
    (setq residual-p t)))

(defmethod smooth ((mg-it <mg-iteration>) mg-data which)
  (let ((smoother (ecase which
		    (:pre pre-smoother_l)
		    (:post post-smoother_l))))
    (when smoother
      (when (slot-value smoother 'iteration::residual-before)
	(ensure-residual mg-it mg-data))
      (funcall (slot-value smoother 'iteration::iterate) sol_l rhs_l res_l)
      (setq residual-p (slot-value smoother 'iteration::residual-before)))))

(defmethod ensure-sol-rhs-res ((mg-it <mg-iteration>) mg-data level)
  (unless (sol_ level)
    (setf (sol_ level) (make-row-vector-for
			(A_ level) (multiplicity (sol_ (1+ level))))))
  (unless (rhs_ level)
    (setf (rhs_ level) (make-column-vector-for
			(A_ level) (multiplicity (rhs_ (1+ level))))))
  (unless (res_ level)
    (setf (res_ level) (make-column-vector-for
			(A_ level) (multiplicity (res_ (1+ level)))))))

(defmethod restrict ((mg-it <mg-iteration>) mg-data)
  "The basic method for restriction restricts the residual, decrements the
level and clears the residual-p flag."
  (let ((l-1 (1- current-level)))
    (ensure-residual mg-it mg-data)
    (ensure-sol-rhs-res mg-it mg-data l-1)
    (if (getbb mg-data :r-vec)
	(gemm! 1.0d0 (R_ l-1) res_l 0.0d0 (res_ l-1))
	(gemm! 1.0d0 (I_ l-1) res_l 0.0d0 (res_ l-1) :tn))
;;    (whereas ((global-i-mat (getbb mg-data :global-i-vec)))
      
	))

(defmethod restrict :after ((mg-it <mg-iteration>) mg-data)
  (decf current-level)
  (setq residual-p nil))

(defmethod restrict :after ((mg-it <correction-scheme>) mg-data)
  (x<-0 sol_l)
  (x<-y rhs_l res_l)
  (setq residual-p t))

(defmethod restrict :after ((mg-it <fas>) mg-data)
  (x<-y rhs_l res_l)
  (gemm! 1.0d0 FAS-R_l (sol_ (1+ current-level)) 0.0d0 sol_l)
  (gemm! 1.0d0 A_l sol_l 1.0d0 rhs_l)
  (setf residual-p t))

;;; Note that no primary method for prolongate is defined for
;;; <mg-iteration>.  Thus, this class has to be merged with
;;; <correction-scheme> or <fas>.
(defmethod prolongate ((mg-it <correction-scheme>) mg-data)
  (let ((l+1 (1+ current-level)))
    (gemm! 1.0d0  I_l sol_l 1.0d0 (sol_ l+1) :nn)))

(defmethod prolongate ((mg-it <fas>) mg-data)
  "This version of FAS prolongation uses the res_ field on the coarser
level for computing the correction to be prolongated."
  (let ((l+1 (1+ current-level)))
    (x<-y res_l sol_l)
    (gemm! -1.0d0 FAS-R_l (sol_ l+1) 1.0d0 res_l)
    (gemm! 1.0d0 I_l res_l 1.0d0 (sol_ l+1)) :nn))

(defmethod prolongate :after ((mg-it <mg-iteration>) mg-data)
  (incf current-level)
  (setq residual-p nil))

(defmethod lmgc ((mg-it <mg-iteration>) mg-data)
  (with-items (&key current-level base-level coarse-grid-it) mg-data
    (cond
      ((= current-level base-level)
       (funcall (slot-value coarse-grid-it 'iterate) sol_l rhs_l res_l)
       (setq residual-p (slot-value coarse-grid-it 'iteration::residual-after)))
      (t
       (smooth mg-it mg-data :pre)
       (restrict mg-it mg-data)
       (loop for i from (cond ((zerop (slot-value mg-it 'gamma)) 0)
			      ((= current-level base-level) 1)
			      (t (slot-value mg-it 'gamma)))
	     downto 1 do
	     (lmgc mg-it mg-data)
	     (unless (= i 1)
	       (ensure-residual mg-it mg-data)))
       (prolongate mg-it mg-data)
       (smooth mg-it mg-data :post)
       ))))

(defmethod f-cycle ((mg-it <mg-iteration>) mg-data)
  (with-items (&key base-level top-level) mg-data
    (loop for level from top-level above base-level do
	  (restrict mg-it mg-data))
    (loop for level from base-level upto top-level do
	  (lmgc mg-it mg-data)
	  (unless (= level top-level)
	    (prolongate mg-it mg-data)))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; make-iterator
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmethod make-iterator ((mg-it <mg-iteration>) (A <sparse-matrix>))
  "This is the method for a multigrid iteration on a uniformly refined
grid, or an algebraic multigrid iteration."
  (dbg :iter "making iterator for <mg-iteration>")
  (with-slots (gamma pre-smooth pre-steps post-smooth post-steps
		     coarse-grid-iteration)
    mg-it
    (let ((mg-data (multilevel-decomposition mg-it A)))
      (with-items (&key a-vec nr-levels coarse-grid-it pre-smooth-vec post-smooth-vec
			current-level base-level top-level sol-vec rhs-vec res-vec)
	  mg-data
	
	;; setup mg-data
	(setq nr-levels (length a-vec))
	(setq base-level (base-level mg-it))
	(setq top-level (1- nr-levels))
	(assert (<= base-level top-level))
	(setq sol-vec (make-array nr-levels :initial-element nil)
	      rhs-vec (make-array nr-levels :initial-element nil)
	      res-vec (make-array nr-levels :initial-element nil))
	(setq pre-smooth-vec (make-array nr-levels)
	      post-smooth-vec (make-array nr-levels))
	;; smoothers on levels above base-level
	(loop for level from top-level above base-level
	      for pre-smooth-l = (if (functionp pre-smooth)
				     (funcall pre-smooth level)
				     pre-smooth)
	      for post-smooth-l = (if (functionp post-smooth)
				      (funcall post-smooth level)
				      post-smooth)
	      for pre-smooth-it = (and (plusp pre-steps)
				       (make-iterator pre-smooth-l (aref a-vec level)))
	      for post-smooth-it = (and (plusp post-steps)
					(or (and (eq pre-smooth-l post-smooth-l) pre-smooth-it)
					    (make-iterator post-smooth-it (aref a-vec level))))
	      do
	      (setf (aref pre-smooth-vec level)
		    (and pre-smooth-it (product-iterator pre-smooth-it pre-steps)))
	      (setf (aref post-smooth-vec level)
		    (and post-smooth-it (product-iterator post-smooth-it post-steps))))
	;; coarse-grid iterator
	(when (or (= top-level base-level) (> gamma 0))
	  (setq coarse-grid-it
		(make-iterator
		 (etypecase coarse-grid-iteration
		   (<linear-iteration> coarse-grid-iteration)
		   (<solver> (make-instance '<solver-iteration> :solver coarse-grid-iteration)))
		  (aref a-vec base-level))))
	
	;; and return the multigrid iterator
	(make-instance
	 '<iterator>
	 :matrix A
	 :residual-before t
	 :initialize nil
	 :iterate
	 #'(lambda (x b r)
	     ;; set top-level solution, rhs, res to given values
	     (setf (sol_ top-level) x
		   (rhs_ top-level) b
		   (res_ top-level) r)
	     ;; initialize solver state
	     (setq residual-p t)
	     (setf current-level top-level)
	     
	     (if (fmg mg-it)
		 (f-cycle mg-it mg-data)
		 (lmgc mg-it mg-data))
	     x)
	 :residual-after
	 (if (= base-level top-level)
	     (slot-value coarse-grid-it 'iteration::residual-after)
	     (aand (aref post-smooth-vec top-level)
		   (slot-value it 'iteration::residual-after))))
	))))
