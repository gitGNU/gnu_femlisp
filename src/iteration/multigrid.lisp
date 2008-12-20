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

(in-package :fl.multigrid)

;;; This file provides an abstract interface for a multigrid iteration.
;;; The standard methods work for multilevel hierarchies obtained by
;;; algebraically coarsening or uniform geometrical refinement.  It
;;; incorporates the correction scheme, FAS, as well as full multigrid
;;; (F-cycle type iterations).

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; <mg-iteration>
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defgeneric select-smoother (mg-it matrix)
  (:documentation "Select a suitable smoother depending on multigrid
iteration and matrix.")
  (:method (mg-it matrix)
    "Returns the Gauss-Seidel method as default smoother."
    (declare (ignore mg-it matrix))
    (make-instance '<gauss-seidel>)))

(defclass <mg-iteration> (<linear-iteration>)
  ((pre-smoother :reader pre-smoother :initarg :pre-smoother :initarg :smoother)
   (pre-steps :reader pre-steps :initform 1 :initarg :pre-steps)
   (post-smoother :reader post-smoother	:initarg :post-smoother :initarg :smoother)
   (post-steps :reader post-steps :initform 1 :initarg :post-steps)
   (gamma :reader gamma :initform 1 :initarg :gamma)
   (base-level :reader base-level :initform 0 :initarg :base-level)
   (coarse-grid-iteration :reader coarse-grid-iteration
			  :initform (make-instance '<lu>)
			  :initarg :coarse-grid-iteration
			  :initarg :coarse-grid-solver)
   (fmg :reader fmg :initform nil :initarg :fmg)
   (combination :reader combination-type :initform :multiplicative :initarg :combination
		:documentation "Switch between additive and multiplicative
combination of corrections from different levels.  The additive version
should be used as a preconditioner."))
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
;;; coarse-grid-it, pre-smoother-vec, post-smoother-vec, sol-vec, rhs-vec,
;;; res-vec, nr-levels, base-level, top-level, current-level, residual-p.

;;; The following macros are anaphoric macros making access shorter.  They
;;; assume that a variable mg-data is available.  They are not exported.

(defmacro A_ (level) `(aref (getbb mg-data :a-vec) ,level))
(defmacro I_ (level) `(aref (getbb mg-data :i-vec) ,level))
(defmacro R_ (level) `(aref (getbb mg-data :r-vec) ,level))
(defmacro sol_ (level) `(aref (getbb mg-data :sol-vec) ,level))
(defmacro rhs_ (level) `(aref (getbb mg-data :rhs-vec) ,level))
(defmacro res_ (level) `(aref (getbb mg-data :res-vec) ,level))
(defmacro pre-smoother_ (level) `(aref (getbb mg-data :pre-smoother-vec) ,level))
(defmacro post-smoother_ (level) `(aref (getbb mg-data :post-smoother-vec) ,level))
(defmacro residual-p_ (level) `(aref (getbb mg-data :residual-p-vec) ,level))


;;; The FAS scheme needs an additional restriction or projection of the
;;; solution vector.
(defmacro FAS-R_ (level) `(aref (getbb mg-data :fas-r-vec) ,level))

(defmacro with-current-level-data (symbols mg-data &body body)
  "Anaphoric, because it works with the MG-DATA symbol."
  `(let ((mg-data ,mg-data))
    (symbol-macrolet ((current-level (getbb mg-data :current-level))
		      ,@(loop for sym in symbols
			      for l-name = (symbol-name sym)
			      for name = (subseq l-name 0 (1- (length l-name)))
			      collect
			      `(,(intern l-name :fl.multigrid)
				(,(intern name :fl.multigrid) current-level))))
	,@body)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Generic function interface for customization
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defgeneric multilevel-decomposition (mg-it A)
  (:documentation "The central generic function constructing the multilevel
hierarchy."))

(defgeneric ensure-mg-residual (mg-it mg-data)
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

(defmethod ensure-mg-residual ((mg-it <mg-iteration>) mg-data)
  (with-current-level-data (residual-p_l A_l sol_l rhs_l res_l) mg-data
    (unless residual-p_l
      (compute-residual A_l sol_l rhs_l res_l)
      (setq residual-p_l t))))

(defmethod smooth ((mg-it <mg-iteration>) mg-data which)
  (with-current-level-data (pre-smoother_l post-smoother_l sol_l rhs_l res_l residual-p_l)
      mg-data
    (let ((smoother (ecase which
		      (:pre pre-smoother_l)
		      (:post post-smoother_l))))
      (when smoother
	(when (slot-value smoother 'fl.iteration::residual-before)
	  (ensure-mg-residual mg-it mg-data))
	(funcall (slot-value smoother 'fl.iteration::iterate) sol_l rhs_l res_l)
	(setq residual-p_l (slot-value smoother 'fl.iteration::residual-after))))))

(defmethod ensure-sol-rhs-res ((mg-it <mg-iteration>) mg-data level)
  (unless (sol_ level)
    (setf (sol_ level) (make-domain-vector-for
			(A_ level) (multiplicity (sol_ (1+ level))))))
  (unless (rhs_ level)
    (setf (rhs_ level) (make-image-vector-for
			(A_ level) (multiplicity (rhs_ (1+ level))))))
  (unless (res_ level)
    (setf (res_ level) (make-image-vector-for
			(A_ level) (multiplicity (res_ (1+ level)))))))

(defmethod restrict ((mg-it <mg-iteration>) mg-data)
  "The basic method for restriction restricts the residual, decrements the
level and clears the residual-p flag."
  (with-current-level-data (res_l) mg-data
    (let* ((l-1 (1- current-level)))
      (ensure-mg-residual mg-it mg-data)
      (ensure-sol-rhs-res mg-it mg-data l-1)
      (if (getbb mg-data :r-vec)
	  (gemm! 1.0 (R_ l-1) res_l 0.0 (res_ l-1))
	  (gemm! 1.0 (I_ l-1) res_l 0.0 (res_ l-1) :tn))
      (dbg :iter "~A" (norm (res_ l-1))))))

(defmethod restrict :after ((mg-it <mg-iteration>) mg-data)
  (with-current-level-data (residual-p_l) mg-data
    (decf current-level)
    (setq residual-p_l nil)))

(defmethod restrict :after ((mg-it <correction-scheme>) mg-data)
  (with-current-level-data (sol_l res_l rhs_l residual-p_l)
      mg-data
    (x<-0 sol_l)
    (copy! res_l rhs_l)
    (setq residual-p_l t)))

(defmethod restrict :after ((mg-it <fas>) mg-data)
  (with-current-level-data (sol_l res_l rhs_l A_l residual-p_l FAS-R_l)
      mg-data
    (copy! res_l rhs_l)
    (gemm! 1.0 FAS-R_l (sol_ (1+ current-level)) 0.0 sol_l)
    (gemm! 1.0 A_l sol_l 1.0 rhs_l)
    (setf residual-p_l t)))

;;; Note that no primary method for prolongate is defined for
;;; <mg-iteration>.  Thus, this class has to be merged with
;;; <correction-scheme> or <fas>.
(defmethod prolongate ((mg-it <correction-scheme>) mg-data)
  (with-current-level-data (sol_l I_l) mg-data
    (gemm! 1.0 I_l sol_l 1.0 (sol_ (1+ current-level)) :nn)))

(defmethod prolongate ((mg-it <fas>) mg-data)
  "This version of FAS prolongation uses the res_ field on the coarser
level for computing the correction to be prolongated."
  (with-current-level-data (sol_l res_l I_l FAS-R_l) mg-data
    (copy! sol_l res_l)
    (gemm! -1.0 FAS-R_l (sol_ (1+ current-level)) 1.0 res_l)
    (gemm! 1.0 I_l res_l 1.0 (sol_ (1+ current-level)) :nn)))

(defmethod prolongate :after ((mg-it <mg-iteration>) mg-data)
  (with-current-level-data (residual-p_l) mg-data
    (incf current-level)
    (setf residual-p_l nil)))

(defmethod lmgc ((mg-it <mg-iteration>) mg-data)
  (with-items (&key base-level coarse-grid-it) mg-data
    (with-current-level-data (sol_l res_l rhs_l residual-p_l)
	mg-data
      (cond
	((= current-level base-level)
	 (funcall (slot-value coarse-grid-it 'iterate) sol_l rhs_l res_l)
	 (setq residual-p_l (slot-value coarse-grid-it 'fl.iteration::residual-after)))
	(t
	 (ecase (combination-type mg-it)
	   (:multiplicative
	    (smooth mg-it mg-data :pre)
	    (restrict mg-it mg-data))
	   (:additive
	    (restrict mg-it mg-data)
	    (incf current-level)
	    (smooth mg-it mg-data :pre)
	    (decf current-level)))
	 (loop for i from (cond ((zerop (slot-value mg-it 'gamma)) 0)
				((= current-level base-level) 1)
				(t (slot-value mg-it 'gamma)))
	       downto 1 do
	       (lmgc mg-it mg-data)
	       (unless (= i 1)
		 (ensure-mg-residual mg-it mg-data)))
	 (ecase (combination-type mg-it)
	   (:additive
	    (incf current-level)
	    (smooth mg-it mg-data :post)
	    (decf current-level)
	    (prolongate mg-it mg-data))
	   (:multiplicative
	    (prolongate mg-it mg-data)
	    (smooth mg-it mg-data :post)))
	 )))))

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
  (declare (optimize debug))
  (dbg :iter "making iterator for <mg-iteration>")
  (with-slots (gamma pre-smoother pre-steps post-smoother post-steps
		     coarse-grid-iteration)
    mg-it
    (let ((mg-data (multilevel-decomposition mg-it A)))
      (with-items (&key a-vec nr-levels coarse-grid-it pre-smoother-vec post-smoother-vec
			current-level base-level top-level sol-vec rhs-vec res-vec
			residual-p-vec)
	  mg-data
	
	;; setup mg-data
	(setq nr-levels (length a-vec))
	(setq base-level (base-level mg-it))
	(setq top-level (1- nr-levels))
	(assert (<= base-level top-level))
	(setq sol-vec (make-array nr-levels :initial-element nil)
	      rhs-vec (make-array nr-levels :initial-element nil)
	      res-vec (make-array nr-levels :initial-element nil)
	      residual-p-vec (make-array nr-levels :initial-element nil))
	(setq pre-smoother-vec (make-array nr-levels)
	      post-smoother-vec (make-array nr-levels))
	;; smoothers on levels above base-level
	(loop for level from top-level above base-level do
	     (let* ((A-l (aref a-vec level))
		    (default-smoother
		     (when (or (not (slot-boundp mg-it 'pre-smoother))
			       (not (slot-boundp mg-it 'post-smoother)))
		       (select-smoother mg-it A-l)))
		    (pre-smoother-l
		     (if (slot-boundp mg-it 'pre-smoother)
			 (if (functionp pre-smoother)
			     (funcall pre-smoother level)
			     pre-smoother)
			 default-smoother))
		    (post-smoother-l
		     (if (slot-boundp mg-it 'post-smoother)
			 (if (functionp post-smoother)
			     (funcall post-smoother level)
			     post-smoother)
			 default-smoother))
		    (pre-smoother-it
		     (and (plusp pre-steps)
			  (make-iterator pre-smoother-l A-l)))
		    (post-smoother-it
		     (and (plusp post-steps)
			  (or (and (eq pre-smoother-l post-smoother-l) pre-smoother-it)
			      (make-iterator post-smoother-l A-l)))))
	       (setf (aref pre-smoother-vec level)
		     (and pre-smoother-it (product-iterator pre-smoother-it pre-steps)))
	       (setf (aref post-smoother-vec level)
		     (and post-smoother-it (product-iterator post-smoother-it post-steps)))))
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
	     (setf (residual-p_ top-level) t)
	     (setf current-level top-level)
	     
	     (if (fmg mg-it)
		 (f-cycle mg-it mg-data)
		 (lmgc mg-it mg-data))
	     x)
	 :residual-after
	 (if (= base-level top-level)
	     (slot-value coarse-grid-it 'fl.iteration::residual-after)
	     (aand (aref post-smoother-vec top-level)
		   (slot-value it 'fl.iteration::residual-after))))
	))))
