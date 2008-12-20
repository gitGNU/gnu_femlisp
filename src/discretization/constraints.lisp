;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;  constraints.lisp
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

(in-package :fl.discretization)

(file-documentation
 "This file contains routines for handling constraints.  Constraints
originate from geometry (identified boundaries), problem (essential
boundary conditions), or discretization (hanging nodes).")

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Essential constraints
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defgeneric essential-boundary-constraints (problem ansatz-space &key level where interface)
  (:documentation "Computation of essential constraints.  Should probably
be incorporated into the ansatz-space definition."))

(defun compute-essential-boundary-constraints
    (ansatz-space &key level (where :surface) interface)
  "Maybe there will be a general solution somewhen, but for the moment we
delegate it to the specialized assembly."
  (essential-boundary-constraints
   (problem ansatz-space) ansatz-space
   :level level :where where :interface interface))

(defun constrained-interpolation-matrix (ansatz-space &key level where imat (type :local))
  "The multigrid algorithm needs an interpolation which satisfies the
constraints like essential or periodic boundary conditions."
  (let ((imat (interpolation-matrix ansatz-space :level level :imat imat :type type)))
    ;; Could this create problems for FAS and non-Dirichlet b.c.??  Should
    ;; prolongation respect the constraints for level+1 and restriction for
    ;; level?  Alternatively, this could be achieved by enforcing the
    ;; constraints after the operation, as in UG.
    (let ((constraints-P (compute-essential-boundary-constraints
			  ansatz-space :level (1+ level) :where where)))
      (remove-projection-range imat constraints-P :row-p t))
    (let ((constraints-P (compute-essential-boundary-constraints
			  ansatz-space :level level :where where)))
      (remove-projection-range imat constraints-P :column-p t))
    imat))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Hanging-node constraints
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun eliminate-hanging-node-constraints-from-matrix
    (A constraints &key (directions '(:left :right))
     &aux constraints-p)
  "Eliminates the constraints destructively from A.  directions is a list
consisting of the two symbols :left and :right which determine if
elimination is done wrt row or column indices.  Returns t if constraints
have been eliminated."
  (declare (optimize (debug 3)))
  (when (member :left directions)
    (for-each-row-key
     #'(lambda (i)
	 (when (and (matrix-row constraints i)
		    (not (matrix-column constraints i)))
	   (setq constraints-p t)
	   (for-each-key-and-entry-in-row
	    #'(lambda (j Aij)
		(for-each-key-and-entry-in-row
		 #'(lambda (k Cik)
		     (gemm! 1.0 Cik Aij 1.0 (mref A k j) :tn))
		 constraints i))
	    A i)
	   (remove-row A i)))
     A))
  (when (member :right directions)
    (for-each-col-key
     #'(lambda (j)
	 (when (and (matrix-row constraints j)
		    (not (matrix-column constraints j)))
	   (setq constraints-p t)
	   (for-each-key-and-entry-in-col
	    #'(lambda (i Aij)
		(for-each-key-and-entry-in-row
		 #'(lambda (k Cjk)
		     (gemm! 1.0 Aij Cjk 1.0 (mref A i k)))
		 constraints j))
	    A j)
	   (remove-column A j)))
     A))
  ;; return if constraints have been eliminated
  constraints-p)

(defun add-local-part! (A B index-table &key (directions '(:left :right)))
  "Adds the parts for the specified indices from B to A.  This operation is
destructive for A!"
  (loop for i being each hash-key of index-table do
	(when (member :right directions)
	  (when (matrix-row B i)
	    (for-each-key-and-entry-in-row
	     #'(lambda (j Bij)
		 (m+! Bij (mref A i j)))
	     B i)))
	(when (member :left directions)
	  (when (matrix-column B i)
	    (for-each-key-and-entry-in-col
	     #'(lambda (j Bji)
		 (unless (gethash j index-table)
		   (m+! Bji (mref A j i))))
	     B i)))))

(defun hanging-node-constraints (ansatz-space &key level ip-type)
  "Returns inter-level constraints between level and level+1.  ip-type=t
sets the constraint matrix to identity for the level-unknowns."
  (let ((constraints-Q
	 (interpolation-matrix
	  ansatz-space
	  :region (refinement-interface (hierarchical-mesh ansatz-space) :level level)
	  :level level))
	;;(constraints-Q (make-ansatz-space-automorphism ansatz-space))
	(constraints-P (make-ansatz-space-automorphism ansatz-space))
	(constraints-r (make-ansatz-space-vector ansatz-space)))
    (for-each-row-key #'(lambda (key) (row<-id constraints-P key)) constraints-Q)
    (when ip-type
      (for-each-col-key
       #'(lambda (key)
	   (unless (matrix-row constraints-Q key)
	     (row<-id constraints-Q key)))
       constraints-Q))
    (values constraints-P constraints-Q constraints-r)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Old version of handling of constraints (used in fedisc.lisp)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun hanging-nodes-constraints (ansatz-space &key level)
  "These are inter-level constraints."
  (let ((imat (interpolation-matrix
	       ansatz-space
	       :region (refinement-interface (hierarchical-mesh ansatz-space) :level level)
	       :level level))
	;;(constraints-Q (make-ansatz-space-automorphism ansatz-space))
	(constraints-P (make-ansatz-space-automorphism ansatz-space))
	(constraints-r (make-ansatz-space-vector ansatz-space)))
    (for-each-row-key #'(lambda (key) (row<-id constraints-P key)) imat)
    (values constraints-P imat constraints-r)))

(defun combined-constraints (old-P old-Q old-r new-P new-Q new-r)
  "Constructs a combined constraint tuple given old and new constraints."
  (let ((extended-P (combined-projection old-P new-P))
	(combined-Q
	 (m+ new-Q
	     (sparse-m*
	      old-Q (extend-by-identity new-Q (column-table old-Q)
					:ignore (row-table new-P) :copy t)
	      :sparsity :A))))
    (values extended-P combined-Q (m+ old-r new-r))))

(defun assemble-constraints (ansatz-space)
  (macrolet ((set-constraints (name P Q r)
	       `(symbol-macrolet ((line (properties ansatz-space)))
		 (setf (getf line ,(intern (concatenate 'string name "-P") :keyword)) ,P)
		 (setf (getf line ,(intern (concatenate 'string name "-Q") :keyword)) ,Q)
		 (setf (getf line ,(intern (concatenate 'string name "-R") :keyword)) ,R))))
	   
    ;; constraints assembly
    (let ((constraints-P (make-ansatz-space-automorphism ansatz-space))
	  (constraints-Q (make-ansatz-space-automorphism ansatz-space))
	  (constraints-r (make-ansatz-space-vector ansatz-space))
	  (h-mesh (hierarchical-mesh ansatz-space)))
      
      ;; hanging nodes constraints
      (loop for level from (1- (top-level h-mesh)) downto 0
	    do
	    (multiple-value-bind (hanging-P hanging-Q hanging-r)
		(hanging-nodes-constraints ansatz-space :level level)
	      (multiple-value-setq (constraints-P constraints-Q constraints-r)
		(combined-constraints constraints-P constraints-Q constraints-r
				      hanging-P hanging-Q hanging-r))
	      #+(or)(break)))
      (set-constraints "HANGING" constraints-P constraints-Q constraints-r)
      
      ;; boundary constraints
      (multiple-value-bind (essential-P essential-Q essential-r)
	  (compute-essential-boundary-constraints
	   ansatz-space :where :surface :interface (refinement-interface h-mesh))
	;; remove restrictions on hanging-node dofs
	(let ((keys (row-keys constraints-P)))
	  (remove-keys essential-P keys)
	  (remove-keys essential-Q keys)
	  (remove-keys essential-r keys))
	;; then combine constraints
	(multiple-value-setq (constraints-P constraints-Q constraints-r)
	  (combined-constraints constraints-P constraints-Q constraints-r
				essential-P essential-Q essential-r))
	#+(or)(break)
	(set-constraints "ESSENTIAL" essential-P essential-Q essential-r))
    
      ;; eliminate the constraints
      (let ((new-P (extend-by-identity constraints-P (column-table constraints-Q) :copy t))
	    (new-Q (extend-by-identity constraints-Q (column-table constraints-Q) :copy t))
	    (new-r constraints-r))
	;; the constraints are incorporated into the ansatz-space
	(set-constraints "CONSTRAINTS" constraints-P constraints-Q constraints-r)
	(set-constraints "IP-CONSTRAINTS" new-P new-Q new-r))
      nil)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Elimination of constraints
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun eliminate-constraints (mat rhs constraints-P constraints-Q constraints-r
			      &key assemble-locally include-constraints)
  "Constraints are given by an equation:  P x = Q x + r

Here x \in V, P is an orthogonal projection on a subspace V_P of V, Q maps
some other space which may have nonempty intersection with V_P to V_P.
With S we denote the mapping Id-P.  This function returns the matrix for a
Galerkin method on the constrained space.  It is used for treating hanging
nodes and essential boundary conditions.  When assemble-locally is t the
sparse structure of mat is used instead of the sparse structure of the
hanging node interface.  When include-constraints is non-nil, the
constraints are included in matrix and rhs."
  (let ((result-mat (and mat (make-analog mat)))
	(result-rhs (and rhs (make-analog rhs)))
	(region (make-hash-table)))
    ;; determine local region of constraint
    (for-each-row-key
     #'(lambda (row-key)
	 (when (matrix-row (if assemble-locally constraints-P mat) row-key)
	   (setf (gethash row-key region) t)))
     (if assemble-locally mat constraints-P))
    ;; for space efficiency, we use a shallow copy of matrix and rhs which
    ;; is deepened in the neighborhood of the constraint.
    (when mat
      (setf (row-table result-mat) (copy-hash-table (row-table mat)))
      (setf (column-table result-mat) (copy-hash-table (column-table mat)))
      (let ((deepen (make-hash-table)))
	(dohash (row-key region)
	  (for-each-key-in-row #'(lambda (col-key) (setf (gethash col-key deepen) t))
			       result-mat row-key)
	  (for-each-key-in-col #'(lambda (row-key) (setf (gethash row-key deepen) t))
			       result-mat row-key))
	(loop for key being each hash-key in deepen do
	      (whereas ((row (matrix-row mat key)))
		(setf (matrix-row result-mat key)
		      (map-hash-table #'(lambda (key val) (values key (copy val)))
				      row)))
	      (whereas ((column (matrix-column mat key)))
		(setf (matrix-column result-mat key)
		      (map-hash-table #'(lambda (key val) (values key (copy val)))
				      column))))))
    (when rhs
      (setf (slot-value result-rhs 'fl.algebra::blocks)
	    (copy-hash-table (slot-value rhs 'fl.algebra::blocks)))
      (dohash (key region)
	(setf (vref result-rhs key) (copy (vref rhs key)))))
    
    ;; do modifications on result-mat: up to now this works only for the
    ;; special case, when the range and the domain of constraints-Q are
    ;; disjoint
    (when mat
      (let* ((Qt*A (sparse-m* constraints-Q mat :job :tn
			      :sparsity (if assemble-locally :B :A)))
	     (Qt*A*Q (sparse-m* Qt*A constraints-Q :job :nn
				:sparsity (if assemble-locally :A :B)))
	     (Qt*A*S (remove-projection-range Qt*A constraints-P :column-p t))
	     (S*A (remove-projection-range result-mat constraints-P :row-p t))
	     (S*A*Q (sparse-m* S*A constraints-Q
			       :sparsity (if assemble-locally :A :B))))
	#+(or)(break)
	;; result-mat = S*A*S
	(remove-projection-range result-mat constraints-P :column-p t)
	(m+! Qt*A*Q result-mat)
	(m+! Qt*A*S result-mat)
	(m+! S*A*Q result-mat)))
    
    ;; do modifications on result-rhs
    (when rhs
      (gemm! -1.0 mat constraints-r 1.0 result-rhs)
      (let ((Qt*rhs (sparse-m* constraints-Q result-rhs :job :tn
			       :sparsity (if assemble-locally :B :A))))
	(remove-projection-range result-rhs constraints-P)
	(m+! Qt*rhs result-rhs)))
    
    (when include-constraints
      (dohash (key region)
	(when mat
	  (m+! (mref constraints-P key key) (mref result-mat key key))
	  (for-each-key-in-row
	   #'(lambda (ck)
	       (m-! (mref constraints-Q key ck) (mref result-mat key ck)))
	   constraints-Q key))
	(when rhs
	  (m+! (vref constraints-r key) (vref result-rhs key)))))

    ;; and return the result
    (values result-mat result-rhs)))
