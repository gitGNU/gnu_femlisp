;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; sparseas.lisp - ansatz spaces and operators between them
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

(in-package :discretization)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; ansatz-space
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass <ansatz-space> ()
  ((fe-class :reader fe-class :initarg :fe-class :type <fe-discretization>)
   (problem :reader problem :initarg :problem :type <problem>)
   (mesh :reader mesh :initarg :mesh :type <mesh>)
   (structure-information :accessor structure-information :initform () :type list))
  (:documentation "An <ansatz-space> is determined by finite element
discretization and problem.  The problem determines constraints like
essential boundary conditions or hanging nodes.  Those constraints are
determined during assembly and stored into the structure-information slot.

This class is used as a mixin for deriving the classes
<ansatz-space-vector> and <ansatz-space-matrix> from <sparse-vector> and
<sparse-matrix>, respectively."))

(defmethod** hierarchical-mesh ((as <ansatz-space>))
  "h-mesh accessor for ansatz-space.  Use it for emphasizing that you work
with a hierarchical mesh."
  (the <hierarchical-mesh> (mesh as)))

(defun make-fe-ansatz-space (fe-class problem mesh)
  "<ansatz-space> constructor."
  (make-instance '<ansatz-space> :fe-class fe-class
		 :problem problem :mesh mesh))

(defgeneric set-constraints (ansatz-space)
  (:documentation "Computes the constraint matrices for this ansatz-space.
Constraints arise partially because of the discretization, e.g. hanging
nodes, and partially because of essential boundary conditions.  Of course,
these matrices change when mesh or discretization are adapted."))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Linear algebra interface 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; The key->size and print-key methods compute the necessary information
;;; for the <sparse-vector> and <sparse-matrix> classes from the
;;; ansatz-space.  Element identification is respected as a property of the
;;; underlying mesh.

(definline cell-key (cell mesh)
  "If cell is identified, its identification is the key."
  (or (cell-identification cell mesh) cell))

(definline representative (obj)
  "Gets a representative cell for some key."
  (if (listp obj) (car obj) obj))

(defmethod** key->size (ansatz-space)
  (let ((fe-class (fe-class ansatz-space)))
    #'(lambda (key)
	(nr-of-inner-dofs
	 (get-fe fe-class (representative key))))))

(defmethod** print-key (ansatz-space)
  (declare (ignore ansatz-space))
  #'(lambda (key)
      (if (vertex? key)
	  (print-vertex key t t)
	  (princ key))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; ansatz-space objects
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass <ansatz-space-object> ()
  ((ansatz-space :reader ansatz-space :initarg :ansatz-space :type <ansatz-space>)
   (discretization-info
    :type list :initform () :accessor discretization-info :documentation
    "Contains additional information obtained when discretizing."))
  (:documentation "Mixin for objects to which an ansatz-space is associated."))

(defmethod mesh ((aso <ansatz-space-object>))
  (mesh (ansatz-space aso)))

(defmethod hierarchical-mesh ((aso <ansatz-space-object>))
  (hierarchical-mesh (ansatz-space aso)))

(defmethod fe-class ((aso <ansatz-space-object>))
  (fe-class (ansatz-space aso)))

(defmethod problem ((aso <ansatz-space-object>))
  (problem (ansatz-space aso)))

(defmethod surface-cells ((aso <ansatz-space-object>))
  (surface-cells-of-highest-dim (mesh aso)))

(defmethod make-analog ((aso <ansatz-space-object>))
  (make-instance (class-of aso) :ansatz-space (ansatz-space aso)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; ansatz-space vectors and matrices
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; These are special sparse vectors and matrices which carry information
;;; about their ansatz-space with them.

(defclass <ansatz-space-vector> (<ansatz-space-object> <sparse-vector>)
  ()
  (:documentation "A sparse vector which is interpreted as the ansatz-space
for a specific fe-class on a given mesh."))

(defmethod initialize-instance :after ((asv <ansatz-space-vector>) &key &allow-other-keys)
  (let ((as (ansatz-space asv)))
    (setf (slot-value asv 'key->size) (key->size as))
    (setf (slot-value asv 'print-key) (print-key as))
    (setf (slot-value asv 'multiplicity) (slot-value (problem as) 'multiplicity))))

(defun make-ansatz-space-vector (as)
  "Deprecated."
  (make-instance '<ansatz-space-vector> :ansatz-space as))

(defun sparse->ansatz-space-vector (svec &key ansatz-space)
  (let ((asv (make-ansatz-space-vector ansatz-space)))
    ;; One could do a check here, if all keys are present in the mesh.
    ;; I'll do that when the first error occurs.
    (setf (slot-value asv 'algebra::blocks)
	  (slot-value svec 'algebra::blocks))
    asv))

(defmethod component ((asv <ansatz-space-vector>) i)
  "Returns an ansatz-space-vector which denotes a component of asv.  The
vector shares part of the values."
  (let* ((comp-fedisc (component (fe-class asv) i))
	 (comp-as (make-fe-ansatz-space comp-fedisc (problem asv) (mesh asv)))
	 (comp-asv (make-instance '<ansatz-space-vector>
				  :ansatz-space comp-as))
	 (indices (coerce (range< 0 (multiplicity asv)) 'vector))
	 (cell->fe (slot-value comp-fedisc 'cell->fe)))
    (for-each-key-and-entry
     #'(lambda (key entry)
	 (let ((fe (funcall cell->fe (representative key))))
	   (when (plusp (nr-of-inner-dofs fe))
	     (setf (vec-ref comp-asv key)
		   (make-instance '<submatrix> :matrix entry
				  :row-keys (inner-dof-indices fe)
				  :col-keys indices)))))
     asv)
    comp-asv))
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; <ansatz-space-morphism>
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass <ansatz-space-morphism> (<sparse-matrix>)
  ((domain-ansatz-space :reader domain-ansatz-space
			:initarg :domain-ansatz-space :type <ansatz-space>)
   (image-ansatz-space :reader image-ansatz-space
		       :initarg :image-ansatz-space :type <ansatz-space>))
  (:documentation "A sparse-matrix which is interpreted as an morphism
between two ansatz-spaces."))

(defun make-ansatz-space-morphism (domain-as image-as)
  (flet ((keys->pattern (row-key col-key)
	   (full-crs-pattern (funcall (key->size image-as) row-key)
			     (funcall (key->size domain-as) col-key))))
    (make-instance
     '<ansatz-space-morphism>
     :domain-ansatz-space domain-as :image-ansatz-space image-as
     :row-key->size (key->size image-as) :print-row-key (print-key image-as)
     :col-key->size (key->size domain-as) :print-col-key (print-key domain-as)
     :keys->pattern #'keys->pattern)))

(defmethod m* ((asm <ansatz-space-morphism>) (asv <ansatz-space-vector>))
  (sparse->ansatz-space-vector (call-next-method)
			       :ansatz-space (image-ansatz-space asm)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; <ansatz-space-automorphism>
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass <ansatz-space-automorphism> (<ansatz-space-object> <sparse-matrix>)
  ()
  (:documentation "A sparse-matrix which is interpreted as an automorphism
for an ansatz-space.  Should probably be made a specialization of
ansatz-space-morphism."))

(defmethod initialize-instance :after ((asa <ansatz-space-automorphism>) &key &allow-other-keys)
  (let* ((as (ansatz-space asa))
	 (key->size (key->size as)))
    (with-slots (ansatz-space row-key->size col-key->size print-row-key print-col-key
			      keys->pattern)
      asa
      (setf row-key->size key->size       col-key->size key->size
	    print-row-key (print-key as)  print-col-key (print-key as)
	    keys->pattern #'(lambda (row-key col-key)
			      (full-crs-pattern (funcall key->size row-key)
						(funcall key->size col-key)))))))

(defun make-ansatz-space-automorphism (as)
  (make-instance '<ansatz-space-automorphism> :ansatz-space as))

(defun sparse->ansatz-space-automorphism (smat &key ansatz-space)
  (declare (type <sparse-matrix> smat))
  (let ((asa (make-ansatz-space-automorphism ansatz-space)))
    ;; One could do a check here, if all keys are present in mesh.  I'll do
    ;; that when the first error occurs.
    (setf (row-table asa) (row-table smat))
    (setf (column-table asa) (column-table smat))
    asa))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; multilevel-ansatz-space vectors and matrices
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; It might be that the following extraction is too costly to be
;;; performed often.  In that case one could consider incorporating it
;;; in a class derived from <ansatz-space-vector>.

(defmethod extract-level ((asv <ansatz-space-vector>) (level fixnum))
  (let* ((sub-vec (make-ansatz-space-vector (ansatz-space asv)))
	 (vblocks (slot-value sub-vec 'algebra::blocks))
	 (level-skel (cells-on-level (hierarchical-mesh asv) level)))
    (for-each-key-and-entry
     #'(lambda (key vblock)
	 (when (member-of-skeleton? (representative key) level-skel)
	   (setf (gethash key vblocks) vblock)))
     asv)
    sub-vec))

(defmethod extended-extract ((asa <ansatz-space-automorphism>) (skel <skeleton>)
			     &key (row? t) (col? t))
  "Extracts a sub-matrix from a sparse matrix.  This routine could be
accelerated by taking member-checks out of the loop."
  (let ((sub-mat (make-ansatz-space-automorphism (ansatz-space asa))))
    (for-each-key-and-entry
     #'(lambda (row-key col-key entry)
	 (when (and (or (not row?) (member-of-skeleton? (representative row-key) skel))
		    (or (not col?) (member-of-skeleton? (representative col-key) skel)))
	   (setf (mat-ref sub-mat row-key col-key) entry)))
     asa)
    sub-mat))

(defmethod extract-level ((asa <ansatz-space-automorphism>) (level fixnum))
  (extended-extract asa (cells-on-level (hierarchical-mesh asa) level)))

(defmethod decompose (as-obj)
  (let* ((nr-levels (nr-of-levels (hierarchical-mesh as-obj)))
	 (decomposed (make-array nr-levels)))
    (loop for k below nr-levels do
	  (setf (aref decomposed k) (extract-level as-obj k)))
    decomposed))

(defmethod make-row-vector-for ((A <ansatz-space-automorphism>) &optional multiplicity)
  (when multiplicity
    (assert (= multiplicity (multiplicity (problem (ansatz-space A))))))
  (make-ansatz-space-vector (ansatz-space A)))

(defmethod make-column-vector-for ((A <ansatz-space-automorphism>) &optional multiplicity)
  (when multiplicity
    (assert (= multiplicity (multiplicity (problem (ansatz-space A))))))
  (make-ansatz-space-vector (ansatz-space A)))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Interpolation
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun interpolation-matrix (ansatz-space &key level region imat)
  "The algorithm works as follows: On each cell of the provided region or
on all cells of the h-mesh a local interpolation matrix computed on the
reference finite element is collected into an interpolation matrix."
  (declare (optimize (debug 3)))
  (unless imat
    (setq imat (make-ansatz-space-automorphism ansatz-space)))
  (let* ((h-mesh (hierarchical-mesh ansatz-space))
	 (level-mesh (if level (cells-on-level h-mesh level) h-mesh))
	 (fe-class (fe-class ansatz-space)))
    (flet ((insert-local-imat (cell)
	     (when (refined-p cell h-mesh)
	       (let* ((refcell (reference-cell cell))
		      (subcells (subcells cell))
		      (children (children cell h-mesh))
		      (local-imat (local-imatrix fe-class refcell)))
		 (dotensor ((i j . mblock) local-imat)
		   (let* ((child (aref children i))
			  (child-id (cell-key (aref children i) h-mesh))
			  (subcell-id (cell-key (aref subcells j) h-mesh)))
		     (when (eq (representative child-id) child)
		       ;; we have to add only one sample for identified cells
		       (m+! mblock (mat-ref imat child-id subcell-id)))))))))
      (skel-for-each-cell #'insert-local-imat (or region level-mesh)))
    imat))

(defun constrained-interpolation-matrix (ansatz-space &key level where imat)
  "The multigrid algorithm needs an interpolation which satisfies the
constraints like essential or periodic boundary conditions."
  (let ((imat (interpolation-matrix ansatz-space :level level :imat imat)))
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
;;;; Projection
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun projection-matrix (ansatz-space &key level region pmat)
  "The algorithm works as follows: On each cell of the provided cell list
or the whole refinement a local projection matrix computed on the reference
finite element is copied into the global projection matrix."
  (unless pmat
    (setq pmat (make-ansatz-space-automorphism ansatz-space)))
  (let* ((h-mesh (hierarchical-mesh ansatz-space))
	 (level-mesh (if level (cells-on-level h-mesh level) h-mesh)))
    (flet ((insert-local-pmat (cell)
	     (when (refined-p cell h-mesh)
	       (let ((refcell (reference-cell cell))
		     (fe-class (fe-class ansatz-space)))
		 (whereas ((local-pmat (local-pmatrix fe-class refcell)))
		   (let ((subcell-children (subcell-children cell h-mesh))
			 (cell-id (cell-key cell h-mesh)))
		     (dotensor ((k . pmat-block) local-pmat)
		       (setf (mat-ref pmat cell-id (cell-key (aref subcell-children k) h-mesh))
			     pmat-block))))))))
      (skel-for-each-cell #'insert-local-pmat (or region level-mesh)))
    pmat))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Transfer between different ansatz spaces
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun c-dirichlet-dof-p (key constraints-P constraints-Q)
  (and (matrix-row constraints-P key)
       (not (matrix-row constraints-Q key))))

(defun c-slave-dof-p (key constraints-P constraints-Q)
  (when (matrix-row constraints-P key)
    (for-each-key-in-row
     #'(lambda (row-key)
	 (unless (eq key row-key)
	   (return-from c-slave-dof-p t)))
     constraints-Q key)))

(defun c-find-master (key mesh constraints-P constraints-Q)
  "Does work only for Lagrange fe.  Otherwise there may be several
masters."
  (if (c-slave-dof-p key constraints-P constraints-Q)
      (whereas ((parent (parent (representative key) mesh)))
	(c-find-master (cell-key parent mesh) mesh constraints-P constraints-Q))
      key))

(defun transfer-matrix (domain-as image-as &key no-slaves)
  "Builds a transfer matrix from domain-as to image-as."
  (let* ((domain-constraints-P (getf (structure-information domain-as) :ip-constraints-P))
	 (domain-constraints-Q (getf (structure-information domain-as) :ip-constraints-Q))
	 (image-constraints-P (getf (structure-information image-as) :ip-constraints-P))
	 (image-constraints-Q (getf (structure-information image-as) :ip-constraints-Q))
	 (mesh (mesh domain-as))
	 (interface (refinement-interface mesh))
	 (domain-disc (fe-class domain-as))
	 (image-disc (fe-class image-as))
	 (prol (make-ansatz-space-morphism domain-as image-as)))
    (assert (eq mesh (mesh image-as)))
    ;;(assert (symmetric-p mat :threshold 1.0e-7 :output t))  ;;!!
    ;; first pass: interpolate on each surface cell
    (doskel (cell mesh)
      (when (or (not (refined-p cell mesh))
		(member-of-skeleton? cell interface))
	(let ((key (cell-key cell mesh)))
	  (when (plusp (funcall (key->size image-as) key))
	    (unless (matrix-row prol key) ; don't handle twice for identified bc
	      (unless (c-dirichlet-dof-p key image-constraints-P image-constraints-Q)
		(let* ((image-fe (get-fe image-disc cell))
		       (domain-fe (get-fe domain-disc cell))
		       (local-prol (local-transfer-matrix domain-fe image-fe))
		       (subcells (subcells cell))
		       (nrows (nr-of-inner-dofs image-fe)))
		  (loop for subcell across subcells
			for subcell-ndofs across (subcell-ndofs domain-fe)
			and l = 0 then (+ l subcell-ndofs)
			unless (zerop subcell-ndofs) do
			(let ((subcell-key (cell-key subcell mesh)))
			  (unless (c-dirichlet-dof-p subcell-key image-constraints-P image-constraints-Q)
			    (m+! (matrix-slice
				   local-prol
				   :from-row 0 :from-col l
				   :nrows nrows :ncols subcell-ndofs)
				 (mat-ref prol key subcell-key))))))))))))
    ;; next pass: eliminate hanging nodes
    (for-each-col-key
     #'(lambda (j)
	 (when (c-slave-dof-p j domain-constraints-P domain-constraints-Q)
	   #+debug (format t "Handling problem key ~A slave=~A~%" j slave-p)
	   (for-each-key-and-entry-in-col
	    #'(lambda (i prol_ij)
		;; eliminate prol_ij (replace j in terms of its representation in domain-constraints-Q)
		#+debug (format t "   column ~A~%" i)
		(for-each-key-and-entry-in-row
		 #'(lambda (k Q_jk)
		     (assert (not (eq j k)))
		     #+debug (format t "     substituting ~A ~A~%" k Q_jk)
		     (gemm! 1.0d0 prol_ij Q_jk 1.0d0 (mat-ref prol i k)))
		 domain-constraints-Q j))
	    prol j)
	   (remove-column prol j)))
     prol)
    (when no-slaves
      (for-each-row-key
       #'(lambda (i)
	   (when (c-slave-dof-p i image-constraints-P image-constraints-Q)
	     (remove-row prol i)))
       prol))
    ;; return the result
    prol))


;;;; Testing: (test-sparseas)
(defun test-sparseas ()
  ;; tests if projection and interpolation with third order finite elements
  ;; yield identity on a quadratic polynomial
  (let* ((dim 1) (domain (n-cube-domain dim))
	 (mesh (uniformly-refined-hierarchical-mesh domain 1))
	 (dummy-problem (make-instance '<problem> :domain domain))
	 (ansatz-space (make-fe-ansatz-space (lagrange-fe 3) dummy-problem mesh))
	 (x (make-ansatz-space-vector ansatz-space))
	 (I (interpolation-matrix ansatz-space))
	 (P (projection-matrix ansatz-space)))
    (set-lagrange-ansatz-space-vector x #'(lambda (coord) #I"coord[0]*(1.0d0-coord[0])"))
    (let ((y (sparse-m* P x :sparsity :B)))
      (let ((z (sparse-m* I y :sparsity :B)))
	(print-svec x)
	(print-svec z)
	(assert (< (norm (m- z x)) 1.0d-10)))))
  (let* ((dim 1) (domain (n-cell-domain dim))
	 (mesh (uniformly-refined-hierarchical-mesh domain 1))
	 (dummy-problem (make-instance '<problem> :domain domain))
	 (ansatz-space (make-fe-ansatz-space (lagrange-fe 2 :nr-comps 2)
					     dummy-problem mesh))
	 (I (interpolation-matrix ansatz-space))
	 (P (projection-matrix ansatz-space)))
    (assert (midentity-p (sparse-m* P I) 1.0e-10)))
  (let* ((dim 1) (domain (n-cell-domain dim))
	 (mesh (uniformly-refined-hierarchical-mesh domain 1))
	 (dummy-problem (make-instance '<problem> :domain domain))
	 (ansatz-space (make-fe-ansatz-space (lagrange-fe 2) dummy-problem mesh))
	 (I (interpolation-matrix ansatz-space))
	 (P (projection-matrix ansatz-space)))
    (assert (midentity-p (sparse-m* P I) 1.0e-10)))
  (let* ((dim 1) (domain (n-cell-domain dim))
	 (mesh (uniformly-refined-hierarchical-mesh domain 1))
	 (dummy-problem (make-instance '<problem> :domain domain))
	 (as1 (make-fe-ansatz-space (lagrange-fe 2) dummy-problem mesh))
	 (as2 (make-fe-ansatz-space (lagrange-fe 3) dummy-problem mesh)))
    (loop repeat 1 do (refine mesh :test (rcurry #'inside-cell? (make-double-vec dim 0.25))))
    (assemble-constraints as1)
    (assemble-constraints as2)
    (let ((tm1->2 (transfer-matrix as1 as2 :no-slaves t))
	  (tm2->1 (transfer-matrix as2 as1 :no-slaves t)))
      (assert (midentity-p (sparse-m* tm2->1 tm1->2) 1.0e-10))))
  )

(tests:adjoin-femlisp-test 'test-sparseas)

