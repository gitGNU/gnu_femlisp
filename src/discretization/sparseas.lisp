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

(in-package :fl.discretization)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; ansatz-space
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass <ansatz-space> (property-mixin)
  ((fe-class :reader fe-class :initarg :fe-class :type <fe-discretization> :documentation
	     "The finite element class for this ansatz space.")
   (problem :reader problem :initarg :problem :documentation
	    "The proplem for this ansatz space which determines essential constraints.")
   (mesh :reader mesh :initarg :mesh :type <mesh> :documentation
	    "The mesh for this ansatz space which determines hanging-node constraints.")
   (multiplicity :reader multiplicity :initarg :multiplicity
		 :documentation "Should be a copy of problem multiplicity."))
  (:documentation "A finite element ansatz space is determined by finite
element discretization, mesh and problem.  The constraints are stored in
the slot @var{properties}."))

(defmethod initialize-instance :after ((as <ansatz-space>) &key &allow-other-keys)
  (unless (slot-boundp as 'multiplicity)
    (setf (slot-value as 'multiplicity)
	  (if (slot-boundp as 'problem)
	      (slot-value (problem as) 'multiplicity)
	      1))))

(defmethod hierarchical-mesh ((as <ansatz-space>))
  "h-mesh accessor for ansatz-space.  Use it for emphasizing that you work
with a hierarchical mesh."
  (the <hierarchical-mesh> (mesh as)))

(defun make-fe-ansatz-space (fe-class problem mesh &optional multiplicity)
  "<ansatz-space> constructor.  Somewhat shorter than with MAKE-INSTANCE."
  (make-instance '<ansatz-space> :fe-class fe-class
		 :problem problem :mesh mesh
		 :multiplicity (or multiplicity
				   (and problem (multiplicity problem))
				   1)))

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

(defmethod key->size (ansatz-space)
  (let ((fe-class (fe-class ansatz-space)))
    #'(lambda (key)
	(nr-of-inner-dofs
	 (get-fe fe-class (representative key))))))

(defmethod print-key (ansatz-space)
  (declare (ignore ansatz-space))
  #'(lambda (key) (princ key)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; ansatz-space objects
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass <ansatz-space-object> (property-mixin)
  ((ansatz-space :reader ansatz-space :initarg :ansatz-space :type <ansatz-space>))
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
    (setf (slot-value asv 'multiplicity) (slot-value as 'multiplicity))))

(defun make-ansatz-space-vector (as)
  "Deprecated."
  (make-instance '<ansatz-space-vector> :ansatz-space as))

(defmethod component ((asv <ansatz-space-vector>) i)
  "Returns an ansatz-space-vector which denotes a component of asv.  The
vector shares part of the values."
  (let* ((comp-fedisc (component (fe-class asv) i))
	 (comp-as (make-instance '<ansatz-space> :fe-class comp-fedisc :problem (problem asv)
				 :mesh (mesh asv) :multiplicity (multiplicity asv)))
	 (comp-asv (make-ansatz-space-vector comp-as))
	 (indices (coerce (range< 0 (multiplicity asv)) 'vector))
	 (cell->fe (slot-value comp-fedisc 'cell->fe)))
    (for-each-key-and-entry
     #'(lambda (key entry)
	 (let ((fe (funcall cell->fe (representative key))))
	   (when (plusp (nr-of-inner-dofs fe))
	     (setf (vref comp-asv key)
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
  (change-class (call-next-method) '<ansatz-space-vector>
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

(defmethod extract-level ((asv <ansatz-space-vector>) level)
  (let* ((sub-vec (make-ansatz-space-vector (ansatz-space asv)))
	 (vblocks (slot-value sub-vec 'fl.algebra::blocks))
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
	   (setf (mref sub-mat row-key col-key) entry)))
     asa)
    sub-mat))

(defmethod extract-level ((asa <ansatz-space-automorphism>) level)
  (extended-extract asa (cells-on-level (hierarchical-mesh asa) level)))

(defmethod decompose (as-obj)
  (let* ((nr-levels (nr-of-levels (hierarchical-mesh as-obj)))
	 (decomposed (make-array nr-levels)))
    (loop for k below nr-levels do
	  (setf (aref decomposed k) (extract-level as-obj k)))
    decomposed))

(defmethod make-domain-vector-for ((A <ansatz-space-automorphism>) &optional multiplicity)
  (when multiplicity
    (assert (= multiplicity (multiplicity (ansatz-space A)))))
  (make-ansatz-space-vector (ansatz-space A)))

(defmethod make-image-vector-for ((A <ansatz-space-automorphism>) &optional multiplicity)
  (when multiplicity
    (assert (= multiplicity (multiplicity (ansatz-space A)))))
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
		       (m+! mblock (mref imat child-id subcell-id)))))))))
      (skel-for-each #'insert-local-imat (or region level-mesh)))
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
		       (setf (mref pmat cell-id (cell-key (aref subcell-children k) h-mesh))
			     pmat-block))))))))
      (skel-for-each #'insert-local-pmat (or region level-mesh)))
    pmat))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Transfer between different ansatz spaces
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun c-dirichlet-dof-p (key constraints-P constraints-Q)
  (and constraints-P
       (matrix-row constraints-P key)
       (not (matrix-row constraints-Q key))))

(defun c-slave-dof-p (key constraints-P constraints-Q)
  (and constraints-P
       (matrix-row constraints-P key)
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
  (let* ((domain-constraints-P (getf (properties domain-as) :ip-constraints-P))
	 (domain-constraints-Q (getf (properties domain-as) :ip-constraints-Q))
	 (image-constraints-P (getf (properties image-as) :ip-constraints-P))
	 (image-constraints-Q (getf (properties image-as) :ip-constraints-Q))
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
				 (mref prol key subcell-key))))))))))))
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
		     (gemm! 1.0 prol_ij Q_jk 1.0 (mref prol i k)))
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

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; GPS choice of solver
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmethod select-linear-solver ((asa <ansatz-space-automorphism>) blackboard)
  "Tries to select a suitable solver depending on the pde problem."
  (select-linear-solver (problem asa) blackboard))


;;;; Testing
(defun test-sparseas ()
  ;; tests if projection and interpolation with third order finite elements
  ;; yield identity on a quadratic polynomial
  (let* ((dim 1) (domain (n-cube-domain dim))
	 (mesh (uniformly-refined-hierarchical-mesh domain 1))
	 (ansatz-space (make-instance '<ansatz-space> :fe-class (lagrange-fe 3) :mesh mesh))
	 (x (make-ansatz-space-vector ansatz-space))
	 (I (interpolation-matrix ansatz-space))
	 (P (projection-matrix ansatz-space)))
    (set-lagrange-ansatz-space-vector x #'(lambda (coord) #I"coord[0]*(1.0-coord[0])"))
    (let ((y (sparse-m* P x :sparsity :B)))
      (let ((z (sparse-m* I y :sparsity :B)))
	(show x)
	(show z)
	(assert (< (norm (m- z x)) 1.0d-10))
	(symbol-package (class-name (class-of z)))
	)))
  (let* ((dim 1) (domain (n-cell-domain dim))
	 (mesh (uniformly-refined-hierarchical-mesh domain 1))
	 (ansatz-space (make-instance
			'<ansatz-space> :fe-class (lagrange-fe 2 :nr-comps 2) :mesh mesh))
	 (I (interpolation-matrix ansatz-space))
	 (P (projection-matrix ansatz-space)))
    (assert (midentity-p (sparse-m* P I) 1.0e-10)))
  (let* ((dim 1) (domain (n-cell-domain dim))
	 (mesh (uniformly-refined-hierarchical-mesh domain 1))
	 (ansatz-space (make-instance '<ansatz-space> :fe-class (lagrange-fe 2) :mesh mesh))
	 (I (interpolation-matrix ansatz-space))
	 (P (projection-matrix ansatz-space)))
    (assert (midentity-p (sparse-m* P I) 1.0e-10)))
  (let* ((dim 1) (domain (n-cell-domain dim))
	 (mesh (uniformly-refined-hierarchical-mesh domain 1))
	 (as1 (make-instance '<ansatz-space> :fe-class (lagrange-fe 2) :mesh mesh))
	 (as2 (make-instance '<ansatz-space> :fe-class (lagrange-fe 3) :mesh mesh)))
    (loop repeat 1 do (refine mesh :test (rcurry #'inside-cell? (make-double-vec dim 0.25))))
    (let ((tm1->2 (transfer-matrix as1 as2 :no-slaves t))
	  (tm2->1 (transfer-matrix as2 as1 :no-slaves t)))
      (assert (midentity-p (sparse-m* tm2->1 tm1->2) 1.0e-10))))
  )

;;; (test-sparseas)
(fl.tests:adjoin-test 'test-sparseas)

