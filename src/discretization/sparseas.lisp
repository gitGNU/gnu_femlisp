;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; sparseas.lisp - ansatz spaces and operators between them
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;
;;; Copyright (C) 2003-2006 Nicolas Neuss, University Heidelberg.
;;; Copyright (C) 2006- Nicolas Neuss, University Karlsruhe.
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
;;; Linear algebra interface 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; The key->size method computes the necessary information for the
;;; <sparse-vector> and <sparse-matrix> classes from the ansatz-space.  Element
;;; identification is respected as a property of the underlying mesh.

(inlining
 (defun cell-key (cell mesh)
   "If cell is identified, its identification is the key."
   (or (cell-identification cell mesh) cell)))

(inlining
 (defun key-cells (key mesh)
   "If cell is identified, its identification is the key."
   (identified-cells (representative key) mesh)))

(defmethod key->size (ansatz-space)
  (lambda (key)
    (nr-of-inner-dofs
     (get-fe ansatz-space (representative key)))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; ansatz-space vectors and matrices
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; These are special sparse vectors and matrices which carry information
;;; about their ansatz-space with them.

(defclass <ansatz-space-vector> (<ansatz-space-object> fl.matlisp::<ht-sparse-vector>)
  ()
  (:documentation "A sparse vector which is interpreted as the ansatz-space
for a specific fe-class on a given mesh."))

(defmethod initialize-instance :after ((asv <ansatz-space-vector>)
                                       &rest args &key multiplicity &allow-other-keys)
  (unless multiplicity
    (setf (slot-value asv 'multiplicity)
          (multiplicity (ansatz-space asv))))
  (apply #'call-hooks 'initialize-ansatz-space-vector asv args))

(defmethod key->size ((asv <ansatz-space-vector>))
  (key->size (ansatz-space asv)))

(defun make-ansatz-space-vector (as &optional multiplicity)
  (make-instance '<ansatz-space-vector>
                 :ansatz-space as
                 :multiplicity multiplicity))

(defun special-ansatz-space-vector (ansatz-space &optional (type :random) (value 1.0))
  "Returns a ansatz space vector for @arg{ansatz-space} filled with
constant or random entries.  Essential constraints are satisfied."
  (with-slots (mesh)
    ansatz-space
    (lret ((asv (make-ansatz-space-vector ansatz-space)))
      (doskel (cell mesh :where :surface)
	(let ((key (cell-key cell mesh)))
	  (when (entry-allowed-p asv key)
	    (ecase type
	      (:random (fill-random! (vref asv key) value))
	      (:constant (fill! (vref asv key) value))))))
      (assemble-constraints ansatz-space)
      (destructuring-bind (&key constraints-P constraints-r &allow-other-keys)
	  (properties ansatz-space)
	(copy! (sparse-m* constraints-P constraints-r) asv))
      )))

(defun random-ansatz-space-vector (ansatz-space)
  "Returns a ansatz space vector for @arg{ansatz-space} filled with random
entries.  Essential constraints are satisfied."
  (special-ansatz-space-vector ansatz-space :random 1.0))

(defgeneric choose-start-vector (ansatz-space)
  (:documentation "Choose a reasonable start vector for some strategy.")
  (:method (as)
    "The default method chooses a random guess for eigenvalue problems and
0 otherwise."
    (if (typep (problem as) '<evp-mixin>)
	(random-ansatz-space-vector as)
	(make-ansatz-space-vector as))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; <ansatz-space-morphism>
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; abstract interface

(defclass <ansatz-space-morphism> ()
  ()
  (:documentation "A mapping between two ansatz-spaces."))

(defgeneric domain-ansatz-space (asm))
(defgeneric image-ansatz-space (asm))

(defmethod row-key->size ((asm <ansatz-space-morphism>))
  (key->size (image-ansatz-space asm)))
(defmethod col-key->size ((asm <ansatz-space-morphism>))
  (key->size (domain-ansatz-space asm)))
(defmethod keys->pattern ((asm <ansatz-space-morphism>))
  (lambda (row-key col-key)
    (full-crs-pattern (funcall (key->size (image-ansatz-space asm)) row-key)
                      (funcall (key->size (domain-ansatz-space asm)) col-key))))

(defclass <domain-image-mixin> ()
  ((domain-ansatz-space :reader domain-ansatz-space
			:initarg :domain-ansatz-space :type <ansatz-space>)
   (image-ansatz-space :reader image-ansatz-space
		       :initarg :image-ansatz-space :type <ansatz-space>)))

(defmethod initialize-instance :after ((asm <ansatz-space-morphism>)
                                       &rest args &key &allow-other-keys)
  (apply #'call-hooks 'initialize-ansatz-space-morphism asm args))

(defun make-ansatz-space-morphism (domain-as image-as)
  (fl.amop::make-programmatic-instance
   '(<ansatz-space-morphism> <domain-image-mixin> fl.matlisp::<ht-sparse-matrix>)
   :domain-ansatz-space domain-as :image-ansatz-space image-as))

;;; Is this necessary?  Check next time...
(defmethod m* ((asm <ansatz-space-morphism>) (asv <ansatz-space-vector>))
  (change-class (call-next-method) '<ansatz-space-vector>
		:ansatz-space (image-ansatz-space asm)))

(defmethod m*-product-instance ((x <ansatz-space-morphism>) (y <ansatz-space-vector>))
  (make-image-vector-for x (multiplicity y)))
(defmethod m*-tn-product-instance ((x <ansatz-space-morphism>) (y <ansatz-space-vector>))
  (make-domain-vector-for x (multiplicity y)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; <ansatz-space-automorphism>
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass <ansatz-space-automorphism> (<ansatz-space-morphism> <ansatz-space-object>)
  ()
  (:documentation "A automorphism of an ansatz space."))

(defmethod domain-ansatz-space ((asa <ansatz-space-automorphism>))
  (ansatz-space asa))

(defmethod image-ansatz-space ((asa <ansatz-space-automorphism>))
  (ansatz-space asa))

(defmethod initialize-instance :after ((asa <ansatz-space-automorphism>)
                                       &rest args &key &allow-other-keys)
  (apply #'call-hooks 'initialize-ansatz-space-automorphism asa args))

(defun make-ansatz-space-automorphism (as)
  (fl.amop::make-programmatic-instance
   '(<ansatz-space-automorphism> fl.matlisp::<ht-sparse-matrix>)
   :ansatz-space as))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; multilevel-ansatz-space vectors and matrices
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; It might be that the following extraction is too costly to be
;;; performed often.  In that case one could consider incorporating it
;;; in a class derived from <ansatz-space-vector>.

(defgeneric extended-extract (mat keys &key row? col?)
  (:documentation "Extract a sub-matrix from a sparse matrix for the given keys."))

(defmethod extended-extract ((asa <ansatz-space-automorphism>) (skel <skeleton>)
			     &key (row? t) (col? t))
  "Extracts a sub-matrix from a sparse matrix.  This routine could be
accelerated by taking member-checks out of the loop."
  (let ((sub-mat (make-ansatz-space-automorphism (ansatz-space asa))))
    (for-each-entry-and-key
     #'(lambda (entry row-key col-key)
	 (when (and (or (not row?) (member-of-skeleton? (representative row-key) skel))
		    (or (not col?) (member-of-skeleton? (representative col-key) skel)))
	   (setf (mref sub-mat row-key col-key) entry)))
     asa)
    sub-mat))

(defmethod make-domain-vector-for ((A <ansatz-space-morphism>) &optional multiplicity)
  (let ((as (domain-ansatz-space A)))
    (make-ansatz-space-vector as multiplicity)))

(defmethod make-image-vector-for ((A <ansatz-space-morphism>) &optional multiplicity)
  (let ((as (image-ansatz-space A)))
    (make-ansatz-space-vector as multiplicity)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Interpolation
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defvar *interpolation-type* :local
  "Default value for interpolation type.  @code{:local} means that
interpolation extends only to the direct children, which is reasonable for
standard interpolation of Lagrange finite elements.  @code{:global}
interpolates from highest-dimensional cells and forms an average.")

(defun interpolation-matrix (ansatz-space &key level region imat (type *interpolation-type*))
  "On each cell of the skeleton @arg{region} or on all cells of level
@arg{level} of the mesh of @arg{ansatz-space}, a local interpolation matrix
is collected into an interpolation matrix.  @arg{type} is the interpolation
type having a default value @var{*interpolation-type*}."
  (declare (optimize (debug 3)))
  (unless imat
    (setq imat (make-ansatz-space-automorphism ansatz-space)))
  (let* ((h-mesh (hierarchical-mesh ansatz-space))
	 (region (or region
		     (if level (cells-on-level h-mesh level) h-mesh)))
	 (dim (dimension h-mesh)))
    (flet ((insert-local-imat (cell)
	     (when (refined-p cell h-mesh)
	       (let* ((subcells (subcells cell))
		      (children
		       (if (eq type :local)
			   (children cell h-mesh)
			   (subcell-children cell h-mesh)))
		      (local-imat (local-interpolation-matrix cell ansatz-space type)))
		 (dotensor ((i j . mblock) local-imat)
		   (let* ((child (aref children i))
			  (child-id (cell-key (aref children i) h-mesh))
			  (subcell-id (cell-key (aref subcells j) h-mesh)))
		     (when (eq (representative child-id) child)
		       ;; we have to add only one sample for identified cells
		       (m+! mblock (mref imat child-id subcell-id)))))))))
      (apply #'skel-for-each #'insert-local-imat region
	     (unless (eq type :local) (list :dimension dim))))
    (unless (eq type :local)
      ;; if interpolation is not local, we have to scale some rows of the
      ;; interpolation matrix by the inverse of the number of supercells
      (let ((nr-of-supercells (make-hash-table)))
	(doskel (cell region :dimension dim)
	  (dovec (child (subcell-children cell h-mesh))
	    (incf (gethash (cell-key child h-mesh) nr-of-supercells 0))))
	(for-each-entry-and-key
	 #'(lambda (entry key)
	     (scal! (/ 1.0 (gethash (car key) nr-of-supercells)) entry))
	 imat)))
    imat))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Projection
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun projection-matrix (ansatz-space &key level region pmat)
  "The algorithm works as follows: On each cell of the provided cell list
or the whole refinement a local projection matrix computed on the reference
finite element is copied into the global projection matrix."
  (ensure pmat (make-ansatz-space-automorphism ansatz-space))
  (let ((h-mesh (hierarchical-mesh ansatz-space)))
    (flet ((insert-local-pmat (cell)
	     (when (refined-p cell h-mesh)
	       (whereas ((local-pmat (local-projection-matrix cell ansatz-space)))
		 (let ((subcell-children (subcell-children cell h-mesh))
		       (cell-id (cell-key cell h-mesh)))
		   (dotensor ((k . pmat-block) local-pmat)
		     (setf (mref pmat cell-id (cell-key (aref subcell-children k) h-mesh))
			   pmat-block)))))))
      (skel-for-each #'insert-local-pmat
		     (or region
			 (if level
			     (cells-on-level h-mesh level)
			     h-mesh))))
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
	 (prol (make-ansatz-space-morphism domain-as image-as)))
    (assert (eq mesh (mesh image-as)))
    ;;(assert (msymmetric-p mat :threshold 1.0e-7 :output t))  ;;!!
    ;; first pass: interpolate on each surface cell
    (doskel (cell mesh)
      (when (or (not (refined-p cell mesh))
		(member-of-skeleton? cell interface))
	(let ((key (cell-key cell mesh)))
	  (when (plusp (funcall (key->size image-as) key))
	    (unless (matrix-row prol key) ; don't handle twice for identified bc
	      (unless (c-dirichlet-dof-p key image-constraints-P image-constraints-Q)
		(let* ((image-fe (get-fe image-as cell))
		       (domain-fe (get-fe domain-as cell))
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
  "Select a suitable solver depending on the pde problem."
  (select-linear-solver (problem asa) blackboard))

#+(or)
(defmethod getrf! ((asa <ansatz-space-automorphism>) &optional ipiv)
  (call-next-method asa (or ipiv (hierarchically-ordered-cells (mesh asa)))))


;;;; Testing
(defun test-sparseas ()
  ;; tests if projection and interpolation with third order finite elements
  ;; yield identity on a quadratic polynomial
  (let* ((dim 1) (domain (n-cube-domain dim))
	 (mesh (uniformly-refined-hierarchical-mesh domain 1))
	 (problem (make-instance '<pde-problem> :domain domain
				 :components '(u) :multiplicity 1))
	 (ansatz-space (make-fe-ansatz-space (lagrange-fe 3) problem mesh))
	 (x (make-ansatz-space-vector ansatz-space))
	 (I (interpolation-matrix ansatz-space))
	 (P (projection-matrix ansatz-space)))
    (time
     (with-workers ((lambda (k)
                      (with-mutual-exclusion (x)
                        (sleep 1.0)
                        k)))
       (work-on 1)
       (work-on 2)))
    (doskel (cell mesh) (print (vref x cell)))
    (set-lagrange-ansatz-space-vector
     x #'(lambda (coord) #I"coord[0]*(1.0-coord[0])"))
    (format t "~&########### x ###########~%")
    (show x)
    (format t "~&########### P ###########~%")
    (show P)
    (let ((y (sparse-m* P x :sparsity :B)))
      (format t "~&########### I ###########~%")
      (show I)
      (format t "~&########### y ###########~%")
      (show y)
      (let ((z (sparse-m* I y :sparsity :B)))
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
    (loop repeat 1 do
	 (refine mesh :indicator (rcurry #'inside-cell? (make-double-vec dim 0.25))))
    (let ((tm1->2 (transfer-matrix as1 as2 :no-slaves t))
	  (tm2->1 (transfer-matrix as2 as1 :no-slaves t)))
      (assert (midentity-p (sparse-m* tm2->1 tm1->2) 1.0e-10))))
  )

;;; (test-sparseas)
(fl.tests:adjoin-test 'test-sparseas)

