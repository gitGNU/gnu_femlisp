;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; sparseif.lisp - Interface between discretization and algebra
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

;;; This file provides the interface between finite-cell
;;; discretization and sparse-matrix representation.

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; interface
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; generation of local vectors

(defun make-local-vec (ansatz-space cell)
  "Generates a local vector for local discretization."
  (let ((fe (get-fe ansatz-space cell))
	(multiplicity (multiplicity ansatz-space)))
    (map 'simple-vector
	 #'(lambda (fe) (make-real-matrix (nr-of-dofs fe) multiplicity))
	 (components fe))))

(defvar *local-mat-pool*
  (make-instance 'parpool :test 'equal)
  "Parallel pool from which local matrices may be extracted.")

(defvar *use-pool-p* nil "T if pool should be used.")
  
;;; (setq *local-mat-pool* (make-instance 'parpool :test 'equal))

(defgeneric make-local-mat (as1 cell1 &optional as2 cell2)
  (:documentation
   "Generates a local matrix discretization for the given ansatz-space(s).")
  (:method (as1 cell1 &optional (as2 as1) (cell2 cell1))
      "Default method.  Allocates a new matrix."
    (let ((fe-1 (get-fe as1 cell1))
          (fe-2 (get-fe as2 cell2)))
      (lret ((result
               (let ((comps-1 (components fe-1))
                     (comps-2 (components fe-2)))
                 (make-filled-array
                  (list (length comps-1) (length comps-2))
                  :initializer (lambda (i j)
                                 (make-real-matrix (nr-of-dofs (aref comps-1 i))
                                                   (nr-of-dofs (aref comps-2 j))))))))
        (array-for-each #'x<-0 result))))
  (:method :around (as1 cell1 &optional (as2 as1) (cell2 cell1))
    "Try using a pool if provided"
    (if (and *use-pool-p* *local-mat-pool*)
        (let* ((fe-1 (get-fe as1 cell1))
               (fe-2 (get-fe as2 cell2))
               (key (list fe-1 fe-2)))
          (with-mutual-exclusion (*local-mat-pool*)
            (aif (get-from-pool *local-mat-pool* key)
                 (progn (array-for-each #'x<-0 it)
                        it)
                 (lret ((local-mat (call-next-method)))
                   (assert local-mat)
                   (register-in-pool *local-mat-pool* key local-mat)))))
        (call-next-method))))

;;; transfer between local and global vector

(defgeneric fill-local-from-global-vec (cell global-vec local-vec)
  (:documentation "Copies the region in global-vec determined by cell to
local-vec."))

(defgeneric get-local-from-global-vec (cell global-vec)
  (:documentation "Maps the region in global-vec determined by cell to a
local vector."))

(defgeneric set-global-to-local-vec (cell global-vec local-vec)
  (:documentation "Sets the region in global-vec determined by cell to the
values of the local vector array."))

(defgeneric increment-global-by-local-vec (cell global-vec local-vec)
  (:documentation "Increments the region in global-vec determined by cell
to the values of the local vector array."))

(defgeneric global-local-vector-operation (global-vec cell local-vec operation)
  (:documentation "Performs an operation interfacing global to local values."))

;;; transfer between local and global matrix

(defgeneric get-local-from-global-mat (global-mat cell &optional domain-cell)
  (:documentation "Maps the region in the global stiffness matrix
determined by cell to a local matrix array."))

(defgeneric set-global-to-local-mat (global-mat local-mat cell &optional domain-cell)
  (:documentation "Sets the region in global-mat determined by cell to the
values of the local matrix array."))

(defgeneric fill-local-from-global-mat (global-mat local-mat cell &optional domain-cell)
  (:documentation "Copies the region in global-mat determined by cell to
local-mat."))

(defgeneric increment-global-by-local-mat (global-mat local-mat cell &optional domain-cell)
  (:documentation "Increments the region in global-mat determined by cell
to the values of local-mat."))

(defgeneric global-local-matrix-operation
    (global-mat local-mat image-cell domain-cell operation
                &key image-subcells domain-subcells)
  (:documentation "Performs some operation interfacing global to local values."))

;;; The actual implementation for <sparse-vector> and <sparse-matrix>.

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; <sparse-vector> interface
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; the following routine should be extended for extracting subvectors, if
;;; partial equations are to be assembled on subdomains

(with-memoization (:type :global)
  (defun fe-secondary-information (fe)
    "Computes for a (vector) finite element the secondary information which
is important when extracting data from ansatz-space vectors and matrices.
For each dof in the fe, the access information to the global matrix in the
form vblock/in-vblock-index and the access information to the local matrix
in the form component-index/in-component-index is computed."
    (memoizing-let ((fe fe))
      (let* ((nr-dofs (reduce #'+ (components fe) :key #'nr-of-dofs))
	     (nr-subcells (length (subcells (reference-cell fe))))
             (nr-components (nr-of-components fe))
	     (component-index (make-array nr-dofs))
	     (in-component-index (make-array nr-dofs))
	     (vblock-index (make-array nr-dofs))
	     (in-vblock-index (make-array nr-dofs))
	     (in-vblock-pos (make-array nr-subcells :initial-element 0))
	     (in-local-start (make-array (list nr-components (1+ nr-subcells)) :initial-element -1))
	     (in-global-start (make-array (list nr-components nr-subcells) :initial-element -1))
	     (k 0))
	(loop+ ((scalar-fe (components fe))
                comp)
          do
          (loop
            for dof across (fe-dofs scalar-fe)
            and in-comp from 0 do
            (let ((sci (dof-subcell-index dof)))
              (when (minusp (aref in-local-start comp sci))
                (setf (aref in-local-start comp sci)
                      in-comp)
                (setf (aref in-global-start comp sci)
                      (aref in-vblock-pos sci)))
              (setf (aref component-index k) comp
                    (aref in-component-index k) in-comp
                    (aref vblock-index k) sci
                    (aref in-vblock-index k) (aref in-vblock-pos sci)
                    k (1+ k))
              (incf (aref in-vblock-pos sci)))
            finally (setf (aref in-local-start comp nr-subcells)
                          (1+ in-comp))))
	(assert (= k nr-dofs))  ; consistency check
	(list :nr-of-components (nr-of-components fe) :nr-dofs nr-dofs
	      :component-index component-index :in-component-index in-component-index
	      :vblock-index vblock-index :in-vblock-index in-vblock-index
              :in-local-start in-local-start :in-global-start in-global-start)))))

(with-memoization (:type :global :test 'equalp)
  (defun fe-extraction-information (fe indices)
    "Computes information for extracting components out of a vector finite
 element."
    (memoizing-let ((components (components fe)) (indices indices))
      (assert (and (every (lambda (i) (< i (length components))) indices)
		   (<= (length indices) (length components))
		   (equalp indices (remove-duplicates indices))))
      (destructuring-bind (&key component-index vblock-index in-vblock-index
				&allow-other-keys)
	  (fe-secondary-information fe)
	(coerce (loop+ ((index indices)) appending
                   (loop+ ((c component-index)
                           (v vblock-index)
                           (iv in-vblock-index))
                      when (and (= c index) (zerop v))
                      collecting iv))
		'vector)))))

(defmethod global-local-vector-operation ((svec <sparse-vector>) (cell <cell>)
                                          local-vec operation)
  (let ((fe (get-fe (ansatz-space svec) cell)))
    (destructuring-bind (&key nr-dofs component-index in-component-index
                              vblock-index in-vblock-index &allow-other-keys)
        (fe-secondary-information fe)
      (let* ((mesh (mesh svec))
             (keys (vector-map (rcurry #'cell-key mesh) (subcells cell))))
        (with-mutual-exclusion (svec)
        ;;(with-region (svec keys)
          (let ((vblocks (extract-value-blocks svec keys))
                (multiplicity (multiplicity svec)))
            (dotimes (i nr-dofs)
              (let ((component-index (aref component-index i))
                    (in-component-index (aref in-component-index i))
                    (vblock-index (aref vblock-index i))
                    (in-vblock-index (aref in-vblock-index i)))
                (dotimes (j multiplicity)
                  (symbol-macrolet
                      ((local (mref (aref local-vec component-index)
                                    in-component-index j))
                       (global (mref (aref vblocks vblock-index)
                                     in-vblock-index j)))
                    (ecase operation
                      (:local<-global (setq local global))
                      (:global<-local (setq global local))
                      (:global+=local (incf global local))
                      (:global-=local (decf global local)))))))))))))

(defmethod fill-local-from-global-vec ((cell <cell>) (svec <sparse-vector>) local-vec)
  (global-local-vector-operation svec cell local-vec :local<-global))

(defmethod get-local-from-global-vec ((cell <cell>) (svec <sparse-vector>))
  (lret ((local-vec (make-local-vec (ansatz-space svec) cell)))
    (fill-local-from-global-vec cell svec local-vec)))

(defmethod set-global-to-local-vec ((cell <cell>) (svec <sparse-vector>) local-vec)
  (global-local-vector-operation svec cell local-vec :global<-local))

(defmethod increment-global-by-local-vec ((cell <cell>) (svec <sparse-vector>) local-vec)
  (global-local-vector-operation svec cell local-vec :global+=local))

(defun set-lagrange-ansatz-space-vector (asv func)
  "Sets an ansatz-space-vector to interpolate a given function.  This is
still a suboptimal implementation for vector functions."
  (let ((ansatz-space (ansatz-space asv)))
    (doskel (cell (mesh asv) :where :surface :dimension :highest)
      (let ((local-vec (make-local-vec ansatz-space cell))
            (fe (get-fe ansatz-space cell)))
        (loop for comp from 0
           and comp-fe across (components fe) do
           (do-dof (dof comp-fe)
             (let ((value (funcall func (local->global cell (dof-gcoord dof)))))
               (setq value (fl.cdr::ensure-1-component-vector value))
               (dotimes (i (multiplicity asv))
                 (setf (mref (aref local-vec comp) (dof-index dof) i)
                       (mref (aref value comp) 0 i))))))
        (set-global-to-local-vec cell asv local-vec)))))

(defun multiple-evaluate-local-fe (local-vec shape-values)
  "Evaluates the vector given in @arg{local-vec} at multiple points.  Here
@arg{local-vec} should be a data vector obtained with
@function{get-local-from-global-vec} and @arg{ip-values} should be a vector
obtained from @function{ip-values}."
  (vector-map (lambda (point-values)
		(vector-map #'m*-tn point-values local-vec))
	      shape-values))

(defun extract-ip-data (cell qrule property-list)
  "Converts all ansatz-space objects in the parameters list into local
value arrays corresponding to the finite element."
  (loop for (key object) on property-list by #'cddr
	when (typep object '<ansatz-space-vector>)
	collect key and
	collect (multiple-evaluate-local-fe
		 (get-local-from-global-vec cell object)
		 (ip-values (get-fe (ansatz-space object) cell) qrule))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; <sparse-matrix> interface
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defparameter *check-lock* (bordeaux-threads:make-recursive-lock))
(defparameter *check-hash-table* (make-hash-table :test 'equal))
(defun lock (mat rk)
  (bordeaux-threads:with-recursive-lock-held (*check-lock*)
    (aif (gethash mat *check-hash-table*)
         (error "Concurrent access to matrix~%~A~%from row~%~A~%and row~%~A~%" mat it rk)
         (setf (gethash mat *check-hash-table*) rk))))
(defun unlock (mat)
  (bordeaux-threads:with-recursive-lock-held (*check-lock*)
    (remhash mat *check-hash-table*)))

(defun delay-possible-p (keys)
  (or (notany #'identification-p keys)
      (let ((ids (remove-if-not #'identification-p keys)))
        (= (length (remove-duplicates ids)) (length ids)))))

(defmethod global-local-matrix-operation ((smat <ansatz-space-automorphism>) local-mat
                                          (image-cell <cell>) (domain-cell <cell>)
                                          operation
                                          &key
                                          (image-subcells t) (domain-subcells t)
                                          (delay-p nil))
  (declare (optimize debug))
  (let ((domain-fe (get-fe (domain-ansatz-space smat) domain-cell))
	(image-fe (get-fe (image-ansatz-space smat) image-cell)))
    (destructuring-bind
        (&key
           ;; ((:nr-dofs nr-dofs-1))
           ;; ((:component-index component-index-1)) ((:in-component-index in-component-index-1))
           ;; ((:vblock-index vblock-index-1)) ((:in-vblock-index in-vblock-index-1))
           ((:in-local-start in-local-start-1)) ((:in-global-start in-global-start-1))
         &allow-other-keys)
	(fe-secondary-information image-fe)
      (destructuring-bind
          (&key
             ;; ((:nr-dofs nr-dofs-2))
             ;; ((:component-index component-index-2)) ((:in-component-index in-component-index-2))
             ;; ((:vblock-index vblock-index-2)) ((:in-vblock-index in-vblock-index-2))
             ((:in-local-start in-local-start-2)) ((:in-global-start in-global-start-2))
           &allow-other-keys)
	  (fe-secondary-information domain-fe)
	(let* ((mesh (mesh smat))
	       (row-keys
                 (if image-subcells
                     (map 'list (rcurry #'cell-key mesh) (subcells image-cell))
                     (list (cell-key image-cell mesh))))
               (col-keys
                 (if (and (eq image-cell domain-cell)
                          (eq image-subcells domain-subcells))
                     ;; shortcut: col-keys:=row-keys
                     row-keys
                     (if domain-subcells
                         (map 'list (rcurry #'cell-key mesh) (subcells domain-cell))
                         (list (cell-key domain-cell mesh))))))
          (let* ((mblocks (with-mutual-exclusion (smat)
                            (extract-value-blocks smat row-keys col-keys)))
                 (operation
                   (ecase operation
                     (:local<-global #'fl.matlisp::matop-x<-y!)
                     (:global<-local #'fl.matlisp::matop-x->y!)
                     (:global+=local #'fl.matlisp::matop-y+=x!)
                     (:global-=local #'fl.matlisp::matop-y-=x!))))
            (flet ((worker (vblock-index-1)
                     (loop
                       for vblock-index-2 below (length col-keys) do
                         (whereas ((global-block (aref mblocks vblock-index-1 vblock-index-2)))
                           (lock global-block (elt row-keys vblock-index-1))
                           (loop for comp1 below (nr-of-components image-fe) do
                             (loop for comp2 below (nr-of-components domain-fe) do
                               (whereas ((local-block (aref local-mat comp1 comp2)))
                                 (funcall
                                  operation
                                  local-block
                                  global-block
                                  (aref in-global-start-1 comp1 vblock-index-1)
                                  (aref in-global-start-2 comp2 vblock-index-2)
                                  (aref in-local-start-1 comp1 vblock-index-1)
                                  (aref in-local-start-2 comp2 vblock-index-2)
                                  (aref in-local-start-1 comp1 (1+ vblock-index-1))
                                  (aref in-local-start-2 comp2 (1+ vblock-index-2))))))
                           (unlock global-block)))))
              (values
               ;; we return a list of thunks
               ;; which can be processed in parallel
               (mapcar (lambda (k) (lambda () (worker k)))
                       (range< 0 (length row-keys)))
               ;; but also a flag indicating if parallel processing is possible
               ;; which may not be the case if the cell contains identified subcells
               (and delay-p
                    (delay-possible-p row-keys)
                    (or (eq row-keys col-keys)
                        (delay-possible-p col-keys))))
              )))))))

(defclass chunk-queue (parqueue)
  ((in-work-count :initform 0)
   (finalizer :initform nil))
  (:documentation "Queue for chunk work.  A chunk is a list of thunks,
which can be enqueued only if the queue is empty."))

(defmethod dequeue :after ((cq chunk-queue))
  (with-mutual-exclusion (cq)
    (with-slots (in-work-count) cq
      (incf in-work-count))))
  
(defun work-done (cq)
  (with-mutual-exclusion (cq)
    (with-slots (in-work-count finalizer) cq
      (decf in-work-count)
      (when (and finalizer (emptyp cq) (zerop in-work-count))
        (funcall finalizer)
        (setf finalizer nil)))))

(defgeneric work-on-queue (work-queue)
  (:method ((cq chunk-queue))
      (loop for work = (with-mutual-exclusion (cq)
                         (unless (emptyp cq)
                           (dequeue cq)))
            while work do
              (funcall work)
              (work-done cq))))

(defgeneric chunk-enqueue (cq chunk &optional finalizer)
  (:method ((cq chunk-queue) (chunk list) &optional finalizer)
      (loop until (with-mutual-exclusion (cq)
                    (with-slots (in-work-count) cq
                      (when (and (emptyp cq)
                                 (zerop in-work-count))
                        ;; we can get rid of the work chunk
                        (dolist (work chunk t)
                          (enqueue work cq))
                        (setf (slot-value cq 'finalizer) finalizer)
                        t)))
            do (work-on-queue cq))))

(defmethod fill-local-from-global-mat ((smat <sparse-matrix>) local-mat
                                       (cell <cell>) &optional (domain-cell cell))
  (global-local-matrix-operation
   smat local-mat cell domain-cell :local<-global))

(defmethod get-local-from-global-mat ((smat <sparse-matrix>) (cell <cell>)
                                      &optional (domain-cell cell))
  (lret* ((as (ansatz-space smat))
          (local-mat (make-local-mat as cell as domain-cell)))
    (fill-local-from-global-mat
     smat local-mat cell domain-cell)))

(defmethod set-global-to-local-mat ((smat <sparse-matrix>) local-mat
                                    (cell <cell>) &optional (domain-cell cell))
  (global-local-matrix-operation
   smat local-mat cell domain-cell :global<-local))

(defmethod increment-global-by-local-mat ((smat <sparse-matrix>) local-mat
                                          (cell <cell>) &optional (domain-cell cell))
  (global-local-matrix-operation
   smat local-mat cell domain-cell :global+=local))

;;;; Testing
(defun test-sparseif ()

  (let* ((fe-class (lagrange-fe 5 :nr-comps 2))
          (fe (get-fe fe-class (n-cube 2))))
         ;;(nr-of-dofs (aref (components fe) 0))
    (fe-secondary-information fe))

  (let ((cq (make-instance 'chunk-queue)))
    (chunk-enqueue cq (list (_ (print 1))  (_ (print 2))) (_ (print "fertig")))
    (work-on-queue cq))

  )

;;; (test-sparseif)
(fl.tests:adjoin-test 'test-sparseif)
