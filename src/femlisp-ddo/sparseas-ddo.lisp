(in-package :ddo-femlisp)

;;; define hooks on MAKE-ANSATZ-SPACE-VECTOR

(add-hook 'fl.discretization::initialize-ansatz-space
          'copy-distributed-properties
          (lambda (as)
            (when (get-property (mesh as) :distributed-p)
              (setf (get-property as :distributed-p) t))))
(add-hook 'fl.discretization::make-analog-aso
          'copy-distributed-properties
          (lambda (aso original-aso)
            (when (get-property original-aso :distributed-p)
              (setf (get-property aso :distributed-p) t))))

#|

;;; for the moment, we refrain from making vectors automatically distributed
;;; until the overall structure has stabilized

(add-hook 'fl.discretization::initialize-ansatz-space-vector
          'copy-distributed-properties
          (lambda (asv &key distributed-p consistent-p &allow-other-keys)
            (setf (get-property asv :distributed-p) distributed-p
                  (get-property asv :consistent-p) consistent-p)))

(add-hook 'fl.discretization::initialize-ansatz-space-automorphism
          'copy-distributed-properties
          (lambda (asa &key distributed-p consistent-p &allow-other-keys)
            (setf (get-property asa :distributed-p) distributed-p
                  (get-property asa :consistent-p) consistent-p)))

(add-hook 'fl.discretization::initialize-ansatz-space-morphism
          'copy-distributed-properties
          (lambda (asm &key distributed-p consistent-p &allow-other-keys)
            (setf (get-property asm :distributed-p) distributed-p
                  (get-property asm :consistent-p) consistent-p)))

|#

(defun make-distributed-in-order
    (items distributed-slots
     &key(key-mapper #'identity) (object-mapper #'identity)
       owners-mapper)
  (let ((tree (make-binary-tree :red-black (compare-lexicographically) :key #'first)))
    (dolist (item items)
      (dbg :ddo-sparseas "Checking distributed: ~A" item)
      (whereas ((key (funcall key-mapper item)))
        (dbg :ddo-sparseas "Key=~A" key)
        (let ((object (funcall object-mapper item))
              (owners (funcall owners-mapper item)))
          (dbg :ddo-sparseas "Object=~A, owners=~A" object owners )
          (if (distributed-p object)
              (assert (set-equal owners (owners object)))
              (insert (list key object owners) tree)))))
    (dotree (entry tree)
      (dbg :ddo-sparseas "Making distributed: ~A" entry)
      (make-distributed-object (second entry) (third entry) distributed-slots))))


(defgeneric addition-merger (object id-value-pairs)
  (:method (object id-value-pairs)
      "The default method does nothing")
  (:method ((vec store-vector) id-value-pairs)
      (dbg :merger "Merging values for vec{~A} from ~A" vec id-value-pairs)
    (assert (distributed-p vec) () "Mergers should be called only on distributed objects")
    (assert (equal (distributed-slots vec) '(store)) () "~A should be '(STORE)" (distributed-slots vec))
    (reduce #'m+! (remove (mpi-comm-rank) id-value-pairs :key #'car)
            :initial-value (store vec) :key #'cadr :from-end t)))

(defgeneric aso-make-distributed (aso &key consistent-p merger)
  (:documentation "Turns aso into a distributed aso")
  (:method :around (aso &key &allow-other-keys)
    (if (mpi-initialized)
        (call-next-method)
        aso))
  (:method :after (aso &key (merger #'addition-merger) &allow-other-keys)
    (dbg :ddo-as "Synchronizing")
    #+(or)
    (dbg :ddo-as "Distributed objects before synchronize: ~A"
         (distributed-data))
    (let ((*synchronization-merger* merger))
      (synchronize))
    #+(or)
    (dbg :ddo-as "Distributed objects after synchronize: ~A" (distributed-data))
    ))

(defun key-local-id (key)
  (local-id (representative key)))
(defun key-owners (key)
  (owners (representative key)))
(defun key-distributed-p (key)
  (distributed-p (representative key)))

(defmethod aso-make-distributed ((asv <ansatz-space-vector>)
                                 &key consistent-p &allow-other-keys)
  (with-slots (ansatz-space) asv
    (when (get-property ansatz-space :distributed-p)
      (setf (get-property asv :distributed-p) t
            (get-property asv :consistent-p) consistent-p)
      (let ((new-distributed-keys ()))
        ;; mark existing distributed entries as changed and find new ones
        (dovec ((vblock key) asv)
          (dbg :ddo-sparseas "Checking key ~A" key)
          (cond
            ((distributed-p vblock)
             (insert-into-changed vblock))
            ((key-distributed-p key)
             (push key new-distributed-keys))))
        ;; generate new distributed entries in a defined order
        (dbg :ddo-sparseas "New distributed keys: ~A" new-distributed-keys)
        (make-distributed-in-order
         new-distributed-keys
         '(store)
         :key-mapper (lambda (key) (list (key-local-id key)))
         :object-mapper (lambda (key) (vref asv key))
         :owners-mapper #'key-owners)
        ))))

(defmethod aso-make-distributed ((asa <ansatz-space-automorphism>)
                                 &key consistent-p &allow-other-keys)
  (with-slots (ansatz-space) asa
    (when (get-property ansatz-space :distributed-p)
      (setf (get-property asa :distributed-p) t
            (get-property asa :consistent-p) consistent-p)
      (let ((new-distributed-keys ()))
        ;; mark existing distributed entries as changed and find new ones
        (dorows (rk asa)
          (when (key-distributed-p rk)
            (dorow (ck asa rk)
              (when (key-distributed-p ck)
                (let ((mblock (mref asa rk ck)))
                  (if (distributed-p mblock)
                      (insert-into-changed mblock)
                      (push (list rk ck) new-distributed-keys)))))))
        ;; generate new distributed entries in a defined order
        (dbg :ddo-sparseas "New distributed keys: ~A" new-distributed-keys)
        (make-distributed-in-order
         new-distributed-keys
         '(store)
         :key-mapper (lambda (key) (mapcar #'key-local-id key))
         :object-mapper (lambda (key) (apply #'mref asa key))
         :owners-mapper (lambda (key)
                          (reduce #'intersection key :key #'key-owners)))))))

(defmethod make-domain-vector-for :around ((A <ansatz-space-morphism>) &optional multiplicity)
  ;; this is provisoric: it is not yet clear to me if it is generally the case
  ;; that domain vectors (solution, corrections) should be distributed and
  ;; image vectors (rhs, residuals) should not
  (declare (ignore multiplicity))
  (lret ((v (call-next-method)))
    (when (get-property (ansatz-space v) :distributed-p)
      (loop for key in (col-keys A) do (vref v key))
      (aso-make-distributed v))))

(defmethod choose-start-vector :around ((as <ansatz-space>))
  ;; this is provisoric as well and may turn out to be not the perfect way how to make
  ;; a calculation distributed
  (lret ((v (call-next-method)))
    (let ((mesh (mesh as)))
      (doskel (cell mesh)
        (let ((key (cell-key cell mesh)))
          (when (entry-allowed-p v key)
            (vref v key)))))
    (aso-make-distributed v)))

  
;;;; for debugging

;;; convert sparse vector to a list

(defun sparse-vector-to-list (svec)
  (let ((mesh (mesh svec)))
    (loop for level from 0 upto (top-level mesh)
          for result = () do
            (doskel (cell (cells-on-level mesh level))
              (let ((key (cell-key cell mesh)))
                (when (entry-allowed-p svec key)
                  (push (list (coerce (midpoint cell) 'list)
                              (entries (vref svec key)))
                        result))))
          collect (reverse result))))

;;; compare two such lists (disregarding the order)

(defun deep-compare (l1 l2 &key atom=)
  (cond ((null l1) (null l2))
        ((atom l1)
         (and (atom l2) (funcall atom= l1 l2)))
        (t (and (consp l2)
                (and (deep-compare (car l1) (car l2) :atom= atom=)
                     (deep-compare (cdr l1) (cdr l2) :atom= atom=))))))

(defun deep-threshold-comparison (threshold)
  (lambda (l1 l2)
    (deep-compare l1 l2 :atom= (_ (mzerop (- _1 _2) threshold)))))

(defun compare-pv-lists (pvs1 pvs2)
  (lret ((p= (deep-threshold-comparison 1e-12))
         (v= (deep-threshold-comparison 1e-12))
         (analysis ()))
    (loop for level from 0
          and pvs1-level in pvs1
          and pvs2-level in pvs2 do
            (format t "~&Level ~D~%" level)
            (loop for (p1 v1) in pvs1-level
                  for flag = nil do
                    (loop for (p2 v2) in pvs2-level do
                      (when (funcall p= p1 p2)
                        (setq flag t)
                        (unless (funcall v= v1 v2)
                          (format t "Mismatch at ~A: ~A <--> ~A~%" p1 v1 v2)
                          (pushnew :mismatch analysis))))
                    (unless flag
                      (format t "Missing key ~A~%" p1)
                      (pushnew :missing analysis))))))

