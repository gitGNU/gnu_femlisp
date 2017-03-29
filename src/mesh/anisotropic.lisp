(in-package :fl.mesh)

(defun copy-refinement-rule (refcell &aux (refcell (reference-cell refcell)))
  "Creates a new copy refinement rule."
  (declare (optimize debug))
  (or (find-if (_ (member :copy (names _))) (refinement-rules refcell))
      (add-new-refinement-rule
       (lret* ((types (make-list (length (factor-simplices refcell)) :initial-element nil))
               (rule (fl.amop:make-programmatic-instance
                      '(refinement-rule <anisotropic-rule-mixin>)
                      :types types
                      :names (list :copy (anisotropic-name types))
                      :reference-cell refcell)))
         (with-slots (boundary-refinement-rules refinement-info) rule
           (setf boundary-refinement-rules
                 (map 'vector #'copy-refinement-rule (boundary refcell)))
           (setf refinement-info
                 (vector (cons refcell
                               (map 'vector (_ (cons (1+ _) 0))
                                    (range< 0 (length (boundary
                                                       refcell))))))))))))

(defun anisotropic-refinement-skeleton (dims types)
  (let ((factors (mapcar (lambda (dim type)
                           (lret ((skel (skeleton (n-simplex dim))))
                             (when type
                               (setf skel (refine skel)))
                             (doskel (cell skel)
                               (unless (zerop (dimension cell))
                                 (setf (skel-ref skel cell)
                                       (list :types (list type)))))))
                         dims types)))
    (reduce (rcurry #'cartesian-product
                    :property-combiner
                    (lambda (props1 props2)
                      (list :types (append (getf props1 :types)
                                           (getf props2 :types)))))
            factors)))

(defun anisotropic-boundary-refinement-rules (refcell types)
  (map 'vector 'anisotropic-refinement-rule
       (boundary refcell)
       (loop for (this . rest) on types
          and factor in (factor-simplices refcell)
          for factor-dim = (dimension factor)
          appending
            (make-list (1+ factor-dim) :initial-element
                       (append previous
                               (when (> factor-dim 1) (list this))
                               rest))
          collect this into previous)
       ))

(defun anisotropic-refinement-info (refcell types)
  (let* ((skel (skeleton refcell))
         (bdry-refined
           (refine skel :indicator
                   (lambda (subcell)
                     (and (not (eql subcell refcell))
                          (loop for side across (boundary refcell)
                                for rule across
                                         (anisotropic-boundary-refinement-rules refcell types)
                                  thereis (induced-refinement-of-subcell-refcells
                                           side rule subcell))))
                   :decouple nil))
         (product-refined (anisotropic-refinement-skeleton
                           (mapcar #'dimension (factor-simplices refcell))
                           types))
         (refp->refb (make-hash-table)))
    ;; fill refp->refb
    (doskel (cell1 bdry-refined)
      (let ((cell2
              (the t (find-cell (_ (and (eql (class-of cell1) (class-of _))
                                        (equalp (midpoint cell1) (midpoint _))))
                                product-refined))))
        (setf (gethash cell2 refp->refb) cell1)))
    ;; collect direct children
    (let ((children
            (safe-sort (find-cells (_ (not (gethash _ refp->refb))) product-refined)
                       #'< :key #'dimension))
          (subcells (subcells refcell)))
      (map 'vector
           (lambda (child)
             (cons child
                   (if (vertex-p child)
                       (vertex-position child)
                       (map 'vector
                            (lambda (child-side)
                              (let ((twin (gethash child-side refp->refb)))
                                (if twin
                                    (let ((parent (the t (get-cell-property twin bdry-refined 'PARENT))))
                                      (cons (the t (position parent subcells))
                                            (the t (position twin (children parent skel)))))
                                    (cons 0 (the t (position child-side children))))))
                            (boundary child)))))
           children))))

(defun anisotropic-refinement-rule (refcell types &aux (refcell (reference-cell refcell)))
  "Create an anisotropic refinement rule."
  (assert (= (length (factor-simplices refcell)) (length types)))
  (when (every #'null types)
    (return-from anisotropic-refinement-rule
      (copy-refinement-rule refcell)))
  (or (find-if (_ (and (anisotropic-rule-p _)
                       (equal types (slot-value _ 'types))))
               (refinement-rules refcell))
      (add-new-refinement-rule
       (lret ((rule (fl.amop:make-programmatic-instance
                     '(refinement-rule <anisotropic-rule-mixin>)
                     :types types
                     :names (list (anisotropic-name types))
                     :reference-cell refcell)))
         (with-slots (boundary-refinement-rules refinement-info) rule
           (setf boundary-refinement-rules
                 (anisotropic-boundary-refinement-rules refcell types))
           (setf refinement-info
                 (anisotropic-refinement-info refcell types)))))))

;;; Initialize some rules
(let ((cell (n-cube 2)))
  (assert (eq (copy-refinement-rule (n-cube 2))
              (copy-refinement-rule (n-cube 2))))
  (anisotropic-refinement-rule (n-cube 2) '(t nil))
  (assert (eq (anisotropic-refinement-rule (n-cube 2) '(t nil))
              (get-refinement-rule cell :anisotropic-tf)))
  (anisotropic-refinement-rule cell '(nil t))
  (assert (eq (get-refinement-rule cell :anisotropic-ft)
              (anisotropic-refinement-rule cell '(nil t)))))

(defun projection-matrix (flags)
  (lret* ((n (length flags))
          (result (eye n)))
    (loop for flag in flags and k from 0
          unless flag do
            (setf (mref result k k) 0.0))))

(defun anisotropic-refinement-indicator (filter-matrix)
  "@arg{filter-matrix} projects orthogonally on the dimensions with
refinement."
  (when (listp filter-matrix)
    (setf filter-matrix (projection-matrix filter-matrix)))
  (lambda (cell)
    (let ((Dphi (local->Dglobal
                 cell (local-coordinates-of-midpoint cell))))
      (anisotropic-name
       (loop for dim in (factor-dimensions cell)
          and from = 0 then (+ from dim)
          collecting
          (not (mzerop (m* filter-matrix
                           (matrix-slice Dphi :from-col from :ncols dim)))))))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Tests
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(in-suite mesh-suite)

(test anisotropic
  (finishes
    (map 'list #'names (refinement-rules (n-cube 2)))
    (anisotropic-refinement-skeleton '(1 2) '(t nil))
    (anisotropic-refinement-info (n-cube 2) '(nil t))
    (refinement-rules (n-cube 2))
    (refinement-rules (n-simplex 2))
    (anisotropic-refinement-rule (n-simplex 2) '(t))
    (anisotropic-refinement-rule (n-simplex 2) '(nil))
    (copy-refinement-rule (n-simplex 2))
    
    (get-refinement-rule (n-simplex 1) :copy)
    #+(or)
    (fl.plot:plot
     (refcell-refinement-skeleton
      (n-cube 2) 1 (anisotropic-name '(t nil))))
    (map 'list #'names (anisotropic-boundary-refinement-rules (n-cube 2) '(t nil)))
    (map 'vector 'midpoint (boundary (n-cube 2)))
    
    (let* ((types '(nil t))
           (indicator (anisotropic-refinement-indicator types)))
      (loop repeat 4
            for mesh = (skeleton (n-cube 2)) then (refine mesh :indicator indicator)
            finally (return mesh)))

    (let* ((mesh (uniform-mesh-on-box-domain (n-cube-domain 2) #(2 1)))
           (ind (anisotropic-refinement-indicator (projection-matrix '(nil t)))))
      (refine mesh :indicator ind))
    
    (let ((mesh (cartesian-product (refine (skeleton (n-cube 1))) (skeleton (n-cube 1)))))
      (refine mesh :indicator (anisotropic-refinement-indicator '(t nil))))
    
    (vector-map #'names (refinement-rules *unit-prism-1-2*))
    
    (let* ((refcell *unit-prism-1-2*)
           (rule (anisotropic-refinement-rule refcell '(nil t))))
      (anisotropic-boundary-refinement-rules refcell '(nil t))
      (refcell-refinement-skeleton refcell 1 rule))
    
    (let* ((cell *unit-prism-1-2*)
           (mesh (skeleton cell))
           (ind (anisotropic-refinement-indicator
                 (diag #d(0.0 1.0 1.0)))))
      (refine mesh :indicator ind))
    
    )
  )
