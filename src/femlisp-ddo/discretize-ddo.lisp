(in-package :ddo-femlisp)

(defmethod assemble-interior :after
    ((as <ansatz-space>) (cells list) &key rhs matrix &allow-other-keys)
  ;; ensure that also newly generated entries in rhs and matrix are distributed
  (when (and rhs
             (get-property rhs :distributed-p)
             (get-property rhs :consistent-p))
    (dbg :ddo-discretize "Handling rhs")
    (aso-make-distributed rhs :consistent-p t))
  (when (and matrix
             (get-property matrix :distributed-p)
             (get-property matrix :consistent-p))
    (dbg :ddo-discretize "Handling matrix")
    (aso-make-distributed matrix :consistent-p t)))

(defmethod increment-global-by-local-vec :around
    ((cell ddo-mixin) (svec <sparse-vector>) local-vec)
  "When the cell is shared, the corresponding rhs contribution
is calculated on several processors.  Thus its incremental effect
should be reduced appropriately."
  (call-next-method cell svec
                    (scal (/ 1.0 (length (owners cell)))
                          local-vec))
  )
  
(defmethod increment-global-by-local-mat :around
    ((smat <sparse-matrix>) local-mat
     (cell ddo-mixin) &optional (domain-cell cell))
  "When the cell is shared, the corresponding matrix contribution
is calculated on several processors.  Thus its incremental effect
should be reduced appropriately."
  (declare (ignorable domain-cell))
  (when domain-cell
    (assert (= (length (owners cell))
               (length (owners domain-cell)))))
  (call-next-method smat
                    (scal (/ 1.0 (length (owners cell)))
                          local-mat)
                    cell domain-cell)
  )
