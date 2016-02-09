(in-package :fl.discretization)

(defun assembly-heap (cells)
  "Build up a dependency graph for each cell and makes a parallel heap for this graph.
The dependency graph is such that assembly of some specific cell C
blocks all subcells of C together with all neighbors
sharing a subcell with C."
  (let* ((subcell-graph (extract-complete-subgraph
                         cells
                         (_ (rest (coerce (subcells _) 'list)))))
         (supercell-graph (invert-graph subcell-graph)))
    (let ((dependents (make-hash-table)))
      (dodic ((cell subcells) subcell-graph)
        (setf (gethash cell dependents)
              (lret ((blocked subcells))
                (loop for cell2 in subcells do
                  (setf blocked (union blocked (dic-ref supercell-graph cell2)))))))
      (make-parallel-heap cells (rcurry #'gethash dependents)))))

;;;; Testing

(defun test-assembly-heap ()
  (let* ((dim 2) (level 1)
         (mesh (uniformly-refined-mesh (n-cube-domain dim) level))
         (cells (cells-of-dim mesh nil)))
    (let ((heap (assembly-heap cells)))
      (loop for cell = (take-object heap) while cell do
        (print cell)
        (drop-object cell heap))))
  (dbg-off :parallel-heap)

  #+(or) (fl.parallel::new-kernel)
  (time
   (let* (;lparallel:*kernel*
          (dim 2) (level 3)
          (mesh (uniformly-refined-mesh (n-cube-domain dim) level))
          (cells (cells-of-dim mesh dim)))
     (process-in-parallel (assembly-heap cells)
         (cell)
       (sleep 0.1)
       (print cell))))
  )
