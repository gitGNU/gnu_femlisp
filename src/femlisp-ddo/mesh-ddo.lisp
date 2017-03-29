(in-package :femlisp-ddo)

(defmethod skel-map :around (func skel)
  (lret ((result (call-next-method)))
    (when (get-property skel :distributed-p)
      (setf (get-property result :distributed-p) t))))

(add-hook 'fl.mesh::skeleton-substance
          'assert-not-on-dd-interface
          (lambda (table)
            (dohash (cell table)
              (assert (not (distributed-p cell))))))

(add-hook 'fl.mesh::substance-boundary-cells
          'assert-non-distributed-substance-boundary
          (lambda (table)
            ;; check that all distributed substance faces have a remote neighbor
            (synchronize 'fl.mesh::substance-boundary-cells)
            (dohash (cell table)
              (when (distributed-p cell)
                (insert-into-changed cell)))
            (let ((*synchronization-merger*
                    (lambda (object id-value-pairs)
                      (let ((nr-neighbors (length id-value-pairs)))
                        (assert (= (length id-value-pairs) 2) ()
                                "Object ~A has ~D distributed neighbors"
                                object nr-neighbors)))))
              (synchronize 'fl.mesh::substance-boundary-cells))
            ;; remove those distributed faces
            (dohash (cell table)
              (when (distributed-p cell)
                (remhash cell table)))))

(defun mesh-graph (mesh)
  ;; we have to take care that the graph does
  ;; not depend at all on the internal cell ordering
  (let ((cells-of-highest-dim
          (safe-sort
           (cells-of-highest-dim mesh)
           (compare-lexicographically)
           :key #'midpoint))
        (table (make-hash-table)))
    (loop for cell in cells-of-highest-dim do
      (loop for side across (boundary cell)
            for key = (or (cell-identification side mesh) side)
            do (push cell (gethash key table))))
    (let ((edges (loop for key being each hash-key of table using (hash-value cells)
                       when (= (length cells) 2)
                         collect (append cells (list 1 key)))))
      (flet ((edge->vec (edge)
               (vector (position (edge-from edge) cells-of-highest-dim)
                       (position (edge-to edge) cells-of-highest-dim))))
        (sort edges (compare-lexicographically) :key #'edge->vec)))))

(defun distribute-mesh (mesh nr-workers)
  (let ((*graph* (mesh-graph mesh))
        (distribution (make-hash-table)))
    ;; fill distribution according to partition of mesh graph
    (loop for part in (partition-graph nr-workers)
          and k from 0 do
            (loop for cell in part do
              (loop for subcell across (subcells cell) do
                (pushnew k (gethash subcell distribution)))))
    ;; We do the same loop again for ensuring a correct order of the
    ;; identified cells.  In this loop we generate the distributed-objects
    ;; and drop those objects which are present only on other processors
    (loop for part in (partition-graph nr-workers)
          and k from 0 do
            (loop for cell in part do
              (loop for subcell across (subcells cell) do
                (unless (distributed-p subcell)
                  (let ((procs (reduce #'union (identified-cells subcell mesh)
                                       :key (rcurry #'gethash distribution)
                                       :initial-value ())))
                    (cond
                      ((member (mpi-comm-rank) procs)
                       (unless (single? procs)
                         (mpi-dbg :partition-mesh "Distributing ~A to ~A" subcell procs)
                         (ddo:ensure-distributed-class (class-of subcell))
                         (ddo:make-distributed-object subcell procs)))
                      (t (mpi-dbg :partition-mesh "Removing ~A" subcell)
                         (remhash subcell (etable mesh (dimension subcell)))))))))))
  (check mesh)
  ;; propagate changes
  (synchronize 'distribute-mesh)
  (setf (get-property mesh :distributed-p) t)
  ;; return a reasonable value
  mesh
  )

(defvar *distribute-n* nil
  "Number of workers to which the mesh should be distributed.")

(defmethod make-mesh-from
    :around ((domain <domain>) &key &allow-other-keys)
  "Distribute mesh if we are inside an MPI distributed calculation."
  (lret ((mesh (call-next-method)))
    (when (and (mpi-initialized) *distribute-n*)
      (assert (<= *distribute-n* (mpi-comm-size)))
      (distribute-mesh mesh *distribute-n*))))

(defmethod fl.mesh::do-refinement! :before
    ((skel <skeleton>) (refined-skel <skeleton>) (task list) &key refined-region)
  "Refines the distributed region first."
  (when (get-property skel :distributed-p)
    (setf (get-property refined-skel :distributed-p) t)
    (let ((distributed-task
            (remove-if-not #'distributed-p task :key #'car)))
      (dbg :ddo-refine "~D: Distributed cells to be refined: ~:{~A belonging to ~A~%~}"
           (mpi-comm-rank)
           (loop repeat 10
                 for (cell . nil) in distributed-task
                 collect (list cell (owners cell))))
      (when distributed-task
        (let ((tree (make-binary-tree
                     :red-black (list-comparison '(0 1))
                     :key (_ (let ((cell (car _)))
                               (list (dimension cell) (local-id cell)))))))
          (dbg :ddo-refine "Inserting in tree: ~D work units" (length distributed-task))
          (loop for work in distributed-task do
            (insert work tree))
          (dbg :ddo-refine "Working on tree")
          (dotree (work tree)
            (destructuring-bind (cell . rule) work
              (fl.mesh::refine-cell! rule cell skel refined-skel refined-region)
              ;; distribute children in the same way as the parent cell
              (let ((procs (owners cell)))
                (assert (> (length procs) 1))
                (for-each (lambda (child)
                            (make-distributed-object child procs))
                          (children cell skel))))))
        (dbg :ddo-refine "Synchronizing")
        (synchronize 'fl.mesh::do-refinement!)))))

