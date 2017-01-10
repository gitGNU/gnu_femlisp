;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; graphs.lisp - graph manipulation routines
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;
;;; Copyright (C) 2016-
;;; Nicolas Neuss, FAU Erlangen-Nuernberg
;;; All rights reserved.
;;; 
;;; Redistribution and use in source and binary forms, with or without
;;; modification, are permitted provided that the following conditions
;;; are met:
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
;;; MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
;;; IN NO EVENT SHALL THE AUTHOR, THE KARLSRUHE INSTITUTE OF TECHNOLOGY,
;;; OR OTHER CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
;;; SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
;;; LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
;;; DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
;;; THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
;;; (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
;;; OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defpackage "NET.SCIPOLIS.GRAPHS"
  (:nicknames "GRAPHS")
  (:use "COMMON-LISP"
        "FL.MACROS" "FL.UTILITIES" "FL.DEBUG")
  (:export "EDGE-FROM" "EDGE-TO" "EDGE-WEIGHT" "*GRAPH*" "ALL-NODES" "PARTITION-GRAPH"
           "REPLACE-NODES-BY-NUMBERS")
  (:documentation
   "This package provides some functions for symmetric graphs."))

(in-package :net.scipolis.graphs)

(defstruct (edge (:type list))
  "An edge in a graph is represented as a structure
internally implemented as a list."
  from to weight)

(defparameter *graph*
  '( (1 2 1) (1 2 2) (2 3 1) (2 3 2) (2 4 3) (1 5 1))
  "A graph is a list of edges.  Edges are qualified by EQ (identity in memory),
therefore several equal edges (same entries) may occur in the graph.")

;; marker strategy: when wrapped inside WITH-MARKERS code can
;;; mark edges in the graph as removed, which is respected by the
;;; routines accessing the graph.  Unmarking and checking the mark
;;; is also possible, but seldom needed.

(defvar *markers* nil
  "If not NIL, this should contain a hash-table for markers.
Is used internally when code is performed within the body of
the macro WITH-MARKERS.")

(defun mark (edge)
  (setf (gethash edge *markers*) t))
(defun unmark (edge)
  (remhash edge *markers*))
(defun markedp (edge)
  "Predicate for determining if an edge has been marked."
  (gethash edge *markers*))

(defmacro with-markers (&body body)
  `(let ((*markers* (make-hash-table :test 'eq)))
     ,@body))

;;; Note: in the following we treat only symmetric graphs for simplicity

;;; elementary interface

(defun all-edges ()
  "Returns all edges in the graph as a list."
  (if *markers*
      (loop for edge in *graph*
            unless (markedp edge) collect edge)
      *graph*))

(defun nodes (edge)
  "Returns both nodes of an edge as a two-element list."
  (list (edge-from edge) (edge-to edge)))

(defun all-nodes ()
  "Returns all nodes in a graph as a list."
  (remove-duplicates
   (mapcan #'nodes (all-edges))))

(defun edges (from &optional to)
  "When called with one argument FROM, returns all edges having FROM as a node.
Otherwise returns all edges having FROM and TO as nodes."
  (if to
      (filter (list from to) (edges from) :test #'set-equal :key #'nodes)
      (filter from (all-edges) :test #'member :key #'nodes)))

;;; higher level interface

(defun is-node-of (node edge)
  "Tests if NODE is a node of EDGE"
  (member node (nodes edge)))

(defun other-node (edge node)
  "Returns the node of EDGE which is not NODE
   or NODE if the edge is a loop."
  (let ((other (first (remove node (nodes edge)))))
    (or other node)))

(defun neighbors (node)
  (remove-duplicates
   (loop for edge in (edges node)
         collect (other-node edge node))))

(defun noneuler-nodes ()
  "Filter out nodes with noneven number of edges."
  (let ((table (make-hash-table)))
    (loop for edge in (all-edges) do
      (loop for node in (nodes edge) do
        (incf (gethash node table 0))))
    (loop for node being each hash-key of table using (hash-value nr)
          when (oddp nr)
            collect node)))

(defun check-euler ()
  "Check if graph is an Euler graph."
  (assert (null (noneuler-nodes)) ()
          "There are nodes with noneven number of edges: ~A"
          (noneuler-nodes)))

;;; Zwiebelschalen-Algorithmus

(defun layer (node)
  (let ((edge (first (edges node))))
    (cond (edge (mark edge)
                (list* node edge
                       (layer (other-node edge node))))
          (t (list node)))))

(defun find-connected-node (path)
  "Find a node in PATH which is connected to the (remaining) graph."
  (when (all-edges)
    (lret ((node (loop for (node edge) on path by #'cddr
                      when (edges node) do (return node)
                      finally (return nil))))
      (assert node () "Graph not connected"))))

(defun try-continue-path (path)
  "Tries to continue path at some node"
  (let ((node (find-connected-node path)))
    (if node
        ;; continue path at node
        (let ((pos (position node path)))
          (try-continue-path
           (append (subseq path 0 pos)
                   (layer node)
                   (subseq path (1+ pos)))))
        ;; nothing to continue
        path)))

(defun euler-path (&optional node)
  (unless node
    (setf node (edge-from (first (all-edges)))))
  (assert (member node (all-nodes)) () "No such node")
  (with-markers
    (try-continue-path (list node))))

;;; Handle Non-Euler-Graphs by augmenting them with additional paths connecting non-eulerian nodes

(defun dijkstra (from to)
  "The Dijkstra algorithm for finding the distance between two nodes"
  (let* ((queue (make-instance 'queue))
         (distances (make-hash-table :test 'eql)))
    (enqueue from queue)
    (setf (gethash from distances)
          (list 0 (list from)))
    (loop for node = (dequeue queue) while node do
      (destructuring-bind (distance path)
          (gethash node distances)
        (unless (aand (gethash to distances)
                      (>= distance (first it)))
          (dolist (edge (edges node))
            (let ((nb (other-node edge node))
                  (nb-distance (+ distance (edge-weight edge))))
              (let ((entry (gethash nb distances)))
                (unless (and entry (<= (first entry) nb-distance))
                  (setf (gethash nb distances)
                        (list nb-distance (list* nb edge path)))
                  (enqueue nb queue))))))))
    (destructuring-bind (distance path)
        (gethash to distances)
      (values distance (reverse path)))))

(defun distance-matrix (nodes)
  "Calculate a distance matrix for the given vector of nodes."
  (lret* ((n (length nodes))
          (matrix (make-array (list n n) :initial-element nil)))
    (dotimes (i n)
      (dotimes (j n)
        (setf (aref matrix i j)
              (dijkstra (elt nodes i) (elt nodes j)))))))

(defun permutations (items)
  (if (null items)
      (list ())
      (loop for elem in items
            appending
            (mapcar (lambda (perm1) (cons elem perm1))
                    (permutations (remove elem items))))))

(defun best-matching (nodes)
  "Caclulates the subdivision of nodes into pairs which
has the lowest value of total length of connecting paths."
  (let ((n (length nodes))
        (distance-matrix (distance-matrix nodes))
        best-value
        best-matching)
    (assert (<= n 5) () "This function has exponential runtime and should be used
only for very small numbers of noneuler-nodes.  Use GOOD-MATCHING instead.")
    (mapcar
     (lambda (permuted-indices)
       (let ((total-length
               (loop for (i j) on permuted-indices by #'cddr
                     sum (aref distance-matrix i j))))
         (unless (aand best-value (<= it total-length))
           (setf best-value total-length
                 best-matching permuted-indices))))
     (permutations (range< 0 n)))
    (values
     (loop for (i j) on best-matching by #'cddr
           collect (list (elt nodes i) (elt nodes j)))
     best-value)))

(defun optimize-pairs (indices distance-matrix)
  "INDICES should be a list of four indices.  The result is
either the same list or, if necessary, a permuted list (i j k l)
such that the value of distance(i,j)+distance(k,l) is lowest."
(flet ((distance (i j)
         (aref distance-matrix i j)))
  (destructuring-bind (i j k l) indices
    (let ((d1 (+ (distance i j) (distance k l)))
          (d2 (+ (distance i k) (distance j l)))
          (d3 (+ (distance i l) (distance j k))))
      (cond ((< d2 (min d1 d3))
             (list i k j l))
            ((< d3 (min d1 d2))
             (list i l j k))
            (t indices))))))

(defun good-matching (nodes)
  "Find a good matching for the given NODES which is a sequence containing n nodes.
The algorithm works by rearranging an index vector #(i1 i2 ...) originally
containing 0 ... n-1 such that ((node(i1) node(i2)) (node(i3) node(i4)) ...)
is a matching which does not change under simple rearrangements which involve only
four positions."
  (let ((n (length nodes))
        (distance-matrix (distance-matrix nodes)))
    (let ((indices (coerce (range< 0 n) 'vector)))
      ;; the indices in this vector determine the pairing because they are
      ;; reordered such that node(i1)-node(i2), node(i3)-node(i4), ...
      ;; are an optimized matching

      ;; we now try two-pair-recombinations until nothing changes
      (loop
        (let ((flag nil))
          (loop for i from 0 below n by 2 do
            (loop for j from (+ i 2) by 2 below n do
              (let* ((old (list (aref indices i)
                                (aref indices (1+ i))
                                (aref indices j)
                                (aref indices (1+ j))))
                     (new (optimize-pairs old distance-matrix)))
                (unless (eq old new)
                  (destructuring-bind (a b c d) new
                    (setf (aref indices i)      a
                          (aref indices (1+ i)) b
                          (aref indices j)      c
                          (aref indices (1+ j)) d))
                  (setf flag :changed)))))
          (unless flag
            (return-from good-matching
              (loop for i from 0 by 2 below n
                    for ind1 = (aref indices i)
                    for ind2 = (aref indices (1+ i))
                    collect (list (elt nodes ind1) (elt nodes ind2))
                      into pairs
                    sum (aref distance-matrix ind1 ind2)
                      into distance
                    finally (return (values pairs distance))))))))))

(defun new-edges-for-shortest-path (a b)
  "Generates new edges (i.e. copies of existing edges)
for a shortest path between nodes A and B."
  (loop for (nil edge) on (nth-value 1 (dijkstra a b)) by #'cddr
        when edge collect (copy-seq edge)))

(defun euler-graph ()
  (append (loop for (a b) in (good-matching (noneuler-nodes))
                appending (new-edges-for-shortest-path a b))
          *graph*))

(defun best-route (&optional node)
  (let ((*graph* (euler-graph)))
    (euler-path node)))

;;;; Graph partitioning

(defun test-and-maybe-swap (table n1 n2)
  (unless (eql n1 n2)
    (flet ((part (node)
             (gethash node table)))
      (let ((part1 (part n1))
            (part2 (part n2)))
        (unless (= part1 part2)
          (let ((c11 (count part1 (neighbors n1) :key #'part :test-not #'=))
                (c12 (count part2 (neighbors n1) :key #'part :test-not #'=))
                (c21 (count part1 (neighbors n2) :key #'part :test-not #'=))
                (c22 (count part2 (neighbors n2) :key #'part :test-not #'=)))
            (let ((delta (- (+ c11 c22) (+ c12 c21) (* 2 (length (edges n1 n2))))))
              (when (plusp delta)
                (dbg :partition-graph "~A -> ~D and ~A -> ~D (improvement ~D)"
                     n1 part2 n2 part1 delta)
                (setf (gethash n1 table) part2
                      (gethash n2 table) part1)))))))))

(defun partition-graph (n)
  (declare (type (integer 1 *) n))
  ;; first, partition the nodes of graph roughly
  (let* ((nodes (coerce (all-nodes) 'vector))
         (table (make-hash-table :test 'eql))
         (nr-nodes (length nodes)))
    (loop for node across nodes
          for i from 0 by (/ n nr-nodes) do
            (setf (gethash node table) (floor i)))
    ;; then try swapping pairs of nodes between partitions
    ;; with the purpose of decreasing interconnections
    (loop
      (let ((flag nil))
        (loop for i from 0 below nr-nodes do
          (loop for j from (1+ i) below nr-nodes do
            (when (test-and-maybe-swap table (aref nodes i) (aref nodes j))
              (setf flag t))))
        (unless flag
          (return-from partition-graph
            (group-by (_ (gethash _ table)) (all-nodes))))))))

(defun replace-nodes-by-numbers (&optional (graph *graph*) keep-id)
  (let* ((*graph* graph)
         (nodes (all-nodes)))
    (flet ((pos (node) (position node nodes)))
      (loop for (a b dist . rest) in (all-edges)
            collect (list* (pos a) (pos b) dist
                           (and keep-id rest))))))

(defun graph-tests ()
  (check-euler)
  (with-markers (layer 1))
  (distance-matrix (all-nodes))
  (good-matching (noneuler-nodes))
  (dijkstra 1 4)
  (best-route 1)
  )
