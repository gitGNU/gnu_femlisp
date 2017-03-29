;;; -*- mode: lisp; fill-column: 75; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; relations.lisp - relations
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;
;;; Copyright (C) 2015-
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

(in-package :net.scipolis.relations)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; enhancements for the trees package
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun find-first-path-satisfying (tree condition &key (direction :forward))
  "Finds the path to the DIRECTION-first node having an item
 satisfying CONDITION.
CONDITION should specify an semi-open interval wrt the trees order!"
  (multiple-value-bind (left right)
      (ecase direction
        (:forward (values #'left #'right))
        (:backward (values #'right #'left)))
    (labels ((find-path (path)
               (let* ((node (first path))
                      (testp (funcall condition (datum node))))
                 (if testp
                     (or (awhen (funcall left node)
                           (find-path (cons it path)))
                         path)
                     (awhen (funcall right node)
                       (find-path (cons it path)))))))
      (aand (root tree) (find-path (list it))))))

(defun find-first-path-geq (tree value)
  (find-first-path-satisfying
   tree
   (lambda (nodeval)
     (not (funcall (pred tree) nodeval value)))))

(defun find-last-path-leq (tree value)
  (find-first-path-satisfying
   tree
   (lambda (nodeval)
     (not (funcall (pred tree) value nodeval)))
   :direction :backward))

(defun next (path &optional (direction :forward))
  (multiple-value-bind (left right)
      (ecase direction
        (:forward (values #'left #'right))
        (:backward (values #'right #'left)))
    (labels ((next-from-subtree (path)
               (let ((node (first path)))
                 (aif (funcall left node)
                      (next-from-subtree (cons it path))
                      path)))
             (up (path)
               (let ((old-node (first path))
                     (new-path (rest path)))
                 (whereas ((new-node (first new-path)))
                   (cond ((eq (funcall left new-node) old-node)
                          new-path)
                         ((eq (funcall right new-node) old-node)
                          (up new-path))
                         (t (error "should not happen")))))))
      (when path
        (let ((old-node (first path)))
          (acond ((funcall right old-node)
                  (next-from-subtree (cons it path)))
                 (t (up path))))))))

(defun inverse-direction (dir)
  (ecase dir
    (:forward :backward)
    (:backward :forward)))

(defun for-each-in-range (worker tree &key from to (direction :forward))
  "Walks tree in the given range and calls worker on each datum."
  (let* ((from-path
           (find-first-path-satisfying
            tree (lambda (nodeval)
                   (not (ecase direction
                          (:forward (funcall (pred tree) nodeval from))
                          (:backward (funcall (pred tree) from nodeval)))))
            :direction direction))
         (to-path
           (find-first-path-satisfying
            tree (lambda (nodeval)
                   (not (ecase direction
                          (:forward (funcall (pred tree) to nodeval))
                          (:backward (funcall (pred tree) nodeval to)))))
            :direction (inverse-direction direction)))
         (path (and to-path from-path))
         (limit-node (first (next to-path direction))))
    (loop for node = (first path)
          while (and node (not (eq node limit-node)))
          do (funcall worker (datum node))
             (setq path (next path direction)))))

(defun tree-leaves (tree)
  (let ((result ()))
    (dotree (item tree)
      (push item result))
    (nreverse result)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; relation indices (binary trees working on tuples)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun list-comparison (order)
  "Returns a comparison function for comparing two lists of numbers
lexicographically, but with elements ordered according to ORDER."
  (lambda (seq1 seq2)
    (loop for k in order
          for a = (elt seq1 k)
          and b = (elt seq2 k)
          do
             (cond ((< a b)
                    (return t))
                   ((< b a)
                    (return nil))))))

(defun make-index (order)
  "Creates a binary tree for indexing the relation in the appropriate order."
  (make-binary-tree :avl (list-comparison order)
                    :key #'identity))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; relation utilities
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun new-index (old-index new-order)
  (lret ((new-index (make-index new-order)))
    (dotree (x old-index)
      (insert x new-index))))

(defun needed-order (select-items)
  "Returns an order in which all fixed variables occur first,
and also the number of fixed and variable positions."
  (loop for k from 0
        and x in select-items
        if (eq x '_)
          collect k into variable
        else
          collect k into fixed
        finally
           (return
             (values (append fixed variable)
                     fixed variable))))

(defun get-index (R order)
  (cdr (assoc order (indices R) :test 'equal)))

(defun ensure-index (R order)
  (or (get-index R order)
      (with-slots (indices) R
        (lret* ((sample-index (cdr (first indices)))
                (index (new-index sample-index order)))
          (push (cons order index)
                indices)))))

(defun extract-from-seq-in-order (seq order)
  (mapcar (_ (elt seq _)) order))

(defun from-key (select-items)
  (substitute 0 '_ select-items))

(defun to-key (select-items)
  (substitute most-positive-fixnum '_ select-items))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; relation class and methods
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass relation ()
  ((arity :reader arity :initarg :arity :type integer)
   (indices :reader indices :initarg :indices :type list
            :documentation "Association list order->index"))
  (:documentation "Relation between numbers"))

(defun make-number-relation (n)
  "Returns a relation between n numbers."
  (let* ((order (loop for k below n collect k))
         (index (make-index order)))
    (make-instance 'relation
                   :arity n
                   :indices (list (cons order index)))))

(defun R-insert (R &rest items)
  (assert (every (_ (< _ most-positive-fixnum)) items))
  (loop for (order . tree) in (indices R)
        do (insert items tree)))


(defun R-select (R &rest select-items)
  "Selects a range of tuples from the relation R."
  (let* ((order (needed-order select-items))
         (tree (ensure-index R order))
         (results ()))
    (for-each-in-range (lambda (x) (push x results))
                       tree
                       :from (from-key select-items)
                       :to (to-key select-items))
    (nreverse results)))

(defun R-some (R &rest select-items)
  "Checks if there exists some entry in R with this specification."
  (let* ((order (needed-order select-items))
         (tree (ensure-index R order)))
    (for-each-in-range (_ (return-from R-some t))
                       tree
                       :from (from-key select-items)
                       :to (to-key select-items))
    nil))

(defun R-remove (R &rest items)
  (let ((entries (apply #'R-select R items)))
    (loop for (order . tree) in (indices R) do
      (loop for entry in entries do
        (trees:delete entry tree)))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;  Testing
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun test-relations ()
  (let ((tree (make-binary-tree :avl #'< :key #'identity)))
    (loop for k below 10 do
      (insert k tree)
      (insert k tree))
    (trees:delete -1 tree)
    (trees:delete 3 tree)
    (trees:delete 3 tree)
    (trees:delete 10 tree)
    (pprint-tree tree)
    (for-each-in-range #'print tree :from 5 :to 0 :direction :backward))

  (let ((tree (make-binary-tree :avl
                                (list-comparison '(0 1))
                                :key #'identity)))
    (loop repeat 10 do
      (insert (list (random 10) (random 10)) tree))
    (loop with count = 0
          repeat 1000 do
            (for-each-in-range (_ (incf count))
                               tree :from '(3 0) :to '(7 -1))
          finally (return count)))
  
  (funcall (list-comparison '(0 1)) '(3 2) '(2 7))
  
  (let ((t1 (make-binary-tree :avl (list-comparison '(0 1)) :key #'identity)))
    (loop repeat 10 do
      (insert (list (random 10) (random 10)) t1))
    (let ((t2 (new-index t1 '(0 1))))
      (pprint-tree t2)
      (dotree (x t2)
        (print x))))

  (let ((R  (make-number-relation 3)))
    (loop repeat 1000 do
      (R-insert R (random 10) (random 10) (random 10)))
    (let ((entries (R-select R 3 '_ 4)))
      (print entries)
      (when entries
        (apply #'R-remove R (first entries)))
      (print (R-select R 3 '_ 4))
      ))

  (let ((tree (make-binary-tree :avl #'< :key #'identity))
        (n 20000))
    (loop for i below n do
      (insert i tree))
    (time
     (loop for i below n do
       (assert (= 1 (lret ((count 0))
                      (for-each-in-range (_ (incf count))
                                         tree
                                         :from i
                                         :to i)))))))

  )

;;;; (test-relations)
