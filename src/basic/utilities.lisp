;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; utilities.lisp - Useful utility functions
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

(in-package :utilities)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Functions
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun required-argument ()
  (error "A required keyword argument was not supplied."))

(defun compose (&rest funcs)
  "Returns the composition of the given functions."
  (destructuring-bind (func1 . rest) (reverse funcs)
    #'(lambda (&rest args)
	(reduce #'(lambda (v f) (funcall f v))
		rest
		:initial-value (apply func1 args)))))

(definline curry (func &rest args)
  #'(lambda (&rest after-args)
      (apply func (append args after-args))))

(definline rcurry (func &rest args)
  #'(lambda (&rest before-args)
      (apply func (append before-args args))))


(defun sans (plist &rest keys)
  "Erik Naggum's sans function (2.12.2002)."
  (let ((sans ()))
    (loop
      (let ((tail (nth-value 2 (get-properties plist keys))))
        ;; this is how it ends
        (unless tail
          (return (nreconc sans plist)))
        ;; copy all the unmatched keys
        (loop until (eq plist tail) do
              (push (pop plist) sans)
              (push (pop plist) sans))
        ;; skip the matched key
        (setq plist (cddr plist))))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Numbers
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(symbol-macrolet ((nn 13))
  (let ((fac-vec (make-array nn :element-type 'fixnum :initial-element 1)))
    (loop for i from 1 below nn do
	  (setf (aref fac-vec i) (* i (aref fac-vec (1- i)))))
    (defun factorial (n)
      "Returns the factorial of a number or a list of numbers (product of
factorials)."
      (declare (type (or fixnum list) n))
      (flet ((fac (n)
	       (declare (type fixnum n))
	       (if (< n nn)
		   (aref fac-vec n)
		   (loop for i from nn upto n
			 and f = (aref fac-vec (1- nn)) then (* f i)
		       finally (return f)))))
	(if (numberp n)
	    (fac n)
	    (loop with f = 1
		  for k of-type fixnum in n
		  do (setq f (* f (fac k)))
		  finally (return f)))))))

(definline square (x)
  "Computes the square of its argument."
  (* x x))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Conversion
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun make-vector (dim func)
  (let ((vec (make-array dim)))
    (dotimes (i dim vec)
      (setf (aref vec i) (funcall func i)))))

(definline vector-map (func vec)
  (map (type-of vec) func vec))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Special typed vectors
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(deftype int () '(signed-byte 32))
(deftype uint () '(unsigned-byte 32))
(deftype int-vec () '(simple-array int (*)))
(deftype uint-vec () '(simple-array uint (*)))

(definline make-uint-vec (dim &optional (init 0))
  (make-array dim :element-type 'uint :initial-element init))

(deftype positive-fixnum () '(and fixnum unsigned-byte))

(deftype fixnum-vec () '(simple-array fixnum (*)))

(definline fixnum-vec (&rest elements)
  (coerce elements 'fixnum-vec))

(definline make-fixnum-vec (dim &optional (init 0))
  (make-array dim :element-type 'fixnum :initial-element init))

(definline vector-slice (vec offset size)
  "Provides a convenient shorthand for constructing a displaced
double-float array."
  (make-array size
	      :element-type 'double-float
	      :displaced-to vec
	      :displaced-index-offset offset))

(defun vector-cut (vec comp)
  (let ((new-vec (make-array (1- (length vec)) :element-type (array-element-type vec))))
    (dotimes (i (1- (length vec)) new-vec)
      (setf (aref new-vec i)
	    (aref vec (if (< i comp) i (1+ i)))))))

;;; vector operations
(definline ivec+ (&rest vecs)
  (apply #'map 'fixnum-vec #'+ vecs))
(definline ivec- (&rest vecs)
  (apply #'map 'fixnum-vec #'- vecs))

(defmethod for-each ((func function) (lst list))
  (mapc func lst))

(defmethod for-each ((func function) (arr array))
  (loop for entry across arr do
	(funcall func entry)))

(definline vector-last (vec)
  (aref vec (1- (length vec))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; arrays
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun map-array (func arr)
  (let ((result (make-array (array-dimensions arr) :initial-element nil)))
    (dotimes (i (car (array-dimensions arr)))
      (dotimes (j (cadr (array-dimensions arr)))
	(setf (aref result i j)
	      (funcall func (aref arr i j)))))
    result))

(defun array-for-each (func &rest arrays)
  (assert (same-p arrays :test #'equalp :key #'array-dimensions))
  (labels ((loop-index (current-indices further-dims)
	   (if (null further-dims)
	       (apply func (mapcar #'(lambda (array)
				       (apply #'aref array current-indices))
				   arrays))
	       (dotimes (k (car further-dims))
		 (loop-index (cons k current-indices) (cdr further-dims))))))
    (loop-index () (reverse (array-dimensions (car arrays))))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; List operations
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun same-p (seq &key (test #'eql) (key #'identity))
  "Returns t if list consists of equal elements."
  (or (= 0 (length seq))
      (let ((item (funcall key (elt seq 0))))
	(not (find-if-not #'(lambda (elem) (funcall test item (funcall key elem)))
			  seq)))))

(defun encapsulate (p-list dim)
  "Example: (encapsulate 1 3) -> (((1)))"
  (if (zerop dim)
      p-list
      (encapsulate (list p-list) (- dim 1))))

(defun filter (item list &key (key #'identity) (test #'eql))
  "Example: (filter 2 '(1 2 3 4))"
  (remove-if-not  #'(lambda (x) (funcall test x item)) list :key key))

(defun flatten (x)
  "Flattens a tree.
Example: (flatten '((1 2) (3) ((4))))  ->  (1 2 3 4)"
  (labels ((rec (x acc)
	     (cond ((null x) acc)
		   ((not (consp x)) (cons x acc))
		   (t (rec (car x) (rec (cdr x) acc))))))
    (rec x ())))

(defun flatten-1 (lst) (apply #'append lst))

(definline single? (lst) (and (consp lst) (null (cdr lst))))

(definline range (k l)
  (loop for i from k upto l
	collect i))

(definline range< (k l)
  (loop for i from k below l
	collect i))

(defun make-set (start size)
  (loop repeat size
	for j from start
	collect j))

(defun take (k lst)
  (loop for x in lst repeat k
	collect x))

(definline thrice (x) (list x x x))
(definline twice (x) (list x x))

(defun splice (args nrins)
  "Example: (splice '(1 2 3 4) '(1 3))  ->  ((1) (2 3 4))"
  (loop for nrin in nrins
	and tail = args then (nthcdr nrin tail)
	collect (take nrin tail)))

(definline mappend (func &rest lists)
  (apply #'append (apply #'mapcar func lists)))

(defun map-product (func list &rest rest-lists)
  "Applies a function to a product of lists.
Example: (map-product #'cons '(2 3) '(1 4)) -> ((2 . 1) (2 . 4) (3 . 1) (3 . 4))"
  (if (null rest-lists)
      (mapcar func list)
      (mappend #'(lambda (x) (apply #'map-product (curry func x) rest-lists))
	       list)))

(defun mklist (obj)
  (if (listp obj) obj (list obj)))

(definline sum-over (lst) (loop for k in lst summing k))

(defun tree-uniform-number-of-branches (tree)
  (if (listp tree)
      (cons (length tree)
	    (tree-uniform-number-of-branches (car tree)))))

(defun tree-uniformp (trees)
  (labels ((rec (trees dims)
	     (if (listp trees) 
		 (and (listp dims) (= (car dims) (length trees))
		      (every #'(lambda (tree) (rec tree (cdr dims))) trees))
		 (null dims))))
    (rec trees (tree-uniform-number-of-branches trees))))

(defun on-leaves (func tree)
  "Executes func on the leaves of a tree."
  (labels ((rec (node)
	     (if (listp node)
		 (mapc #'rec node)
		 (funcall func node))))
    (rec tree)))

(defun map-tree (func tree)
  "Maps a tree.  (map-tree #'1+ '((1 (2)))) => ((2 (3)))"
  (mapcar
   #'(lambda (item)
       (if (listp item)
	   (map-tree func item)
	   (funcall func item)))
   tree))

(defun rfind-if (func tree)
  "From Graham's book"
  (if (atom tree)
      (and (funcall func tree) tree)
      (or (rfind-if func (car tree))
          (if (cdr tree) (rfind-if func (cdr tree))))))

(defun rfind (item tree &key (key #'identity) (test #'eql))
  (rfind-if #'(lambda (item2) (funcall test item (funcall key item2)))
	    tree))

(defun rmember-if (func tree)
  "From Graham's book"
  (if (atom tree)
      (funcall func tree)
      (or (rmember-if func (car tree))
          (if (cdr tree) (rmember-if func (cdr tree))))))

(defun rmember (item tree &key (key #'identity) (test #'eql))
  (rmember-if #'(lambda (item2) (funcall test item (funcall key item2)))
	    tree))

(defun check-properties (place properties)
  "Checks if all of the properties are in place."
  (every (curry #'getf place) properties))

(defun sort-lexicographically (list &key (key #'identity))
  "Sorts a list of coordinate vectors lexicographically."
  (flet ((my-< (pos1 pos2)
	   (loop for x across pos1 and y across pos2
		 when (< x y) do (return t)
		 finally (return nil))))
    (sort list #'my-< :key key)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Queues (from Graham's book)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun make-queue ()
  "Makes a queue represented by a pair of pointers to start and end."
  (cons nil nil))

(defun enqueue (obj q)
  "Puts an object into the queue."
  (if (null (car q))
      (setf (cdr q) (setf (car q) (list obj)))
    (setf (cdr (cdr q)) (list obj)
	  (cdr q) (cdr (cdr q)))))

(defun dequeue (q)
  "Pops an object from the queue."
  (pop (car q)))

(defun queue->list (q) (car q))

(defun peek-first (q) (caar q))
(defun peek-last (q) (cadr q))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Doubly linked lists
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defstruct (dll-item (:conc-name dli-))
  (object nil)
  (succ nil)
  (pred nil))

(defstruct dll
  (first nil)
  (last nil))

(defun dll-front-insert (obj dll)
  (let ((new (make-dll-item :object obj))
	(first (dll-first dll)))
    (when first
      (setf (dli-pred first) new)
      (setf (dli-succ new) first))
    (unless first
      (setf (dll-last dll) new))
    (setf (dll-first dll) new)))

(defun dll-rear-insert (obj dll)
  (let ((new (make-dll-item :object obj))
	(last (dll-last dll)))
    (when last
      (setf (dli-succ last) new)
      (setf (dli-pred new) last))
    (unless last (setf (dll-first dll) new))
    (setf (dll-last dll) new)))

(defun dll-check-consistency (item dll)
  (let ((predecessor (dli-pred item))
	(successor (dli-succ item)))
    (if (null predecessor)
	 (assert (eq item (dll-first dll)))
	 (assert (eq item (dli-succ predecessor))))
    (if (null successor)
	 (assert (eq item (dll-last dll)))
	 (assert (eq item (dli-pred successor))))))

(defun dll-remove-item (item dll)
  (dll-check-consistency item dll)
  (let ((predecessor (dli-pred item))
	(successor (dli-succ item)))
    (if (null predecessor)
	(setf (dll-first dll) successor)
	(setf (dli-succ predecessor) successor))
    (if (null successor)
	(setf (dll-last dll) predecessor)
	(setf (dli-pred successor) predecessor))))

(defun dll-find (key dll)
  (loop for item = (dll-first dll) then (dli-succ item) until (null item)
	when (eq key (dli-object item)) do (return item)))

(defun dll-remove (obj dll)
  (dll-remove-item (dll-find obj dll) dll))

(defun dll-pop-first (dll)
  (let ((entry (dll-first dll)))
    (when entry
      (dll-remove-item entry dll)
      (dli-object entry))))

(defun dll-pop-last (dll)
  (let ((entry (dll-last dll)))
    (when entry
      (dll-remove-item entry dll)
      (dli-object entry))))

(defun dll->list (dll)
  (loop for item = (dll-first dll) then (dli-succ item) until (null item)
	collect item))

(defun dll-empty-p (dll)
  (null (dll-first dll)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Additional hash-table operations
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmacro dohash ((looping-var hash-table) &body body)
  "Loops through a hash-table.  If looping-var is an atom it loops through
the keys, if it is a list of one element it loops through the values, if it
is a list of two items it loops through key and value."
  (flet ((single? (lst) (and (consp lst) (null (cdr lst)))))
    (let ((key (if (single? looping-var)
		   (gensym)
		   (if (atom looping-var) looping-var (car looping-var))))
	  (value (if (atom looping-var)
		     (gensym)
		     (car (last looping-var)))))
      `(maphash
	#'(lambda (,key ,value)
	    ,@(when (atom looping-var) `((declare (ignore ,value))))
	    ,@(when (single? looping-var) `((declare (ignore ,key))))
	    ,@body)
	,hash-table))))

(defun map-hash-table (func ht)
  "Calls func given in the first argument on each key of the
hash-table from the second argument.  func must be a two-valued
function returning the new key and the new value.  Those pairs are
stored in a new hash-table."
  (let ((new-ht (make-hash-table :test (hash-table-test ht))))
    (dohash ((key val) ht)
      (multiple-value-bind (new-key new-val)
	  (funcall func key val)
	(setf (gethash new-key new-ht) new-val)))
    new-ht))

(defun copy-hash-table (hash-table)
  "Copies a hash-table."
  (let ((new-hash-table (make-hash-table :test (hash-table-test hash-table))))
    (dohash ((key val) hash-table)
      (setf (gethash key new-hash-table) val))
    new-hash-table))

(defun get-arbitrary-key-from-hash-table (ht)
  "Get an arbitrary key from hash-table."
  (dohash (key ht)
    (return-from get-arbitrary-key-from-hash-table key)))

(defun display-ht (hash-table)
  "Displays hash-table in the form key1 --> value1 ..."
  (dohash ((key val) hash-table)
    (format t "~A --> ~A~%" key val)))

(defun hash-table-keys (hash-table)
  "Collects the keys of hash-table into a list."
  (loop for key being the hash-keys of hash-table
	collect key))

(defun hash-table-values (hash-table)
  "Collects the values of hash-table into a list."
  (loop for val being the hash-values of hash-table
	collect val))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Memoizing
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun memoize-1 (func)
  "Primitive memoizer.  Does handle only one-argument and
non-recursive functions."
  (let ((table (make-hash-table)))
    #'(lambda (key)
	(multiple-value-bind (value found)
	    (gethash key table)
	  (if found
	      value
	      (setf (gethash key table)
		    (funcall func key)))))))

(defun memoize (func)
  "Primitive memoizer.  Handles multi-argument non-recursive functions."
  (let ((table (make-hash-table :test #'equalp)))
    #'(lambda (&rest args)
	(princ args)
	(multiple-value-bind (value found)
	    (gethash args table)
	  (if found
	      value
	      (setf (gethash args table)
		    (apply func args)))))))

(defun memoize-1-symbol (funsym)
  "Memoizes 1-argument functions named by the given symbol."
  (let ((unmemoized (symbol-function funsym)) 
	(table (make-hash-table)))
    (setf (symbol-function funsym)
	  #'(lambda (key)
	      (multiple-value-bind (value found)
		  (gethash key table)
		(if found
		    value
		    (setf (gethash key table)
			  (funcall unmemoized key))))))))

(defun memoize-symbol (funsym &key (test 'equal))
  "Memoizes multi-argument functions named by the given symbol."
  (let ((unmemoized (symbol-function funsym)) 
	(table (make-hash-table :test test)))
    (setf (get funsym :memoization-table) table)
    (setf (symbol-function funsym)
	  #'(lambda (&rest args)
	      (multiple-value-bind (value found)
		  (gethash args table)
		(if found
		    value
		    (setf (gethash args table)
			  (apply unmemoized args))))))))

(defmacro defmemo-1 (name &rest rest)
  `(progn
    (defun ,name ,@rest)
    (memoize-1-symbol ',name)
    ))

(defmacro defmemo (name &rest rest)
  (let ((unmemoized (intern (format nil "~A-UNMEMOIZED" name))))
    `(progn
      (defun ,name ,@rest)
      (memoize-symbol ',name)
      (defun ,unmemoized ,@rest)
      )))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Association lists
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(declaim (inline geta))
(defun geta (alist key)
  "An analog to GETF for association lists."
  (cdr (assoc key alist)))

(define-setf-expander geta (alist key &environment env)
  "An analog to (SETF GETF) for association lists."
  (multiple-value-bind (temps vals stores store-form access-form)
      (get-setf-expansion alist env)
    (with-gensyms (keytemp result record maybe-new-record)
      ;; Return the setf expansion for geta
      (values
       `(,@temps ,keytemp ,record ,maybe-new-record)
       `(,@vals
         ,key
         (assoc ,keytemp ,access-form)
         (or ,record (cons ,key nil)))
       (list result)
       `(if ,record
         (setf (cdr ,record) ,result)
         (let ((,(first stores) (cons ,maybe-new-record ,access-form)))
           (setf (cdr ,maybe-new-record) ,result)
           ,store-form
           ,result))
       `(cdr ,maybe-new-record)))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Blackboards
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(declaim (inline make-blackboard getbb))

(defun blackboard (&rest items)
  "Makes items into an blackboard.  Copies if necessary to ensure that no
literal list is modified."
  (if (eq (car items) :blackboard)
      items
      (list* :blackboard t (copy-seq items))))

(defun getbb (blackboard key &optional default)
  "Gets an item from the blackboard.  If there is no such item it returns
the default value."
  (getf (cddr blackboard) key default))

(defun (setf getbb) (value blackboard key)
  "Setter for an blackboard."
  (setf (getf (cddr blackboard) key) value))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Sets
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun set-p (lst)
  (or (null lst)
      (and (not (member (car lst) (cdr lst)))
	   (set-p (cdr lst)))))

(defun maximally-connected (connected disconnected &key (test #'eql) (combine #'adjoin))
  "Finds a maximally connected set by taking the union of the elements in
connected with the sets of disconnected-sets.  Returns the maximally
connected sets and the remaining disconnected ones.

Example: (maximally-connected '(1 2) '((3 4) (2 3) (5 6)) :test #'intersection :combine #'union))
-> values '(1 2 3 4) '((5 6))"
  (loop while
	(dolist (elem disconnected)
	  (when (funcall test elem connected)
	    (setq disconnected (remove elem disconnected))
	    (setq connected (funcall combine connected elem))
	    (return t)))
	finally (return (values connected disconnected))))

;;;(maximally-connected '(1 2) '((3 4) (2 3) (5 6)) :test #'intersection :combine #'union)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Ordered sets
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun ordered-set-difference (set1 set2 &key (test #'eql))
  (loop for elem in set1
	unless (member elem set2 :test test)
	collect elem))

(defun ordered-intersection (set1 set2 &key (test #'eql))
  (loop for elem in set1
	when (member elem set2 :test test)
	collect elem))

(defun ordered-union (set1 set2 &key (test #'eql))
  (append set1 (ordered-set-difference set2 set1 :test test)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Subsets
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; Example: (k-subsets '(1 2 3) 2)  -> ((1 2) (1 3) (2 3))
(defun k-subsets (set k)
  (cond ((minusp k) '())
	((zerop k) (list '()))
	((= k 1) (mapcar #'list set))
	(t (loop for elems on set
		 nconc (mapcar #'(lambda (set) (cons (car elems) set))
				 (k-subsets (cdr elems) (- k 1)))))))

;;; Example: (k->l-subsets '(1 2 3) 1 2)
;;; -> ((1) (2) (3) (1 2) (1 3) (2 3))
(defun k->l-subsets (set k l)
  (loop for i from k to l
	nconc (k-subsets set i)))

(defun subsets (set) (k->l-subsets set 0 (length set)))
(defun nonempty-subsets (set) (k->l-subsets set 1 (length set)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Ordered partitions of natural numbers
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun n-partitions-of-k (n k)
  "Returns a list of all ordered partitions of k into n natural numbers.
Example: (n-partitions-of-k 2 3) -> ((0 3) (1 2) (2 1) (3 0))"
  (cond ((zerop n) (if (zerop k) (list ()) ()))
	((= n 1) (list (list k)))
	(t (loop for i upto k
		 nconc (loop for partition in (n-partitions-of-k (1- n) (- k i))
			     collect (cons i partition))))))

(defun positive-n-partitions-of-k (n k)
  "Returns a list of all positive ordered partitions of k into n natural numbers.
Example: (positive-n-partitions-of-k 2 3) -> ((1 2) (2 1))"
  (cond ((or (< n 1) (zerop k)) ())
	((= n 1) (list (list k)))
	(t (loop for i from 1 to k
		 nconc (loop for partition in (positive-n-partitions-of-k (1- n) (- k i))
			     collect (cons i partition))))))

(defun positive-partitions-of-k (k)
  "Returns a list of all positive ordered partitions of k.
Example: (positive-partitions-of-k 3) -> ((1 2) (2 1))"
  (loop for n from 1 upto k
	nconc (positive-n-partitions-of-k n k)))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Permutations
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; We characterize a permutation pi by a vector consisting of indices from
;;; 0,..., N-1.  

(defun permutation-p (perm)
  "Checks if the argument is a possible permutation vector."
  (loop for i below (length perm)
	unless (find i perm) do (return nil)
	finally (return t)))

(defun identity-permutation-p (perm)
  "Checks if the permutation is the identity."
  (dotimes (i (length perm) t)
    (unless (= i (aref perm i))
      (return nil))))

(defun permute-into (perm index new-index)
  "A permutation acts on some vector v consisting of indices by permuting it
according to result[i] = v[perm[i]]."
  (dotimes (k (length index) new-index)
    (setf (aref new-index k) (aref index (aref perm k)))))

(defun permute (perm index)
  "See permute-into."
  (map 'fixnum-vec (curry #'aref index) perm))

(defun permutation-inverse (perm-vec)
  "Returns the inverse of a permutation (which is a vector of numbers
from 0,...,n-1)."
  (loop with inv-vec = (make-fixnum-vec (length perm-vec) 0)
	for i below (length perm-vec) do
	(setf (aref inv-vec (aref perm-vec i)) i)
	finally (return inv-vec)))

(defun permutation-shifted-inverse (perm-vec)
  (loop with inv-vec = (make-fixnum-vec (length perm-vec) 0)
	for i below (length perm-vec) do
	(setf (aref inv-vec (1- (aref perm-vec i))) (1+ i))
	finally (return inv-vec)))

(defun permutation-signum (perm)
  (let ((perm (copy-seq perm))
	(result 1))
    (dotimes (i (length perm) result)
      (loop until (= (aref perm i) i) do
	    (rotatef (aref perm i) (aref perm (aref perm i)))
	    (setq result (- result))))))
(declaim (ftype (function (array) (member -1 1)) permutation-signum))

;;; Testing
(defun test-utilities ()
  (assert (= 1 (permutation-signum #(5 3 4 0 2 1))))
  (assert (tree-uniformp '((1 2) (3 4) (5 6))))
  (on-leaves #'princ '((1 2) 3))
  ;; test doubly-linked lists
  (let ((x (make-dll)))
    (dll-front-insert 1 x)
    (dll-rear-insert 2 x)
    (dll-rear-insert 3 x)
    (dll-remove 1 x)
    (dll-remove 2 x)
    (dll-pop-first x)
    (dll-pop-first x))
  (let ((ht (make-hash-table)))
    (setf (gethash 1 ht) 1)
    (setf (gethash 2 ht) 4)
    (dohash ((key value) ht)
	    (format t "~&~A -> ~A~%" key value))
    (dohash (key ht) (format t "~A~%" key))
    (dohash ((value) ht) (format t "~A~%" value)))
  (let ((alist ()))
    (setf (geta alist :test) 1)
    (setf (geta alist :hello) 2)
    (setf (geta alist :test) 3)
    (geta alist :test))
  )

;; (test-utilities)
(tests::adjoin-femlisp-test 'test-utilities)

