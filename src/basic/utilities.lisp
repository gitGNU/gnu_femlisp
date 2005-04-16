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

(in-package :fl.utilities)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Functions
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun required-argument ()
  "Calling this function results in an error.  Such a call may be used as
default form when an argument should be supplied."
  (error "A required argument was not supplied."))

#+(or)
(defun compose (&rest functions)
  "Returns the composition of @arg{functions}."
  (if (null functions)
      #'identity
      (destructuring-bind (func1 . rest) (reverse functions)
	#'(lambda (&rest args)
	    (reduce #'(lambda (v f) (funcall f v))
		    rest
		    :initial-value (apply func1 args))))))

(definline curry (func &rest args)
  "Supplies @arg{args} to @arg{func} from the left."
  #'(lambda (&rest after-args)
      (apply func (append args after-args))))

(definline rcurry (func &rest args)
  "Supplies @arg{args} to @arg{func} from the right."
  #'(lambda (&rest before-args)
      (apply func (append before-args args))))

(defun sans (plist &rest keys)
  "Removes the items marked by @arg{keys} from the property list
@arg{plist}.  This function was posted at 2.12.2002 to the
@emph{comp.lang.lisp} newsgroup by Erik Naggum."
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
      "Compute the factorial of @arg{n}.  @arg{n} can also be a list of
numbers in which case the product of the factorials of the components is
computed."
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
  "Return the square of @arg{x}."
  (* x x))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Boxes
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun box (object)
  "Boxes an object."
  (list object))

(defun unbox (box)
  "Getter for a boxed object."
  (car box))

(defun (setf unbox) (value box)
  "Setter for a boxed object."
  (setf (car box) value))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Sequences
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defgeneric for-each (func collection)
  (:documentation "Applies @arg{func} to each element of @arg{collection}.")
  (:method ((func function) (seq sequence))
	   "Applies @arg{func} to each element of the sequence @arg{seq}."
	   (map nil func seq)))

(defun partial-sums (seq)
  "Returns a sequence of the same type as @arg{seq} consisting of its
partial sums."
  (let ((sum 0))
    (map (type-of seq) #'(lambda (x) (incf sum x)) seq)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Vectors
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(definline vector-map (func &rest vecs)
  "Map @arg{vec} with @arg{func} to a vector of the same type."
  (apply #'map (type-of (elt vecs 0)) func vecs))

(deftype positive-fixnum ()
  "Positive fixnum tpye."
  '(and fixnum unsigned-byte))

(deftype fixnum-vec ()
  "Vector with elements of type @code{fixnum}."
  '(simple-array fixnum (*)))

(definline fixnum-vec (&rest elements)
  "Returns a @symbol{FIXNUM-VEC} constructed from the parameters."
  (coerce elements 'fixnum-vec))

(definline make-fixnum-vec (dim &optional (init 0))
  "Construct a @symbol{FIXNUM-VEC} of size @arg{dim} initialized by
@arg{init}."
  (make-array dim :element-type 'fixnum :initial-element init))

(defun vector-cut (vec comp)
  (let ((new-vec (make-array (1- (length vec)) :element-type (array-element-type vec))))
    (dotimes (i (1- (length vec)) new-vec)
      (setf (aref new-vec i)
	    (aref vec (if (< i comp) i (1+ i)))))))

(definline vector-last (vec)
  "Returns the last element of @arg{vec}."
  (aref vec (1- (length vec))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Arrays
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun array-for-each (func &rest arrays)
  "Calls @arg{func} on all element tuples of the array arguments."
  (assert (same-p arrays :test #'equalp :key #'array-dimensions))
  (labels ((loop-index (current-indices further-dims)
	   (if (null further-dims)
	       (apply func (mapcar #'(lambda (array)
				       (apply #'aref array current-indices))
				   arrays))
	       (dotimes (k (car further-dims))
		 (loop-index (cons k current-indices) (cdr further-dims))))))
    (loop-index () (reverse (array-dimensions (car arrays))))))

(defun undisplace-array (array)
  "Return the fundamental array and the start and end positions into it of
the displaced array.  (Erik Naggum, c.l.l. 17.1.2004)"
  (loop with start = 0 do
	(multiple-value-bind (to offset)
	    (array-displacement array)
	  (if to
	      (setq array to
		    start (+ start offset))
	      (return (values array start (+ start (length array))))))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; List operations
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun same-p (sequence &key (test #'eql) (key #'identity))
  "Returns t if @arg{sequence} consists of equal elements."
  (or (= 0 (length sequence))
      (let ((item (funcall key (elt sequence 0))))
	(not (find-if-not #'(lambda (elem) (funcall test item (funcall key elem)))
			  sequence)))))

(defun flatten (tree)
  "Flatten a tree.  Example:
@lisp
  (flatten '((1 2) (3) ((4))))  @result{}  (1 2 3 4)
@end lisp"
  (labels ((rec (tree acc)
	     (cond ((null tree) acc)
		   ((not (consp tree)) (cons tree acc))
		   (t (rec (car tree) (rec (cdr tree) acc))))))
    (rec tree ())))

(defun flatten-1 (lst) (apply #'append lst))

(definline single? (lst) (and (consp lst) (null (cdr lst))))

(definline range<= (k l)
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

(defun splice (items lengths)
  "Breaks the list @arg{items} in pieces of lengths determined by
@arg{nrins}.  Example:
@lisp
  (splice '(1 2 3 4) '(1 3)) @result{} ((1) (2 3 4))
@end lisp"
  (loop for l in lengths
	and tail = items then (nthcdr l tail)
	collect (take l tail)))

(definline mappend (func &rest lists)
  "Map @function{func} over @arg{lists} while appending the results."
  (apply #'append (apply #'mapcar func lists)))

(defun map-product (func list &rest rest-lists)
  "Applies @arg{func} to a product of lists.  Example:
@lisp
  (map-product #'cons '(2 3) '(1 4))
  @result{} ((2 . 1) (2 . 4) (3 . 1) (3 . 4))
@end lisp"
  (if (null rest-lists)
      (mapcar func list)
      (mappend #'(lambda (x) (apply #'map-product (curry func x) rest-lists))
	       list)))

(defun mklist (obj)
  "Wraps @arg{obj} in a list, if it is not already a list."
  (if (listp obj) obj (list obj)))

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
  "Executes @arg{func} on the leaves of @arg{tree}."
  (labels ((rec (node)
	     (if (listp node)
		 (mapc #'rec node)
		 (funcall func node))))
    (rec tree)))

(defun map-tree (func tree)
  "Maps @arg{tree} using @arg{func}.  Example:
@lisp
  (map-tree #'1+ '((1 (2)))) @result{} ((2 (3)))
@end lisp"
  (if (atom tree)
      (funcall func tree)
      (mapcar
       #'(lambda (item)
	   (if (listp item)
	       (map-tree func item)
	       (funcall func item)))
       tree)))

(defun check-properties (place properties)
  "Checks if all of the @arg{properties} are in the property list
@arg{place}."
  (every (curry #'getf place) properties))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Queues (from Graham's book)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun make-queue ()
  "Creates and returns an empty queue."
  (cons nil nil))

(defun enqueue (object queue)
  "Puts @arg{object} into the @arg{queue}."
  (if (null (car queue))
      (setf (cdr queue) (setf (car queue) (list object)))
    (setf (cdr (cdr queue)) (list object)
	  (cdr queue) (cdr (cdr queue)))))

(defun dequeue (queue)
  "Pops an object from @arg{queue}."
  (pop (car queue)))

(defun queue->list (queue)
  "Transforms @arg{queue} to a list."
  (car queue))

(defun list->queue (list)
  "Transforms @arg{list} to a queue."
  (let ((copy (copy-seq list)))
    (cons copy (last copy))))

(defun peek-first (queue)
  "Returns the first item in @arg{queue} (without popping it)."
  (caar queue))
(defun peek-last (queue)
  "Returns the last item in @arg{queue}."
  (cadr queue))

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
	collect (dli-object item)))

(defun list->dll (list)
  (let ((dll (make-dll)))
    (dolist (item list dll)
      (dll-rear-insert item dll))))

(defun dll-empty-p (dll)
  (null (dll-first dll)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Additional hash-table operations
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmacro dohash ((looping-var hash-table) &body body)
  "Loops through @arg{hash-table}.  If @arg{looping-var} is an atom
@emph{key}, loop through the keys; if it is a list of the form
@emph{(value)} loop through the values; if it is a list of the form
@emph{(key value)} loop through key and value."
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

(defun map-hash-table (func hash-table)
  "Call @arg{func} given in the first argument on each key of
@arg{hash-table}.  @arg{func} must return the new key and the new value as
two values.  Those pairs are stored in a new hash-table."
  (let ((new-ht (make-hash-table :test (hash-table-test hash-table))))
    (dohash ((key val) hash-table)
      (multiple-value-bind (new-key new-val)
	  (funcall func key val)
	(setf (gethash new-key new-ht) new-val)))
    new-ht))

(defun copy-hash-table (hash-table)
  "Copy @arg{hash-table}."
  (map-hash-table #'(lambda (key value) (values key value)) hash-table))

(defun get-arbitrary-key-from-hash-table (hash-table)
  "Get an arbitrary key from @arg{hash-table}."
  (dohash (key hash-table)
    (return-from get-arbitrary-key-from-hash-table key)))

(defun display-ht (hash-table)
  "Display @arg{hash-table} in the form key1 -> value1 ..."
  (dohash ((key val) hash-table)
    (format t "~A -> ~A~%" key val)))

(defun hash-table-keys (hash-table)
  "Collect the keys of @arg{hash-table} into a list."
  (loop for key being the hash-keys of hash-table
	collect key))

(defun hash-table-values (hash-table)
  "Collect the values of @arg{hash-table} into a list."
  (loop for val being the hash-values of hash-table
	collect val))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Memoizing
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun memoize-1 (func &key (test 'eql))
  "Memoizes the function @arg{func} which should be a non-recursive
function of one argument."
  (let ((table (make-hash-table :test test)))
    #'(lambda (key)
	(multiple-value-bind (value found)
	    (gethash key table)
	  (if found
	      value
	      (setf (gethash key table)
		    (funcall func key)))))))

(defun memoize-symbol (funsym &key (test 'equal))
  "Memoizes multi-argument functions named by the symbol @arg{funsym}."
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

(defmacro memoize ((defun funsym &rest rest))
  "Defines a function and memoizes it."
  (assert (eq defun 'defun))
  `(progn (defun ,funsym ,@rest)
	  (memoize-symbol ',funsym)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Association lists
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun geta (alist key)
  "An analog to @code{GETF} for association lists."
  (cdr (assoc key alist)))

(define-setf-expander geta (alist key &environment env)
  "An analog to @code{(SETF GETF)} for association lists."
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

(defclass blackboard ()
  ((items :initform () :initarg :items :documentation
	  "A property list of items on the blackboard."))
  (:documentation
   "A blackboard where data items can be put and extracted using the
function @code{GETBB}."))

(defun blackboard (&rest items)
  "Make the property list supplied in @arg{items} into an blackboard.
Copies @arg{items} to make sure that no literal list is modified."
  (make-instance 'blackboard :items (copy-seq items)))

(defmethod describe-object ((blackboard blackboard) stream)
  (format stream "~&~A is a blackboard containing the following items:~%~{~S -> ~S~%~}"
	  blackboard (slot-value blackboard 'items)))

(defun getbb (blackboard key &optional default)
  "Get the item for @arg{key} from @arg{blackboard}.  If there is no such
@arg{key} return @arg{default}."
  (getf (slot-value blackboard 'items) key default))

(defun (setf getbb) (value blackboard key)
  "Setter corresponding to @code{GETBB}."
  (setf (getf (slot-value blackboard 'items) key) value))

(eval-when (:compile-toplevel :load-toplevel :execute)
  (defun with-items-expander (prop-vars prop-names prop-defaults plist body)
    "Helper function for @code{with-items}."
    `(progn
      #+(or)
      ,@(mapcar #'(lambda (prop default)
		    `(setf (getf ,plist ',prop) (getf ,plist ',prop ,default)))
		prop-names prop-defaults)
      (symbol-macrolet
	    ,(mapcar #'(lambda (var prop default)
			 `(,var (getf ,plist ',prop ,default)))
		     prop-vars prop-names prop-defaults)
	  ,@body))))

(defmacro with-items (properties blackboard-form &body body)
  "Work with the items on @arg{blackboard} corresponding to
@arg{properties}.  If some property is a list, the second element is the
default value and the third is an alias to be used to refer to this
parameter.  Example:
@lisp
  (with-items (&key sol (rhs nil rhs-high)) blackboard
     (setq sol rhs-high))
@end lisp"
  (with-gensyms (blackboard) 
    (let (key-p prop-vars prop-defaults prop-names)
      (dolist (prop properties)
	(cond
	  ((eq prop '&key) (assert (not key-p))
	   (setq key-p t))
	  (t (let ((prop-var
		    (if (atom prop)
			prop
			(or (nth 2 prop) (nth 0 prop))))
		   (prop-default (and (listp prop) (nth 1 prop)))
		   (prop-name (if (atom prop) prop (nth 0 prop))))
	       ;;(format t "~A ~A ~A" prop-var prop-default prop-name)
	       (push prop-var prop-vars)
	       (push prop-default prop-defaults)
	       (push (if key-p
			 (intern (symbol-name prop-name) :keyword)
			 prop-name)
		     prop-names)))))
      `(let ((,blackboard ,blackboard-form))
	,(with-items-expander prop-vars prop-names prop-defaults
			      `(slot-value ,blackboard 'items) body)))))

(defun transfer-bb (from-bb to-bb items &key ensure)
  "Transfer @arg{items} between the blackboards @arg{from-bb} and
@arg{to-bb}.  When @arg{ensure} is set, an existing item is not modified."
  (dolist (item items)
    (destructuring-bind (from to)
	(if (consp item) item (list item item))
      (if ensure
	  (ensure (getbb to-bb to)) (getbb from-bb from))
	  (setf (getbb to-bb to) (getbb from-bb from)))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Sets
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun set-p (list)
  "Checks if @arg{list} is a set, i.e. if no members occur twice."
  (or (null list)
      (and (not (member (car list) (cdr list)))
	   (set-p (cdr list)))))

(defun maximally-connected (connected disconnected &key (test #'eql) (combine #'adjoin))
  "Finds a maximally connected set by taking the union of the elements in
connected with the sets of disconnected-sets.  Returns the maximally
connected sets and the remaining disconnected ones.  Example:

@lisp
  (maximally-connected '(1 2) '((3 4) (2 3) (5 6))
                       :test #'intersection :combine #'union)
  @result{} (1 2 3 4), ((5 6))
@end lisp"
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

(defun k-subsets (set k)
  "Returns all subsets of @arg{set} with length @arg{k}.  Example:
@lisp
  (k-subsets '(1 2 3) 2)
  @result{} ((1 2) (1 3) (2 3))
@end lisp"
  (cond ((minusp k) '())
	((zerop k) (list '()))
	((= k 1) (mapcar #'list set))
	(t (loop for elems on set
		 nconc (mapcar #'(lambda (set) (cons (car elems) set))
				 (k-subsets (cdr elems) (- k 1)))))))

(defun k->l-subsets (set k l)
  "Returns all subsets of @arg{set} with length between @arg{k} and
@arg{l}.  Example:
@lisp
  (k->l-subsets '(1 2 3) 1 2) @result{}
  ((1) (2) (3) (1 2) (1 3) (2 3))
@end lisp"
  (loop for i from k to l
	nconc (k-subsets set i)))

(defun subsets (set)
  "Returns a list of all subsets of @arg{set}."
  (k->l-subsets set 0 (length set)))

(defun nonempty-subsets (set)
  "Returns a list of all nonempty subsets of @arg{set}."
  (k->l-subsets set 1 (length set)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Ordered partitions of natural numbers
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun n-partitions-of-k (n k)
  "Returns a list of all ordered partitions of @arg{k} into @arg{n} natural
numbers.  Example:
@lisp
  (n-partitions-of-k 2 3)
   @result{} ((0 3) (1 2) (2 1) (3 0))
@end lisp"
  (cond ((zerop n) (if (zerop k) (list ()) ()))
	((= n 1) (list (list k)))
	(t (loop for i upto k
		 nconc (loop for partition in (n-partitions-of-k (1- n) (- k i))
			     collect (cons i partition))))))

(defun positive-n-partitions-of-k (n k)
  "Returns a list of all positive ordered partitions of @arg{k} into
@arg{n} natural numbers.  Example:
@lisp
  (positive-n-partitions-of-k 2 3)
  @result{} ((1 2) (2 1))
@end lisp"
  (cond ((or (< n 1) (zerop k)) ())
	((= n 1) (list (list k)))
	(t (loop for i from 1 to k
		 nconc (loop for partition in (positive-n-partitions-of-k (1- n) (- k i))
			     collect (cons i partition))))))

(defun positive-partitions-of-k (k)
  "Returns a list of all positive ordered partitions of @arg{k}.
Example:
@lisp
  (positive-partitions-of-k 3) @result{} ((1 2) (2 1))
@end lisp"
  (loop for n from 1 upto k
	nconc (positive-n-partitions-of-k n k)))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Permutations
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun permutation-p (perm)
  "Checks if @arg{perm} is a possible permutation vector.  A permutation pi
is characterized by a vector containing the indices from 0,...,
@function{length}(@arg{perm})-1 in some order."
  (loop for i below (length perm)
	unless (find i perm) do (return nil)
	finally (return t)))

(defun identity-permutation-p (perm)
  "Checks if the permutation is the identity."
  (dotimes (i (length perm) t)
    (unless (= i (aref perm i))
      (return nil))))

(defun permute-into (perm v result)
  "A permutation @arg{perm} acts on the vector @arg{v} by permuting it
according to @math{result[i] = v[perm[i]]}."
  (dotimes (k (length v) result)
    (setf (aref result k) (aref v (aref perm k)))))

(defun permute (perm index)
  "Functional version of @function{permute-into}."
  (map 'fixnum-vec (curry #'aref index) perm))

(defun permutation-inverse (perm-vec)
  "Returns the inverse of the permutation given by @arg{perm-vec}."
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
  (range<= 1 5)
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
  (let ((list '(1 2 3)))
    (equalp (dll->list (list->dll list)) list))
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
  (let ((blackboard (blackboard)))
    (with-items (&key test hello) blackboard
      (describe blackboard)
      (setf test 2)
      (describe blackboard)
      (setf test 3)
      (describe blackboard))
      (getbb blackboard :test))
  (let ((a (box 1)))
    (list
     (fluid-let (((unbox a) 3))
       (unbox a))
     (unbox a)))
  (partial-sums #(1 2 3))
  )

;; (fl.utilities::test-utilities)
(fl.tests:adjoin-test 'test-utilities)

