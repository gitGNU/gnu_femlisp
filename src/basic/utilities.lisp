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

 (defgeneric evaluate (f x)
  (:documentation "Generic evaluation of functions on an argument.  Numbers and
arrays are treated as constants.  Special evaluation is defined for multivariate
polynomials on vectors and for <function> objects.")
  (:method (object x)
    "The default method treats object as a constant function.  This is a
dubious feature on which one should probably not rely."
    (declare (ignore x))
    object)
  (:method ((f function) x)
    "Lisp functions are evaluated by @function{funcall}."
    (funcall f x)))

(defgeneric compose-2 (f g)
  (:documentation "Composes two function objects @arg{f} and @arg{g}.")
  (:method (f g)
    #'(lambda (x) (evaluate f (evaluate g x)))))

(defun compose (&rest functions)
  "Returns the composition of @arg{functions}."
  (cond ((null functions) #'identity)
	((single? functions) (car functions))
	(t (compose-2 (car functions) (apply #'compose (cdr functions))))))

(inlining
 (defun curry (func &rest args)
   "Supplies @arg{args} to @arg{func} from the left."
   #'(lambda (&rest after-args)
       (apply func (append args after-args)))))

(inlining
 (defun rcurry (func &rest args)
   "Supplies @arg{args} to @arg{func} from the right."
   #'(lambda (&rest before-args)
       (apply func (append before-args args)))))

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
;;;; Booleans
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun xor (a b)
  (not (eql (not a) (not b))))

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
  (vector object))

(defun unbox (box)
  "Getter for a boxed object."
  (aref box 0))

(defun (setf unbox) (value box)
  "Setter for a boxed object."
  (setf (aref box 0) value))

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

(defmacro mapf (var func)
  "Replaces the value of the variable @arg{var} by setting it to its map
with @arg{func}."
  `(setf ,var (map (type-of ,var) ,func ,var)))

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

(inlining
 (defun vector-last (vec)
   "Returns the last element of @arg{vec}."
   (aref vec (1- (length vec)))))

(inlining
 (defun constant-vector (dim value)
   "Returns a uniform constant vector of which all elements are @arg{value}."
   (make-array dim :element-type (type-of value) :initial-element value)))

(inlining  ; useful for propagating type information
 (defun zero-vector (dim element-type)
  "Returns a uniform vector for the given element type filled with zeros."
  (make-array dim :element-type element-type
	      :initial-element (coerce 0 element-type))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Arrays
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmacro with-array ((access array) &body body)
  "Posted by Kent Pitman at cll,8.5.2007 as 'array-access'."
  (let ((temp (gensym (symbol-name access))))
   `(let ((,temp ,array))
     (flet ((,access (&rest indexes)
              (apply #'aref ,temp indexes))
           ((setf ,access) (val &rest indexes)
              (setf (apply #'aref ,temp indexes) val)))
      (declare (inline #',access))
      ,@body))))

(defun for-each-tuple (func limits)
  "Calls @arg{func} on each tuple greater or equal to (0 ... 0) and below
@arg{dims}."
  (labels ((helper (tuple limits)
	     (if (null limits)
		 (funcall func tuple)
		 (dotimes (i (car limits))
		   (helper (cons i tuple) (cdr limits))))))
    (helper () (reverse limits))))

(defmacro dotuple ((index limits &optional result) &body body)
  "Loops through each tuple below @arg{limits}."
  (assert (symbolp index))
  `(progn
    (for-each-tuple (lambda (,index) ,@body) ,limits)
    ,result))

(defun array-for-each (func &rest arrays)
  "Calls @arg{func} on all element tuples of the array arguments."
  (assert (same-p arrays :test #'equalp :key #'array-dimensions))
  (dotuple (indices (array-dimensions (car arrays)))
    (apply func (mapcar #'(lambda (array)
			    (apply #'aref array indices))
			arrays))))

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

(defun modify (lst &key add remove)
  (union add (set-difference lst remove)))

(definline single? (lst) (and (consp lst) (null (cdr lst))))

(definline range<= (k l)
  (loop for i from k upto l
	collect i))

(definline range< (k l)
  (loop for i from k below l
	collect i))

(defun take (k lst)
  (loop for x in lst repeat k
	collect x))

(definline thrice (x) (list x x x))
(definline twice (x) (list x x))

(defun split-by-length (items lengths)
  "Breaks the list @arg{items} in pieces of lengths determined by
@arg{nrins}.  Example:
@lisp
  (split-by-length '(1 2 3 4) '(1 3)) @result{} ((1) (2 3 4))
@end lisp"
  (loop for l in lengths
	and tail = items then (nthcdr l tail)
	collect (take l tail)))

(defun mappend (func &rest lists)
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

(defun check-properties (place properties)
  "Checks if all of the @arg{properties} are in the property list
@arg{place}."
  (every (curry #'getf place) properties))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Tree operations
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun on-leaves (func tree)
  "Executes @arg{func} on the leaves of @arg{tree}."
  (labels ((rec (node)
	     (if (listp node)
		 (mapc #'rec node)
		 (funcall func node))))
    (rec tree)))

(defun find-leaf-if (test tree &key (key #'identity))
  "Finds a leaf in @arg{tree} where @arg{test} returns true."
  (on-leaves (lambda (leaf)
	       (when (funcall test (funcall key leaf))
		 (return-from find-leaf-if leaf)))
	     tree)
  nil)

(defun find-leaf (atom tree &key (key #'identity) (test #'eql))
  "Finds atom in @arg{tree}."
  (find-leaf-if (curry test atom) tree :key key))

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

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Queues
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass queue ()
  ((head :accessor head :initform nil)
   (tail :accessor tail :initform nil))
  (:documentation "Queue class accessed with @function{enqueue} and
@function{dequeue}."))

(defgeneric enqueue (object queue)
  (:documentation "Puts @arg{object} into the @arg{queue}.")
  (:method (object (queue queue))
    (with-slots (head tail) queue
      (if (null head)
	  (setf tail (setf head (list object)))
	  (setf (cdr tail) (list object)
		tail (cdr tail)))
      nil)))

(defgeneric emptyp (queue)
  (:documentation "Tests if @arg{queue} is empty.")
  (:method ((queue queue))
    (not (head queue))))

(defgeneric dequeue (queue)
  (:documentation "Pops an object from @arg{queue}.  Returns as second
value T if the queue was empty.")
  (:method ((queue queue))
    (let ((emptyp (emptyp queue)))
      (values (pop (head queue)) emptyp))))

(defgeneric finish (queue)
  (:documentation "Finishes the queue.  Nothing can be written to it
  afterwards.  This function is mostly useful when different threads write
  and read from the queue."))

(defun queue->list (queue)
  "Transforms @arg{queue} to a list."
  (head queue))

(defun list->queue (list)
  "Transforms @arg{list} to a queue."
  (let ((queue (make-instance 'queue)))
    (setf (head queue) (copy-seq list))
    (setf (tail queue) (last (head queue)))
    queue))

(defgeneric dequeue-all (queue)
  (:documentation "Clears @arg{queue} and returns content as a list.")
  (:method ((queue queue))
    (prog1 (head queue) 
      (setf (head queue) nil (tail queue) nil))))

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

(defun dll-front-insert (obj dll &optional insert-item-p)
  "Inserts @arg{obj} in @arg{dll}.  It returns the newly created
@class{dll-item}."
  (let ((new (if (or insert-item-p (not (typep obj 'dll-item)))
		 (make-dll-item :object obj)
		 obj))
	(first (dll-first dll)))
    (when first
      (setf (dli-pred first) new)
      (setf (dli-succ new) first))
    (unless first
      (setf (dll-last dll) new))
    (setf (dll-first dll) new)))

(defun dll-rear-insert (obj dll &optional insert-item-p)
  (let ((new (if (or insert-item-p (not (typep obj 'dll-item)))
		 (make-dll-item :object obj)
		 obj))
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
      `(block nil
	(maphash
	 #'(lambda (,key ,value)
	     ,@(when (atom looping-var) `((declare (ignore ,value))))
	     ,@(when (single? looping-var) `((declare (ignore ,key))))
	     ,@body)
	 ,hash-table)))))

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

(defun collect-in-hash-table (list identifier &optional (type 'eql))
  "Puts the items from @arg{list} in a hash-table identified by
@arg{identifier}."
  (lret ((table (make-hash-table :test type)))
    (loop for item in list do
	  (setf (gethash (funcall identifier item) table) item))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; iteration (loop+)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass range ()
  ((from :initform 0 :initarg :from :documentation
	 "Start of range, defaults to 0.")
   (to :initform nil :initarg :to :documentation
       "Inclusive end of range, defaults to infinity.")
   (below :initform nil :initarg :below :documentation
       "Exclusive end of range, defaults to infinity.")
   (by :initform 1 :initarg :by :documentation
       "Step size."))
  (:documentation "Range of numbers for iteration."))

(defun range (&rest args &key to below)
  "Constructor for a range of numbers."
  (when (and to below)
    (error "Only one limit should be given."))
  (apply #'make-instance 'range args))

(defgeneric iterator (x)
  (:documentation "Returns an iterator for @arg{x}.")
  (:method ((x range)) (slot-value x 'from))
  (:method ((x vector)) 0)
  (:method ((x list)) x))

(defgeneric iterator-next (vec iterator)
  (:documentation "Returns an incremented @arg{iterator}.")
  (:method ((range range) i) (+ i (slot-value range 'by)))
  (:method ((vec vector) i) (1+ i))
  (:method ((vec list) tail) (cdr tail)))

(defgeneric iterator-end-p (vec iterator)
  (:method ((range range) i)
    (with-slots (to below) range
      (cond (to (> i to))
	    (below (>= i below)))))
  (:method ((vec vector) i) (>= i (length vec)))
  (:method ((vec list) tail) (null tail)))
  
(defgeneric reference (vec iterator)
  (:documentation "Reader for the element of @arg{vec} referenced by @arg{iterator}.")
  (:method ((range range) i) i)
  (:method ((vec vector) i) (aref vec i))
  (:method ((vec list) tail) (car tail)))
  
(defgeneric (setf reference) (value vec iterator)
  (:documentation "Setter for the element of @arg{vec} referenced by @arg{iterator}.")
  (:method (value (vec vector) i)
    (setf (aref vec i) value))
  (:method (value (vec list) tail)
    (setf (car tail) value)))

(defmacro loop+ (items &body body)
  "Iterates @arg{body} over @arg{items}.  Example:
@lisp
  (let ((x (make-array 10))
	(y (make-list 10 :initial-element 1)))
    (loop+ ((xc x) (yc y) i) doing
       (setf xc (+ i yc))
       finally (return x)))
@end lisp"
  (let ((vectors (loop for i below (length items)
		      collect (gensym (format nil "V~D" i))))
	(iterators (loop for i below (length items)
		      collect (gensym (format nil "I~D" i)))))
    `(let ,(mapcar #'(lambda (vector item)
		       (list vector (if (listp item) (cadr item) `(range))))
		   vectors items)
       (symbol-macrolet ,(mapcar #'(lambda (vector iterator item)
				     `(,(if (listp item) (car item) item)
					(reference ,vector ,iterator)))
				 vectors iterators items)
	   (loop ,@(loop for iterator in iterators and vector in vectors appending
			 `(for ,iterator = (iterator ,vector)
			   then (iterator-next ,vector ,iterator)))
		 ,@(loop for iterator in iterators and vector in vectors appending
			 `(until (iterator-end-p ,vector ,iterator)))
	    ,@body
	 )))))

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

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Association lists
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun geta (alist key)
  "An analog to @code{GETF} for association lists."
  (cdr (assoc key alist)))

#-gcl
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
  "Make the property list supplied in @arg{items} into a blackboard.
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
  "Returns the sign of @arg{perm}."
  (let ((perm (copy-seq perm))
	(result 1))
    (dotimes (i (length perm) result)
      (loop until (= (aref perm i) i) do
	    (rotatef (aref perm i) (aref perm (aref perm i)))
	    (setq result (- result))))))
(declaim (ftype (function (array) (member -1 1)) permutation-signum))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; safe sorting
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun safe-sort (seq &rest args)
  (apply #'stable-sort (copy-seq seq) args))

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
    (dohash ((value) ht) (format t "~A~%" value)
	    (return)))
  (let ((alist ()))
    (setf (geta alist :test) 1)
    (setf (geta alist :hello) 2)
    (setf (geta alist :test) 3)
    (geta alist :test))
  (let ((blackboard (blackboard)))
    (with-items (&key test) blackboard
      (describe blackboard)
      (setf test 2)
      (describe blackboard)
      (setf test 3)
      (describe blackboard))
      (getbb blackboard :test))
  (maximally-connected '(1 2) '((3 4) (2 3) (5 6))
                       :test #'intersection :combine #'union)
  (let ((a (box 1)))
    (list
     (fluid-let (((unbox a) 3))
       (unbox a))
     (unbox a)))
  (partial-sums #(1 2 3))
  (let ((x (list->queue '(1 2 3))))
    (dequeue x))
  )

;; (fl.utilities::test-utilities)
#+(or)(fl.tests:adjoin-test 'test-utilities)

