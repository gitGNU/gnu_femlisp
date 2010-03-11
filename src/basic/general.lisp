;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; general.lisp - globally useful definitions
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;
;;; Copyright (C) 2003-2006 Nicolas Neuss, University of Heidelberg.
;;; Copyright (C) 2006- Nicolas Neuss, University of Karlsruhe.
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
;;; NO EVENT SHALL THE AUTHOR, THE UNIVERSITIES HEIDELBERG AND KARLSRUHE OR
;;; OTHER CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
;;; SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
;;; LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
;;; DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
;;; THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
;;; (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
;;; OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(in-package :fl.utilities)

;;;; This files provides definitions which have influence across several
;;;; other packages.

(defgeneric make-analog (obj)
  (:documentation "Generate an analogous but empty data structure.")
  (:method (obj) (make-instance (class-of obj))))

(defun file-documentation (docstring)
  "If the manual is sorted by file, the string handed to this function
describes the use of the respective file."
  (declare (ignore docstring)))

(defun concept-documentation (docstring)
  "Documents a certain concept."
  (declare (ignore docstring)))

(defclass property-mixin ()
  ((properties :accessor properties :initform () :initarg :properties
	       :type list :documentation
    "A property list which is used to store unstructured information about
this object."))
  (:documentation "A mixin which adds a slot of properties to the class."))

;;; old interface: deprecated

#+(or)
(defun property-set-p (object property)
  "Returns T if @arg{property} is found in the object's properties."
  (get-properties (slot-value object 'properties) (list property)))

#+(or)
(defmacro with-properties (properties object &body body)
  "Work with @arg{properties} on the property list of @arg{object}."
  (with-gensyms (obj)
    `(let ((,obj ,object))
       (symbol-macrolet ,(loop for prop in properties collect
			      `(,prop (getf (properties ,obj) ',prop)))
	 ,@body))))

;;; old/new interface

(defun get-property (object property)
  "Gets @arg{property} for @arg{object}.  Returns NIL also if
@arg{property} is not available."
  (getf (slot-value object 'properties) property))

(defun (setf get-property) (value object property)
  "Sets the property @arg{property} of @arg{problem} to @arg{value}."
  (setf (getf (slot-value object 'properties) property) value))

;;; new interface using CLOS
(defmethod shared-initialize :after ((object property-mixin) slot-names
				     &rest initargs &key &allow-other-keys)
  (declare (ignore slot-names))
  (let* ((class (class-of object))
	 (ordinary-slot-initargs (mappend #'fl.amop:slot-definition-initargs
					  (fl.amop:compute-slots class))))
    (loop for (key value) on initargs by #'cddr
	 unless (member key ordinary-slot-initargs) do
	 (setf (slot-value object (intern (symbol-name key) *package*))
	       value))))

(defmethod slot-missing (class (object property-mixin) slot-name
			 (operation (eql 'slot-value)) &optional new-value)
  (declare (ignore class new-value))
  (getf (slot-value object 'properties) slot-name))

(defmethod slot-missing (class (object property-mixin) slot-name
			 (operation (eql 'setf)) &optional new-value)
  (declare (ignore class))
  (setf (getf (slot-value object 'properties) slot-name) new-value))

(defmethod slot-missing (class (object property-mixin) slot-name
			 (operation (eql 'slot-boundp)) &optional new-value)
  (declare (ignore class slot-name new-value))
  t)

(defmethod slot-missing (class (object property-mixin) slot-name
			 (operation (eql 'slot-makunbound)) &optional new-value)
  (declare (ignore class new-value))
  (error "Property slot '~A' cannot be unbound" slot-name))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; generic dictionary access
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defgeneric dic-ref (dictionary key)
  (:documentation "Returns the value of the entry for @arg{key} in
@arg{dictionary}.")
  (:method ((table list) key)
    (let ((entry (assoc key table)))
      (values (cdr entry) entry)))
  (:method ((table hash-table) key)
    (gethash key table))
  (:method ((vec vector) index)
    (aref vec index))
  (:method ((func function) key)
    (funcall func key)))

(defgeneric (setf dic-ref) (value dictionary key)
  (:documentation "Sets the entry @arg{key} in @arg{dictionary} to
@arg{value}.")
  (:method (value (table hash-table) key)
    (setf (gethash key table) value))
  (:method (value (vec vector) index)
    (setf (aref vec index) value)))

(defgeneric dic-for-each (function dictionary)
  (:documentation "Applies @arg{function} to each key-value pair in
  @arg{dictionary}.")
  (:method (func (dic list))
    (loop for (key . value) in dic
       do (funcall func key value)))
  (:method (func (dic hash-table))
    (maphash func dic))
  (:method (func (vec vector))
    (loop for i from 0 and x across vec do
      (funcall func i x))))

(defgeneric dic-for-each-key (function dictionary)
  (:documentation "Applies @arg{function} to each key in @arg{dictionary}.")
  (:method (func dic)
    (dic-for-each (lambda (key value)
                    (declare (ignore value))
                    (funcall func key))
                  dic)))

(defgeneric dic-for-each-value (function dictionary)
  (:documentation "Applies @arg{function} to each value in @arg{dictionary}.")
  (:method (func dic)
    (dic-for-each (lambda (key value)
                    (declare (ignore key))
                    (funcall func value))
                  dic)))

(defgeneric dic-remove (key dictionary)
  (:documentation "Removes @arg{key} from @arg{dictionary}.")
  (:method (key (dic hash-table))
    (remhash key dic)))

(defgeneric dic-empty-p (dictionary)
  (:documentation "Tests if @arg{dictionary} is empty.")
  (:method (dic)
    (not (mapper-some (constantly t) dic)))
  (:method ((dic hash-table))
    (zerop (hash-table-count dic))))

;;; derived
(defgeneric keys (dic)
  (:documentation "Returns a list of all keys of @arg{dic}.")
  (:method (dic)
    (mapper-collect #'dic-for-each-key dic)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Dictionaries and memoization
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; small-cache-dictionary
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass small-cache-dictionary ()
  ((size :reader size :initform 1 :initarg :size)
   (test :initform 'eql :initarg :test)
   (store :reader store :initform ()))
  (:documentation "Implementation with lists."))

(defmethod initialize-instance :after ((dic small-cache-dictionary) &key &allow-other-keys)
  (assert (plusp (size dic))))

(defmethod dic-ref ((dic small-cache-dictionary) key)
  (quickly
    (with-slots (store test) dic
      (loop for tail on store
	 and prev-tail = nil then tail do
	 (when (slowly ; gets rid of an SBCL optimization note
		 (funcall test (caar tail) key))
	   (when prev-tail
	     (setf (cdr prev-tail) (cdr tail)
		   (cdr tail) store
		   store tail))
	   (return (values (cdar store) t)))))))

(defmethod (setf dic-ref) (value (dic small-cache-dictionary) key)
  (with-slots (store size test) dic
    (loop for tail on store
	  and prev-tail = nil then tail
	  and i from 1 do
	  (when (funcall test (caar tail) key)
	    ;; eliminate entry from store
	    (cond (prev-tail
		   (setf (cdr prev-tail) (cdr tail))
		   (return))
		  (t (return-from dic-ref
		       (setf (cdr tail) value)))))
	  (when (null (cdr tail))
	    (when (= i size)
	      ;; full: drop last entry
	      (if prev-tail
		  (setf (cdr prev-tail) nil)
		  (setf store nil)))
	    (return)))
    (pushnew (cons key value) store)
    value))

(defmethod dic-for-each (func (dic small-cache-dictionary))
  (dic-for-each func (store dic)))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; sorted-hash-table
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass sorted-hash-table ()
  ((test :initform #'eql :initarg :test)
   (store :initform (make-dll))
   (item-store))
  (:documentation "This is a hash-table where entries are guaranteed to be
  accessed in the order in which they are interned."))

(defmethod initialize-instance :after
    ((dic sorted-hash-table) &key &allow-other-keys)
  (with-slots (test item-store size) dic
    (setf item-store (make-hash-table :test test))))

(defmethod dic-ref ((dic sorted-hash-table) key)
  (with-slots (store item-store) dic
    (whereas ((item (gethash key item-store)))
      (values (cdr (dli-object item)) t))))

(defmethod (setf dic-ref) (value (dic sorted-hash-table) key)
  (with-slots (store item-store size test) dic
    (let ((item (gethash key item-store)))
      (if item
          (setf (cdr (dli-object item)) value)
          (setf (gethash key item-store)
                (dll-rear-insert (cons key value) store)))))
  value)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; cache-dictionary
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass cache-dictionary ()
  ((size :reader size :initform 1 :initarg :size)
   (test :initform #'eql :initarg :test)
   (store :initform (make-dll))
   (item-store))
  (:documentation "This is a dictionary which can only hold a certain
  number of items.  If more are inserted, only the most recently accessed
  items are kept."))

(defmethod initialize-instance :after ((dic cache-dictionary) &key &allow-other-keys)
  (with-slots (test item-store size) dic
    (assert (or (null size) (plusp size)))
    (setf item-store (make-hash-table :test test :size size))
    ))

(defmethod dic-ref ((dic cache-dictionary) key)
  (with-slots (store item-store) dic
    (whereas ((item (gethash key item-store)))
      ;; when in interior of list, put it at front
      (when (dli-pred item)
	(dll-remove-item item store)
	(dll-front-insert item store))
      (values (cdr (dli-object item)) t))))

(defmethod (setf dic-ref) (value (dic cache-dictionary) key)
  (with-slots (store item-store size test) dic
    (let ((item (gethash key item-store)))
      (cond
        (item
         ;; set value of item
         (setf (cdr (dli-object item)) value)
         ;; ensure that the item is the first
         (when (dli-pred item)
           (dll-remove-item store item)
           (dll-front-insert item store)))
        (t
         (setf (gethash key item-store)
               (dll-front-insert (cons key value) store))
         ;; ensure the size constraint
         (when (and size (> (hash-table-count item-store) size))
           (remhash (car (dll-pop-last store)) item-store))))))
  value)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; computed-value-dictionary
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass computed-value-dictionary ()
  ((keys :reader keys :type sequence :initarg :keys)
   (compute :reader compute :type function :initarg :compute
            :documentation "A (usually memoized) function computing the values")))

(defmethod dic-for-each-key (func (dic computed-value-dictionary))
  (for-each func (keys dic)))

(defmethod dic-for-each (func (dic computed-value-dictionary))
  (for-each (lambda (key)
              (funcall func key (funcall (compute dic) key)))
            (keys dic)))

(defmethod dic-ref ((dic computed-value-dictionary) key)
  (funcall (compute dic) key))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Memoization
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmacro memoizing-let (bindings &body body)
  "The @arg{body} is memoized for the keys given as a list of
@arg{bindings}."
  (declare (ignore bindings body))
  (error "This global macro should be used only within the context of
@macro{with-memoization}."))

(defmacro memoizing (defun-expr)
  "As a global macro it memoizes the following function definition.  Inside a
@function{with-memoization} environment, it memoizes its function argument."
  (assert (eq (first defun-expr) 'defun))
  `(prog1
       ,defun-expr
     (memoize-symbol ',(second defun-expr))))

(defun memoization-table (size test)
  (if size
      (make-instance (if (< size 100) 'small-cache-dictionary 'cache-dictionary)
                     :size size :test test)
      (make-hash-table :test test)))

(defmacro with-memoization ((&key (type :global) size id (debug t) (test ''equal)) &body body)
  "Sets up a memoization environment consisting of a table, and a captured
symbol @symbol{memoizing-let} memoizing its body depending on the
arguments.  Example of usage:
@lisp
  (with-memoization (:type :local :id 'test)
    (defun test (n)
      (memoizing-let ((k (* 2 n)))
	(sleep 1)
	(* k k))))
@end lisp
If @arg{type} is :global, the table is thread-safe and the same for all
threads, if @arg{type} is :local, it is special for each thread."
  (when (eq type :local) (assert id))
  (with-gensyms (mutex table key key-exprs value-body func foundp value)
    `(let ((,mutex (fl.multiprocessing:make-mutex))
	   (,table ,(when (eq type :global)
                          `(memoization-table ,size ,test))))
      (declare (ignorable ,mutex ,table))
      (flet ((,table ()
	       ,(ecase type
		       (:global table)
		       (:local
			`(or (getf fl.multiprocessing:*thread-local-memoization-table* ,id)
			  (setf (getf fl.multiprocessing:*thread-local-memoization-table* ,id)
                                (memoization-table ,size ,test)))))))
	;; the memoizing-let symbol is captured, 
	(macrolet ((memoizing (,func)
		     `(lambda (&rest ,',key)
		       (fl.multiprocessing:with-mutex (,',mutex)
			 (multiple-value-bind (,',value ,',foundp)
			     (dic-ref ,'(,table) ,',key)
			   (if ,',foundp
			       ,',value
			       (progn
				 ,',(when debug `(dbg :memoize "Memoizing: id=~A, key=~A" ,id ,key))
				 (setf (dic-ref ,'(,table) ,',key)
				       (apply ,,func ,',key))))))))
		   (memoizing-let (,key-exprs &body ,value-body)
		     `(let ,,key-exprs
		       ;; declares should be inserted here from value-body
		       (fl.multiprocessing:with-mutex (,',mutex)
			 (let ((,',key (list ,@(mapcar #'car ,key-exprs))))
			   (multiple-value-bind (,',value ,',foundp)
			       (dic-ref ,'(,table) ,',key)
			     (if ,',foundp
				 ,',value
				 (progn
				   ,',(when debug `(dbg :memoize "Memoizing: id=~A, key=~A" ,id ,key))
				   (setf (dic-ref ,'(,table) ,',key)
					 (progn ,@,value-body))))))))))
	  ,@body)))))


;;;; Testing

(defun test-general ()
  (loop+ ((i (range :below 0))) do (error "should not be reached: ~D" i))
  (let ((x (make-array 10 :initial-element 0))
	(y (make-list 10 :initial-element 1)))
    (loop+ ((xc x) (yc y) i) doing
       (setf xc (+ i yc))
       finally (return x)))
  (keys (make-hash-table))
  (loop+ ((k '(1 2 3)) i) collecting (+ k i))
  (let ((dic (make-instance 'small-cache-dictionary :size 2)))
    (describe dic)
    (setf (dic-ref dic 1) 1)
    (describe dic)
    (print (dic-ref dic 1))
    (describe dic)
    (setf (dic-ref dic 2) 4)
    (describe dic)
    (print (dic-ref dic 2))
    (describe dic)
    (print (dic-ref dic 1))
    (describe dic)
    (setf (dic-ref dic 3) 9)
    (describe dic)
    (print (dic-ref dic 1))
    (describe dic)
    (dic-for-each (lambda (key value)
                    (print (cons key value)))
                  dic)
    (keys dic))

  (print
   (macroexpand '(with-memoization ()
		  (defun test (n)
		    (memoizing-let ((k (* 2 n)))
		      (sleep 1)
		      (* k k))))))
  (with-memoization ()
    (flet ((test (n)
	     (memoizing-let ((k (* 2 n))
			     (l (* 3 n)))
	       (sleep 1)
	       (* k l))))
      (time (test 5))
      (time (test 5))))
  
  (with-memoization (:type :local :size 2 :id 'test-memoization)
    (flet ((test (n)
	     (memoizing-let ((k (* 2 n))
			     (l (* 3 n)))
	       (sleep 1)
	       (* k l))))
      (time (test 5))
      (time (test 5))))

  (dbg-on :memoize)
  (with-memoization (:id 'test)
    (let* ((f (lambda (x) (sleep x) ))
	   (memoized-f (memoizing f)))
      (loop repeat 10 do (print (funcall memoized-f 1.0)))))
  (with-memoization (:id 'test)
    (let ((f (memoizing (lambda (x) (sleep x)))))
      (loop repeat 10 do (print (funcall f 1.0)))))
  (dbg-off)

  #-cmu
  (let ((x (make-instance 'property-mixin :properties (list :a 1))))
    (assert (= (slot-value x :a) 1))
    (setf (slot-value x :a) 3)
    (assert (= (slot-value x :a) 3)))

  (let ((x (make-instance 'property-mixin :properties '(:a 1)))
	(n 100))
    (time (loop repeat n do (get-property x :a)))
    #-cmu (time (loop repeat n do (slot-value x :a)))
    )
  (let ((x (make-instance 'property-mixin)))
    (describe x)
    (with-slots (x) x
      x))
  (dic-for-each-value #'print #(1 2 3))

  (let* ((count 0)
         (f (with-memoization (:type :global :id 'test)
              (memoizing (lambda (k)
                           (incf count)
                           (* k k))))))
    (loop repeat 100 do
         (fl.mp:make-thread
          (lambda ()
            (sleep (random 1.0))
            (funcall f (random 10)))))
    (sleep 2.0)
    count)

  )

;;; (test-general)
(fl.tests:adjoin-test 'test-general)