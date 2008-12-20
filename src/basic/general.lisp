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
  (:documentation "Generate an analogous but empty data structure."))

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
  (:method ((table hash-table) key)
    (gethash key table)))

(defgeneric (setf dic-ref) (value dictionary key)
  (:documentation "Sets the entry @arg{key} in @arg{dictionary} to
@arg{value}.")
  (:method (value (table hash-table) key)
    (setf (gethash key table) value)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Memoization
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; Cache-dictionary list

(defclass cd-item ()
  ((key :initarg :key)
   (value :initarg :value)
   (pred :initform nil)
   (succ :initform nil)))

(defclass cd-list ()
  ((first :initform nil)
   (last :initform nil)))

(defmethod emptyp ((cdl cd-list))
  (with-slots (first) cdl
    (null first)))

(defmethod push-front ((cdi cd-item) (cdl cd-list))
  (with-slots (first last) cdl
    (with-slots (pred succ) cdi
      (setf succ first pred nil))
    (if first
	(setf (slot-value first 'pred) cdi)
	(setf last cdi))
    (setf first cdi)))

(defmethod delete-from ((cdi cd-item) (cdl cd-list))
  (with-slots (pred succ) cdi
    (with-slots (first last) cdl
      (if pred
	  (setf (slot-value pred 'succ) succ)
	  (setf first succ))
      (if succ
	  (setf (slot-value succ 'pred) pred)
	  (setf last pred)))
    (setf pred nil succ nil))
  cdi)

(defmethod pop-front ((cdl cd-list))
  (lret ((first (slot-value cdl 'first)))
    (with-slots (pred succ) first
      (setf (slot-value cdl 'first) succ)
      (if succ
	  (setf (slot-value succ 'pred) nil)
	  (setf (slot-value cdl 'last) nil))
      (setf pred nil succ nil))))

(defmethod pop-rear ((cdl cd-list))
  (lret ((last (slot-value cdl 'last)))
    (with-slots (pred succ) last
      (setf (slot-value cdl 'last) pred)
      (if pred
	  (setf (slot-value pred 'succ) nil)
	  (setf (slot-value cdl 'first) nil))
      (setf pred nil succ nil))))

(defgeneric coerce-to (seq type)
  (:documentation "Generic version of @arg{coerce}.")
  (:method (seq type) (coerce seq type)))

(defmethod coerce-to ((cdl cd-list) type)
  (coerce (loop for item = (slot-value cdl 'first)
	     then (slot-value item 'succ) while item
	     collect (slot-value item 'object))
	  type))

;;; Cache dictionary

(defclass cache-dictionary ()
  ((size :reader size :initform 1 :initarg :size)
   (test :initform #'eql :initarg :test)
   (store :reader store :initform ())))

(defmethod initialize-instance :after ((dic cache-dictionary) &key &allow-other-keys)
  (assert (plusp (size dic))))

(defmethod dic-ref ((dic cache-dictionary) key)
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

(defmethod (setf dic-ref) (value (dic cache-dictionary) key)
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

(defmacro memoizing-let (bindings &body body)
  "The @arg{body} is memoized for the keys given as a list of
@arg{bindings}."
  (declare (ignore bindings body))
  (error "This global macro should be used only within the context of
@macro{with-memoization}."))

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
	   (,table
	    ,(if size
		 `(make-instance 'cache-dictionary :size ,size :test ,test)
		 `(make-hash-table :test ,test))))
      (declare (ignorable ,mutex ,table))
      (flet ((,table ()
	       ,(ecase type
		       (:global table)
		       (:local
			`(or (getf fl.multiprocessing:*thread-local-memoization-table* ,id)
			  (setf (getf fl.multiprocessing:*thread-local-memoization-table* ,id)
			   ,(if size
				`(make-instance 'cache-dictionary :size ,size :test ,test)
				`(make-hash-table :test ,test))))))))
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
  (loop+ ((k '(1 2 3)) i) collecting (+ k i))
  (let ((dic (make-instance 'cache-dictionary :size 2)))
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
    (describe dic))
  
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
  
  )

(with-memoization (:test 'equalp)
  (defun lagrange-fe (order &key (nr-comps 1) (type :uniform))
    (memoizing-let ((order order) (nr-comps nr-comps) (type type))
      (list order nr-comps type))))

;;; (test-general)
(fl.tests:adjoin-test 'test-general)