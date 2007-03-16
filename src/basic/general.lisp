;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; general.lisp - globally useful definitions
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

(defun property-set-p (object property)
  "Returns T if @arg{property} is found in the object's properties."
  (get-properties (slot-value object 'properties) (list property)))

(defun get-property (object property)
  "Gets @arg{property} for @arg{object}."
  (getf (slot-value object 'properties) property))

(defun (setf get-property) (value object property)
  "Sets the property @arg{property} of @arg{problem} to @arg{value}."
  (setf (getf (slot-value object 'properties) property) value))

(defmacro with-properties (properties object &body body)
  "Work with @arg{properties} on the property list of @arg{object}."
  (with-gensyms (obj)
    `(let ((,obj ,object))
       (symbol-macrolet ,(loop for prop in properties collect
			      `(,prop (getf (properties ,obj) ',prop)))
	 ,@body))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Runtime compilation
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun runtime-compile (source)
  "Calls compile on the provided @arg{source}.  When :compile is activated
for debugging, the source code is printed."
  (let ((*print-circle* nil))
    (dbg :compile "Compiling source: ~%~S~%" source))
  (funcall (if (dbg-p :compile) #'compile #'fl.port:compile-silently)
	   nil source))

(defun compile-and-eval (source)
  "Compiles and evaluates the given @arg{source}.  This should be an ANSI
compatible way of ensuring method compilation."
  (dbg :compile "Compiling and evaluating: ~%~S~%" source)
  (funcall (funcall (if (dbg-p :compile) #'compile #'fl.port:compile-silently)
		    nil `(lambda () ,source))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; iteration
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
    (with-slots (from to below) range
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

(defclass cache-dictionary ()
  ((size :reader size :initform 1 :initarg :size)
   (test :initform #'eql :initarg :test)
   (store :reader store :initform ())))

(defmethod initialize-instance :after ((dic cache-dictionary) &key &allow-other-keys)
  (assert (plusp (size dic))))

(defmethod dic-ref ((dic cache-dictionary) key)
  (declare (optimize speed))
  (with-slots (store test) dic
    (loop for tail on store
	  and prev-tail = nil then tail do
	  (when (funcall test (caar tail) key)
	    (when prev-tail
	      (setf (cdr prev-tail) (cdr tail)
		    (cdr tail) store
		    store tail))
	    (return (values (cdar store) t))))))

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
  (loop+ ((i (range :below 0))) do (error "should not be reached"))
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
  )

;;; (test-general)
(fl.tests:adjoin-test 'test-general)