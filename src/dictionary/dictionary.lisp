;;; -*- mode: lisp; fill-column: 75; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; __FILENAME__.lisp - Sample header for a Femlisp source file
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;
;;; Copyright (C) 2003-2006 Nicolas Neuss, University of Heidelberg.
;;; Copyright (C) 2006-2011 Nicolas Neuss, KIT Karlsruhe.
;;; Copyright (C) 2011-
;;; Nicolas Neuss, Friedrich-Alexander-Universitaet Erlangen-Nuernberg
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

(in-package :fl.dictionary)

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

(defgeneric dic-for-each (function dictionary &key &allow-other-keys)
  (:documentation "Applies @arg{function} to each key-value pair in
  @arg{dictionary}.")
  (:method :around (function dictionary &key parallel &allow-other-keys)
    (if (and parallel lparallel:*kernel*)
        (with-workers (function)
          ;;; we do not want nested parallelization!
          (call-next-method #'work-on dictionary :parallel nil))
        (call-next-method)))
  (:method (func (dic list) &key &allow-other-keys)
    (loop for (key . value) in dic
              do (funcall func key value)))
  (:method (func (dic hash-table) &key &allow-other-keys)
    (maphash func dic))
  (:method (func (vec vector) &key &allow-other-keys)
    (loop for i from 0 and x across vec do
      (funcall func i x))))
  
(defgeneric dic-for-each-key (function dictionary &key &allow-other-keys)
  (:documentation "Applies @arg{function} to each key in @arg{dictionary}.")
  (:method (func dic &rest args &key &allow-other-keys)
    (apply #'dic-for-each
           (lambda (key value)
             (declare (ignore value))
             (funcall func key))
           dic args)))

(defgeneric dic-for-each-value (function dictionary &key &allow-other-keys)
  (:documentation "Applies @arg{function} to each value in @arg{dictionary}.")
  (:method (func dic &rest args)
    (apply #'dic-for-each
           (lambda (key value)
             (declare (ignore key))
             (funcall func value))
           dic args)))

(defgeneric dic-remove (key dictionary)
  (:documentation "Removes @arg{key} from @arg{dictionary}.")
  (:method (key (dic hash-table))
    (remhash key dic)))

(defgeneric dic-empty-p (dictionary)
  (:documentation "Tests if @arg{dictionary} is empty.")
  (:method (dic)
    (not (mapper-some (constantly t) dic)))
  (:method ((list list))
    (null list))
  (:method ((dic hash-table))
    (zerop (hash-table-count dic))))

;;; derived
(defgeneric keys (dic)
  (:documentation "Returns a list of all keys of @arg{dic}.")
  (:method (dic)
    (mapper-collect #'dic-for-each-key dic)))

(defmacro dodic ((looping-var dic &key parallel) &body body)
  "Loops through @arg{dic}.  If @arg{looping-var} is an atom
@emph{key}, loop through the keys; if it is a list of the form
@emph{(value)} loop through the values; if it is a list of the form
@emph{(key value)} loop through key and value."
  (let ((key (if (single? looping-var)
                 (gensym)
                 (if (atom looping-var) looping-var (car looping-var))))
        (value (if (atom looping-var)
                   (gensym)
                   (car (last looping-var)))))
    `(block nil
       (dic-for-each
        (lambda (,key ,value)
          ,@(when (atom looping-var) `((declare (ignore ,value))))
          ,@(when (single? looping-var) `((declare (ignore ,key))))
          ,@body)
        ,dic :parallel ,parallel))))

(defun get-key-from-dic (dic)
  "Returns a key from dic if nonempty."
  (unless (dic-empty-p dic)
    (dodic (key dic)
      (return-from get-key-from-dic
        (values key t)))))
      
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Dictionaries and memoization
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass dictionary ()
  ()
  (:documentation "A superclass for user-defined dictionaries."))

(defclass sorted-dictionary (dictionary)
  ()
  (:documentation "A superclass for sorted dictionaries."))

(defgeneric dic-pop (sorted-dic)
  (:documentation "A pop routine for sorted dictionaries or heaps."))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; small-cache-dictionary
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass small-cache-dictionary (dictionary)
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

(defmethod dic-for-each (func (dic small-cache-dictionary) &rest args)
  (apply #'dic-for-each func (store dic) args))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; sorted-hash-table
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass sorted-hash-table (sorted-dictionary)
  ((test :initform #'eql :initarg :test)
   (store :initform (make-dll))
   (insertion-order :initform :FIFO :initarg :insertion-order)
   (item-store))
  (:documentation "This is a hash-table where new entries are interned
  FIFO-like (insertion-order=:fifo) or
  heap-like (insertion-order=:heap)."))

(defmethod initialize-instance :after
    ((dic sorted-hash-table) &key &allow-other-keys)
  (with-slots (test item-store size) dic
    (setf item-store (make-hash-table :test test))))

(defmethod dic-ref ((dic sorted-hash-table) key)
  (with-slots (store item-store) dic
    (whereas ((item (gethash key item-store)))
      (values (cdr (dli-object item)) t))))

(defmethod (setf dic-ref) (value (dic sorted-hash-table) key)
  (with-slots (store item-store size test insertion-order) dic
    (let ((item (gethash key item-store)))
      (if item
          (setf (cdr (dli-object item)) value)
          (setf (gethash key item-store)
                (funcall (ecase insertion-order
                           (:fifo #'dll-rear-insert)
                           (:heap #'dll-front-insert))
                         (cons key value) store)))))
  value)

(defmethod dic-for-each (func (dic sorted-hash-table) &key (direction :forward) &allow-other-keys)
  (with-slots (store) dic
    (dll-for-each (_ (funcall func (car _) (cdr _))) store direction)))

(defmethod dic-remove (key (dic sorted-hash-table))
  (with-slots (store item-store) dic
    (awhen (gethash key item-store)
      (dll-remove-item it store)
      (remhash key item-store))))

(defmethod dic-empty-p ((dic sorted-hash-table))
  (with-slots (store) dic
    (dll-empty-p store)))

(defmethod dic-pop ((dic sorted-hash-table))
  (with-slots (store item-store) dic
    (destructuring-bind (key . val)
        (dll-pop-first store)
      (remhash key item-store)
      (values key val))))

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

(defmethod dic-for-each-key (func (dic computed-value-dictionary)
                             &key &allow-other-keys)
  (for-each func (keys dic)))

(defmethod dic-for-each (func (dic computed-value-dictionary)
                         &key &allow-other-keys)
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
  (warn "@macro{memoizing-let} should only be used within the context of
@macro{with-memoization}.")
  `(let ,bindings ,@body))

(defmacro memoizing-let* (bindings &body body)
  "The @arg{body} is memoized for the keys given as a list of
@arg{bindings}."
  (warn "@macro{memoizing-let*} should only be used within the context of
@macro{with-memoization}.")
  `(let ,bindings ,@body))

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

;;; non-thread-safe memoization

#+(or)
(defmacro with-memoization ((&key type size id (debug t) (test ''equal)) &body body)
  "Sets up a memoization environment consisting of a table, and a captured
symbol @symbol{memoizing-let} memoizing its body depending on the
arguments.  Example of usage:
@lisp
  (with-memoization (:size 4 :id 'test)
    (defun test (n)
      (memoizing-let ((k (* 2 n)))
	(sleep 1)
	(* k k))))
@end lisp

Note that this definition is later -in the file
@path{parallel-adaptions.lisp}- superseded by one
handling also thread-local memoization."
  (declare (ignore type))
  (with-gensyms (table key key-exprs value-body func foundp value)
    `(let ((,table (memoization-table ,size ,test)))
       ;; the memoizing-let symbol is captured, 
       (macrolet ((memoizing (,func)
                    `(lambda (&rest ,',key)
                       (multiple-value-bind (,',value ,',foundp)
                           (dic-ref ,',table ,',key)
                         (if ,',foundp
                             ,',value
                             (progn
                               ,',(when debug `(dbg :memoize "Memoizing: id=~A, key=~A" ,id ,key))
                               (setf (dic-ref ,',table ,',key)
                                     (apply ,,func ,',key)))))))
                  (memoizing-let-generator (let-symbol ,key-exprs &body ,value-body)
                    (let ((,key-exprs
                           (mapcar (lambda (entry)
                                     (if (symbolp entry)
                                         (list entry entry)
                                         entry))
                                   ,key-exprs)))
                      `(,let-symbol ,,key-exprs
                         ;; declares should be inserted here from value-body
                         (let ((,',key (list ,@(mapcar #'car ,key-exprs))))
                           (multiple-value-bind (,',value ,',foundp)
                               (dic-ref ,',table ,',key)
                             (if ,',foundp
                                 ,',value
                                 (progn
                                   ,',(when debug `(dbg :memoize "Memoizing: id=~A, key=~A" ,id ,key))
                                   (setf (dic-ref ,',table ,',key)
                                         (progn ,@,value-body)))))))))
                  (memoizing-let (&rest rest) `(memoizing-let-generator let ,@rest))
                  (memoizing-let* (&rest rest) `(memoizing-let-generator let* ,@rest))
                  )
         ,@body))))

;;; Thread-safe memoization

(defparameter *thread-local-memoization-table* nil
  "Dynamic variable specialized for each worker thread separately.  Should
not be used globally.")

(defmacro with-memoization ((&key (type :global) size id debug (test ''equal)) &body body)
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
    `(let ((,mutex (bordeaux-threads:make-recursive-lock))
	   (,table ,(ecase type
                           (:global `(memoization-table ,size ,test))
                           (:local
                            `(ensure (getf *thread-local-memoization-table* ,id)
                                     (memoization-table ,size ,test))))))
	;; the memoizing-let symbol is captured, 
       (macrolet ((memoizing (,func)
                    `(lambda (&rest ,',key)
		       (bordeaux-threads:with-recursive-lock-held (,',mutex)
			 (multiple-value-bind (,',value ,',foundp)
			     (dic-ref ,',table ,',key)
			   (if ,',foundp
			       ,',value
			       (progn
				 ,',(when debug `(dbg :memoize "Memoizing: id=~A, key=~A" ,id ,key))
				 (setf (dic-ref ,',table ,',key)
				       (apply ,,func ,',key))))))))
                  (memoizing-let-generator (let-symbol ,key-exprs &body ,value-body)
                    (let ((,key-exprs
                            (mapcar (lambda (entry)
                                      (if (symbolp entry)
                                          (list entry entry)
                                          entry))
                                    ,key-exprs)))
                      `(,let-symbol ,,key-exprs
                         ;; declares should be inserted here from value-body
                         (bordeaux-threads:with-recursive-lock-held (,',mutex)
                           (let ((,',key (list ,@(mapcar #'car ,key-exprs))))
                             (multiple-value-bind (,',value ,',foundp)
                                 (dic-ref ,',table ,',key)
                               (if ,',foundp
                                   ,',value
                                   (progn
                                     ,',(when debug `(dbg :memoize "Memoizing: id=~A, key=~A" ,id ,key))
                                     (setf (dic-ref ,',table ,',key)
                                           (progn ,@,value-body))))))))))
                  (memoizing-let (&rest rest) `(memoizing-let-generator let ,@rest))
                  (memoizing-let* (&rest rest) `(memoizing-let-generator let* ,@rest))
                  )
         ,@body))))


;;;; Testing

(defun test-general ()
  (let ((table (make-instance 'sorted-hash-table :insertion-order :heap)))
    (setf (dic-ref table 'a) 1)
    (setf (dic-ref table 'b) 2)
    (setf (dic-ref table 'c) 3)
    (dic-ref table 'b)
    (dic-for-each (_ (format t "~A->~A~%" _1 _2)) table
                  :direction :forward))
  (keys (make-hash-table))
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
	     (memoizing-let (n
                             (k (* 2 n))
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

  (dic-for-each-value #'print #(1 2 3))

  (if lparallel:*kernel*
      (let* ((count 0)
             (f (with-memoization (:type :global :id 'test)
                  (memoizing (lambda (k)
                               (sleep 1.0)
                               (incf count)
                               (* k k))))))
        (fl.parallel:pwork
         f (make-array (lparallel:kernel-worker-count) :initial-element '(1)))
        (time (funcall f 1))
        count))

  (with-memoization (:type :local :size 2 :id 'test-memoization)
    (flet ((test (n)
	     (memoizing-let ((k (* 2 n))
			     (l (* 3 n)))
	       (sleep 1)
	       (* k l))))
      (time (test 5))
      (time (test 5))))
  
  (with-memoization (:size 2)
    (flet ((test (n)
	     (memoizing-let* ((k (* 2 n))
                              (l (* k n)))
	       (sleep 1)
	       (* k l))))
      (time (test 5))
      (time (test 5))))

  (let* ((count 0)
         (f (with-memoization (:type :global :id 'test)
              (memoizing (lambda (k)
                           (incf count)
                           (* k k))))))
    (loop repeat 100 do
         (bordeaux-threads:make-thread
          (lambda ()
            (sleep (random 1.0))
            (funcall f (random 10)))))
    (sleep 2.0)
    count)
  )

;;; (test-general)
(fl.tests:adjoin-test 'test-general)
