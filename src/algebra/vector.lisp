;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; vector.lisp - Vector class and interface
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

(in-package :algebra)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; double-vec
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(deftype double-vec () '(simple-array double-float (*)))

(definline list->double-vec (lst)
  (map 'double-vec
       #'(lambda (x) (if (typep x 'double-float) x (float x 1.0d0)))
       lst))

(definline double-vec (&rest comps)
  (list->double-vec comps))

(definline make-double-vec (dim &optional (init 0.0d0))
  "make-double-vec: double-vec constructor"
  (make-array dim :element-type 'double-float :initial-element (float init 0.0d0)))

(definline unit-vector (dim i)
  (let ((vec (make-double-vec dim)))
    (setf (aref vec i) 1.0d0)
    vec))

(defgeneric multiplicity (vec)
  (:documentation "We allow multiple vectors, for solving linear problems
in parallel."))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Vector operation interface
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; random access to cells
(defgeneric* vec-ref (vec index))
(defgeneric* (setf vec-ref) (val vec index))

;;; loop over indices and entries
(defgeneric* for-each-key (func arr))
(defgeneric* for-each-entry (func vec))
(defgeneric* for-each-key-and-entry (func arr))
(defgeneric* for-each-entry-of-vec1 (func arr1 arr2))
(defgeneric* for-each-entry-of-vec2 (func arr1 arr2))

(defmacro dovec ((loop-vars vec) &body body)
  "Usage: (dovec ((key) vec) ...)
          (dovec (entry vec) ...)
          (dovec ((key entry) vec) ...)"
  (let* ((loop-vars (if (consp loop-vars) loop-vars (list nil loop-vars)))
	 (vec-for-each (if (null (car loop-vars))
			   (if (null (cdr loop-vars))
			       (error "no loop variable")
			       'for-each-entry)
			   (if (or (null (cdr loop-vars)) (null (cadr loop-vars)))
			       'for-each-key
			       'for-each-key-and-entry))))
    `(,vec-for-each #'(lambda ,(remove nil loop-vars) ,@body) ,vec)))

(defmethod total-entries (obj)
  "Default method for total entries of vector or matrix is recursively
defined."
  (let ((entries 0))
    (for-each-entry
     #'(lambda (entry)
	 (incf entries (total-entries entry)))
     obj)
    entries))

(defmethod total-entries ((obj number)) 1)

;;; destructive operations
(defgeneric* x<-0 (x))
(defgeneric* x<-s (x s))
(defgeneric* x<-random (x s))

(defgeneric* x<-y (x y))
(defgeneric* x+=y (x y))
(defgeneric* x-=y (x y))
(defgeneric* copy! (x y))

(defgeneric* x<-s*y (x s y))
(defgeneric* x+=s*y (x s y))
(defgeneric* x-=s*y (x s y))

;;; non-destructive operations
(defgeneric* vec+ (vec1 vec2))
(defgeneric* vec-s* (val vec))
(defgeneric* vec- (vec1 vec2))
(defgeneric* vec-sp (vec1 vec2))
(defgeneric* zero-vector (x))

;;; copying
(defgeneric* make-analog (vec))  
(defgeneric* nr-of-entries (vec))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Default definitions
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Definitions for the double-vec class
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; PCL cannot dispatch on double-vec, so we use array.  But these methods
;;; are inlined by the compiler, thus no speed disadvantage should arise.

(with-fast-clos ()
(defmethod** vec-ref ((vec array) (index integer))
  (aref vec index))
(defmethod** (setf vec-ref) (val (vec array) (index integer))
  (setf (aref vec index) val))
)

(defmethod** for-each-key ((func function) (arr array))
  (dotimes (i (length arr)) (funcall func i)))
(defmethod** for-each-entry ((func function) (arr array))
  (cond ((vectorp arr)
	 (dotimes (i (length arr)) (funcall func (aref arr i))))
	(t (array-for-each func arr))))

(defmethod** for-each-key-and-entry ((func function) (arr array))
  (dotimes (i (length arr))
    (funcall func i (aref arr i))))

(defmethod** for-each-entry-of-vec1 ((func function) (arr1 array) (arr2 array))
  (if (vectorp arr1)
      (dotimes (i (length arr1))
	(funcall func (aref arr1 i) (aref arr2 i)))
      (array-for-each func arr1 arr2)))
(defmethod** for-each-entry-of-vec1 ((func function) (lst1 list) (lst2 list))
  (loop for el1 in lst1 and el2 in lst2 do
	(funcall func el1 el2)))
(defmethod** for-each-entry-of-vec2 ((func function) (arr1 array) (arr2 array))
  (if (vectorp arr2)
      (dotimes (i (length arr2))
	(funcall func (aref arr1 i) (aref arr2 i)))
      (array-for-each func arr1 arr2)))
(defmethod** for-each-entry-of-vec2 ((func function) (lst1 list) (lst2 list))
  (loop for el1 in lst1 and el2 in lst2 do
	(funcall func el1 el2)))


(defmethod** x<-0 ((x array))
  (fill x 0.0d0))
(defmethod** x<-s ((x array) (s double-float))
  (fill x s))
(defmethod** x<-random ((x array) s)
  (dotimes (k (length x) x)
    (setf (aref x k) (random s))))

(defmethod** x<-y ((x array) (y array))
  (dotimes (k (length x) x)
    (setf (aref x k) (aref y k))))
(defmethod** copy! ((x array) (y array))
  (dotimes (k (length x) x)
    (setf (aref y k) (aref x k))))
(defmethod** x+=y ((x array) (y array))
  (dotimes (k (length x) x)
    (incf (aref x k) (aref y k))))
(defmethod** x-=y ((x array) (y array))
  (dotimes (k (length x) x)
    (decf (aref x k) (aref y k))))

(defmethod** x<-s*y ((x array) (s number) (y array))
  (dotimes (k (length x) x)
    (setf (aref x k) (* s (aref y k)))))
(defmethod** x+=s*y ((x array) (s number) (y array))
  (if (vectorp x)
      (dotimes (k (length x) x)
	(incf (aref x k) (* s (aref y k))))
      (call-next-method)))
(defmethod** x-=s*y ((x array) (s number) (y array))
  (dotimes (k (length x) x)
    (decf (aref x k) (* s (aref y k)))))

;;; later we will need also the following for application to blocks
(defmethod** x<-s ((x array) (s array))
  (dotimes (k (length x) x)
    (setf (aref x k) (aref s k))))
(defmethod** x<-s*y ((x array) (s array) (y array))
  (dotimes (k (length x) x)
    (setf (aref x k) (* (aref s k) (aref y k)))))
(defmethod** x+=s*y ((x array) (s array) (y array))
  (dotimes (k (length x) x)
    (incf (aref x k) (* (aref s k) (aref y k)))))
(defmethod** x-=s*y ((x array) (s array) (y array))
  (dotimes (k (length x) x)
    (decf (aref x k) (* (aref s k) (aref y k)))))

(defmethod** vec+ ((vec1 array) (vec2 array))
  (map (type-of vec1) #'+ vec1 vec2))
(defmethod** vec-s* ((val number) (vec array))
  (map (type-of vec) #'(lambda (x) (* val x)) vec))
(defmethod** vec- ((vec1 array) (vec2 array))
  (map 'double-vec #'- vec1 vec2))

(defmethod** vec-sp ((vec1 array) (vec2 array))
  (reduce #'vec+ (map 'simple-vector #'vec-sp vec1 vec2)))

(defmethod** vec-sp ((vec1 list) (vec2 list))
  (reduce #'vec+ (mapcar #'vec-sp vec1 vec2)))

(defmethod** zero-vector ((x array))
  (map (type-of x) #'zero-vector x))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; the non-destructive operations make sense also for numbers
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmethod** vec+ ((x number) (y number)) (+ x y))
(defmethod** vec- ((x number) (y number)) (- x y))
(defmethod** vec-s* ((x number) (y number)) (* x y))
(defmethod** vec-sp ((x number) (y number)) (* x y))
(defmethod** zero-vector ((x number)) (coerce 0 (type-of x)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; recursive definitions are OK for some compound classes
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmethod x<-0 (x)  (for-each-entry #'x<-0 x))
(defmethod x<-s (x s) (for-each-entry (rcurry #'x<-s s) x))
(defmethod x<-random (x s)  (for-each-entry (rcurry #'x<-random s) x))

(defmethod x<-y (x y) (for-each-entry-of-vec2 #'x<-y x y))
(defmethod x+=y (x y) (for-each-entry-of-vec2 #'x+=y x y))
(defmethod x-=y (x y) (for-each-entry-of-vec2 #'x-=y x y))
(defmethod copy! (x y) (for-each-entry-of-vec1 #'copy! x y))

(defmethod x<-s*y (x (s number) y) (for-each-entry-of-vec2 #'(lambda (x y) (x<-s*y x s y)) x y))
(defmethod x+=s*y (x (s number) y) (for-each-entry-of-vec2 #'(lambda (x y) (x+=s*y x s y)) x y))
(defmethod x-=s*y (x (s number) y) (for-each-entry-of-vec2 #'(lambda (x y) (x-=s*y x s y)) x y))

(defmethod vec- (vec1 vec2)
  (vec+ vec1 (vec-s* -1 vec2)))
(defmethod vec-sp (vec1 vec2)
  (let ((result nil))
    (for-each-entry-of-vec1
     #'(lambda (x y)
	 (if result
	     (x+=y result (vec-sp x y))
	     (setf result (vec-sp x y))))
     vec1 vec2)
    result))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; BLAS method generation
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; The following is an old version of blas command generation which is
;;; used only for the handling of matlisp matrices in matlisp.lisp.

(defparameter *vector-blas-command-list*
  '((x<-0 x) (x<-s x s) (x<-y x y) (copy! x y) (x+=y x y) (x-=y x y)
    (x<-s*y x s y) (x+=s*y x s y) (x-=s*y x s y)))
       
(defun initialize-vector-blas-methods
  (vec-class vec-scalar blas-caller
   &key blas-1-caller (vector-transformer #'identity) (data-transformer #'identity))
  "Standard template for initializing the blas methods for the commands from
*vector-blas-command-list*.  For examples see the application to
matlisp-matrices below, and for sparse matrices in sparse.lisp.  The optional
blas-1-caller is used for example in the case of sparse vectors where we loop
over the y-entries instead of the x-entries."
  (funcall
   (compile
    nil
    `(lambda ()
       ,@(mapcar
	  #'(lambda (proc.args)
	      (destructuring-bind (proc . args) proc.args
		`(defmethod ,proc
		   ,(mapcar #'(lambda (obj)
				(ecase obj
				  ((x) (list 'x vec-class))
				  ((y) (list 'y vec-class))
				  ((s) (list 's vec-scalar))))
			    args)
		   (,@(if (and blas-1-caller (member 'y args))
			  (funcall blas-1-caller proc)
			  (funcall blas-caller proc))
		      ,@(mapcar #'(lambda (obj)
				    (ecase obj
				      ((x) (funcall vector-transformer 'x))
				      ((y) (funcall vector-transformer 'y))
				      ((s) (funcall data-transformer 's))))
				args)))))
	  *vector-blas-command-list*)))))

(defun test-vector ()
  (let ((x (make-double-vec 2))
	(y (make-double-vec 2 1.0d0)))
    (x<-0 x)
    (x+=y x y)
    (x+=s*y x 2.0d0 y)
    (x-=y x y)
    (x-=s*y x 3.0d0 y)
    (assert (= (aref x 0) -1.0d0)))
  (let ((x #(0.0d0 0.0d0)))
    (x+=s*y x 0.5d0 #(1.0d0 0.0d0))
    (x+=s*y x 0.5d0 #(0.0d0 1.0d0)))
  (let* ((x (unit-vector 2 0))
	 (y (unit-vector 2 1))
	 (z (vec+ x y)))
    (list x y z (vec-s* 0.5d0 z)))
  )

;;; (test-vector)
(tests::adjoin-femlisp-test 'test-vector)
