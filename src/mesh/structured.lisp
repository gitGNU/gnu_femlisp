;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; structured.lisp - Not in use yet
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


(deftype int () '(signed-byte 32))
(deftype uint () '(unsigned-byte 32))
(deftype int-vec () '(simple-array int (*)))
(deftype uint-vec () '(simple-array uint (*)))
(deftype positive-fixnum () '(and fixnum unsigned-byte))
(deftype double-vec () '(simple-array double-float (*)))  
(deftype fixnum-vec () '(simple-array fixnum (*)))

(declaim (inline make-double-vec))
(defun make-uint-vec (dim &optional (init 0))
  (make-array dim :element-type 'uint :initial-element init))
(defun make-fixnum-vec (dim &optional (init 0))
  (make-array dim :element-type 'fixnum :initial-element init))
(defun make-double-vec (dim &optional (init 0.0d0))
  "make-double-vec: double-vec constructor"
  (make-array dim :element-type 'double-float :initial-element init))
(defun fixnum-vec (&rest cells)
  (coerce cells 'fixnum-vec))

(defclass structured-grid ()
  ((dimensions :accessor dimensions :initarg :dimensions :type fixnum-vec)
   (offset0 :accessor offset0 :initarg :offset0 :initform 0 :type fixnum)
   (offsets :accessor offsets :initarg :offsets :type fixnum-vec)))

(defun vec-to-pos (vec grid)
  (loop for shift across vec and
	offset across (offsets grid)
	summing (* shift offset)))

(defclass structured-grid-vector ()
  ((grid :accessor grid :initarg :grid :type structured-grid)
   (entries :accessor entries :initarg :entries :type double-vec)))

(defclass stencil ()
  ((entries :accessor entries :initarg :entries :type list))
  (:documentation "A stencil line consists of a double-value followed by a series of relative
offsets.  E.g. the usual 5-point stencil would have the entries
  '((4.0d0 #(0 0)) (-1.0d0 #(0 -1) #(-1 0) #(1 0) #(0 1)))"))

(defun normalize-entries (entries)
  "Throws out zeros, collects equal entries."
  (loop with result = ()
	for entry in entries
	for value = (car entry)
	for result-entry = (assoc value result)
	unless (zerop value) do
	(if result-entry
	    (setf (cdr result-entry) (union (cdr entry) (cdr result-entry)))
	    (push entry result))
	finally (return result)))


(defun make-stencil (entries)
  (make-instance
   'stencil
   :entries (normalize-entries entries)))

(defun stencil-application (stencil grid pos values result-values)
  `(setf (aref ,result-values pos)
    (+
     ,@(mapcar
	#'(lambda (entry-line)
	    `(* ,(car entry-line)
	      (+ ,@(loop for vec in (cdr entry-line)
			 collect `(aref ,values (the fixnum (+ ,pos ,(vec-to-pos vec grid))))))))
	(entries stencil)))))

(stencil-application
 (make-stencil '((-0.25d0 #(-1 -1) #(-1 0) #(-1 1) #(0 -1) #(0 0) #(0 1) #(1 -1) #(1 0) #(1 1))))
 (make-instance 'structured-grid :dimensions #(1000 1000) :offsets #(1 1000))
 'pos
 'entries
 'entries)


1. zuviel Typdeklaration

2. Hauptproblem: make-double-vec war nicht inline deklariert!!

3. vielleich (compilation-speed 0)


(defun test2 ()
  (let* ((n 1000)
	 (mat (make-array (* n n) :element-type 'double-float
			  :initial-element 1.0d0)))
    (declare (type (simple-array double-float (*)) mat)
	     (fixnum n)
	     (optimize (speed 3) (debug 0) (safety 0)))
    (loop for i of-type fixnum below 10 do
	    (loop for pos1 of-type fixnum from 1 below (1- n) do
		    (loop for pos2 of-type fixnum from (+ pos1 n)
			  below (the fixnum (+ pos1 (the fixnum (* (the fixnum (1- n)) n))))
			  by n do
			    (let ((a0 (aref mat (the fixnum (- pos2 1001))))
				  (a1 (aref mat (the fixnum (- pos2 1))))
				  (a2 (aref mat (the fixnum (+ pos2 999))))
				  (a3 (aref mat (the fixnum (- pos2 1000))))
				  (a4 (aref mat (the fixnum (- pos2 0))))
				  (a5 (aref mat (the fixnum (+ pos2 1000))))
				  (a6 (aref mat (the fixnum (- pos2 999))))
				  (a7 (aref mat (the fixnum (+ pos2 1))))
				  (a8 (aref mat (the fixnum (+ pos2 1001)))))
			      (setf (aref mat pos2)
				    (* 0.1111111111d0
				       (+ a0 a1 a2 a3 a4 a5 a6 a7 a8)))))))))
(defun test ()
  (let* ((dim 2) (n 1000)
	 (entries (make-double-vec (expt n dim) 1.0d0)))
    (declare (optimize (speed 3) (safety 0) (debug 0) (compilation-speed 0)))
    (loop for pos1 of-type fixnum from 1 below (- n 1)
	  do
	  (loop for pos2 of-type fixnum from (+ pos1 n) by n below (+ pos1 (* n (- n 1)))
		do
		(setf (aref entries pos2)
		      (let ((t0 (aref entries (+ pos2 -1001)))
			    (t1 (aref entries (+ pos2 -1)))
			    (t2 (aref entries (+ pos2 999)))
			    (t3 (aref entries (+ pos2 -1000)))
			    (t4 (aref entries (+ pos2 0)))
			    (t5 (aref entries (+ pos2 1000)))
			    (t6 (aref entries (+ pos2 -999)))
			    (t7 (aref entries (+ pos2 1)))
			    (t8 (aref entries (+ pos2 1001))))
			(* 0.1111111111 (+ t0 t1 t2 t3 t4 t5 t6 t7 t8))))))))

(defun test ()
  (let* ((dim 2) (n 1000)
	 (entries (make-double-vec (expt n dim) 1.0d0)))
    (declare (optimize (speed 3) (safety 0) (debug 0) (compilation-speed 0)))
    (dotimes (i 10)
      (loop
       for pos1 of-type fixnum from 1 below (- n 1) do
       (loop
	for pos2 of-type fixnum from (+ pos1 n) by n below (+ pos1 (* n (- n 1)))
	do
	#+(and)
	(setf (aref entries pos2)
	      (+ (* 0.1111111111d0 (aref entries (+ pos2 -1001)))
		 (* 0.1111111111d0 (aref entries (+ pos2 -1)))
		 (* 0.1111111111d0 (aref entries (+ pos2 999)))
		 (* 0.1111111111d0 (aref entries (+ pos2 -1000)))
		 (* 0.1111111111d0 (aref entries (+ pos2 0)))
		 (* 0.1111111111d0 (aref entries (+ pos2 1000)))
		 (* 0.1111111111d0 (aref entries (+ pos2 -999)))
		 (* 0.1111111111d0 (aref entries (+ pos2 1)))
		 (* 0.1111111111d0 (aref entries (+ pos2 1001)))))
	#+(or)
	(setf (aref entries pos2)
	      (* 0.1111111111d0
		 (+ (aref entries (+ pos2 -1001))
		    (aref entries (+ pos2 -1))
		    (aref entries (+ pos2 999))
		    (aref entries (+ pos2 -1000))
		    (aref entries (+ pos2 0))
		    (aref entries (+ pos2 1000))
		    (aref entries (+ pos2 -999))
		    (aref entries (+ pos2 1))
		    (aref entries (+ pos2 1001))))))))))
(disassemble 'test)

(time (test))

(defun apply-stencil (stencil sgv)
  (let* ((sg (grid sgv))
	 (entries sgv)
	 (dims (dimensions grid))
	 (offsets (dimensions grid))
	 (offset0 (offset0 grid))
	 (apply-locally (apply-locally stencil grid)))
    
    ;; loop through internal positions
    (let ((offsets (make-fixnum-vec 1 1))
	  (dims (make-fixnum-vec 1 10))
	  (entries (make-double-vec 10))
	  (offset0 0))
    (declare (type fixnum offset0)
	     (type fixnum-vec dims offsets)
	     (type double-vec entries))
    (loop with offset of-type fixnum = (aref offset 0)
	  for pos of-type fixnum from (the fixnum (+ offset0 offset)) by offset
	  until (the fixnum (+ offset0 (the fixnum (* offset (aref dims 0)))))
	  summing (aref entries pos)))
    
	  (loop with offset of-type fixnum = (aref offset 1)
		for pos of-type fixnum from (the fixnum (+ offset0 offset1)) by offset1
		until (the fixnum (+ offset0 (the fixnum (* offset1 (aref dims 0)))))
		do
	  (loop with offset of-type fixnum = (aref offset 1)
		initially (incf pos offset)
		repeat (- (aref dimensions 0) 2) do
		(funcall apply-stencil-locally entries pos)
		finally (incf pos offset))
	  
	  (loop for j from 1  below (1- (aref dimensions 1)) do
		
	  (apply-filter filter entries i)))))

(defun index->offset (index offsets)
  (declare (values fixnum))
  (loop for i across index
	and offset across offsets
	sum (* i offset)))

(defun apply-stencil-at-interior-pos-function (stencil offsets entries i)
  `(setf (aref ,entries ,i)
    (+ ,@(loop for stencil-line in stencil collecting
	       `(*
		 ,(car filter-line)
		 (+
		  ,@(loop for index in (cdr filter-line) collecting
			  `(aref ,entries ,(+ i (index->offset index offsets))))))))))

(apply-filter-at-pos '((4.0d0 #(0 0)) (-1.0d0 #(0 -1) #(-1 0) #(1 0) #(0 1)))
		     #(1 10)
		     'entries
		     8)

