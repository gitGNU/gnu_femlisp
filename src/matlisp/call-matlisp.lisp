;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; call-matlisp.lisp - call external Matlisp functions
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;
;;; Copyright (C) 2007 Nicolas Neuss, University of Karlsruhe.
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

(in-package :fl.matlisp)

#+(or)(load (concatenate 'string cl-user::*cl-home* "matlisp/start.lisp"))

;;; (load-matlisp)

(defun complex-to-real-vector (cvec)
  (declare (type (simple-array (complex double-float) *) cvec))
  (let* ((n (length cvec))
	 (n2 (* 2 n)))
    (let ((dvec (make-double-vec n2)))
      (dotimes (i n)
	(let ((c (aref cvec i))
	      (k (* 2 i)))
	  (setf (aref dvec k) (realpart c))
	  (setf (aref dvec (1+ k)) (imagpart c))))
      dvec)))

(defun real-to-complex-vector (dvec)
  #+lispworks (declare (optimize (float 0)))
  (declare (type double-vec dvec))
  (very-quickly
    (let ((n2 (length dvec)))
      (assert (evenp n2))
      (let ((n  (/ n2 2)))
	(declare (type fixnum n))
	(let ((cvec (make-array n :element-type '(complex double-float))))
	  (dotimes (i n)
	    (setf (aref cvec i)
		  (let ((k (* 2 i)))
		    (complex (aref dvec k) (aref dvec (1+ k))))))
	  cvec)))))

(defun matlisp-symbol (string)
  (let ((matlisp-package (find-package "MATLISP")))
    (unless matlisp-package (error "Matlisp is not available"))
    (intern string matlisp-package)))

(defun fl-matlisp->matlisp (obj)
  (if (typep obj 'fl.matlisp:standard-matrix)
      (let ((element-type (element-type obj)))
	(cond
	  ((eq element-type 'double-float)
	   (make-instance (matlisp-symbol "REAL-MATRIX")
			  :nrows (nrows obj) :ncols (ncols obj) :store (store obj)))
	  ((equal element-type '(complex double-float))
	   (make-instance (matlisp-symbol "COMPLEX-MATRIX")
			  :nrows (nrows obj) :ncols (ncols obj) :store
			  (complex-to-real-vector (store obj))))
	  (t (error "Matlisp cannot handle ~A -matrices" element-type))))
      obj))

(defun matlisp->fl-matlisp (obj)
  (if (typep obj (matlisp-symbol "STANDARD-MATRIX"))
      (let ((element-type
	     (cond ((typep obj (matlisp-symbol "REAL-MATRIX")) 'double-float)
		   ((typep obj (matlisp-symbol "COMPLEX-MATRIX")) '(complex double-float))
		   (t (error "Unknown Matlisp matrix.")))))
	(make-instance (fl.matlisp:standard-matrix element-type)
		       :nrows (funcall (matlisp-symbol "NROWS") obj)
		       :ncols (funcall (matlisp-symbol "NCOLS") obj)
		       :store (let ((store (funcall (matlisp-symbol "STORE") obj)))
				(if (eq element-type 'double-float)
				    store
				    (real-to-complex-vector store)))))
      obj))

(defun translate-list (translator list)
  (let (table)
    (loop for obj in list
       unless (assoc obj table) do
	 (push (cons obj (funcall translator obj)) table))
    (mapcar (curry #'geta table) list)))

(defun matlisp (name &rest args)
  (let ((func (symbol-function (intern (string-upcase name) "MATLISP"))))
    (apply #'values
	   (translate-list
	    #'matlisp->fl-matlisp
	    (multiple-value-list
	     (apply func (translate-list #'fl-matlisp->matlisp args)))))))

;;; testing

(defun test-call-matlisp ()
  
  (unless (find-package "MATLISP")
    (format t "Matlisp is not loaded - skipping tests~%")
    (return-from test-call-matlisp))
  
  (matlisp->fl-matlisp (fl-matlisp->matlisp (eye 2)))
  (translate-list #'fl-matlisp->matlisp (list (eye 2)))
  (matlisp "EIG" (eye 4))
  )
