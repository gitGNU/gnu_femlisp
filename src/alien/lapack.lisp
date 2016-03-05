;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; lapack.lisp - generating interfaces to LAPACK functions
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
;;; NO EVENT SHALL THE AUTHOR, THE UNIVERSITY OF KARLSRUHE, OR OTHER
;;; CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
;;; EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
;;; PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
;;; PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
;;; LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
;;; NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
;;; SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(in-package :cl-user)

(defpackage :fl.lapack
  (:use :common-lisp :fl.utilities :fl.port :fl.debug :fl.macros)
  (:export "CL->LAPACK-TYPE" "CALL-LAPACK" "CALL-LAPACK-WITH-ERROR-TEST" "LAPACK"))

(in-package :fl.lapack)

(defparameter *lapack-types*
  '(:float :double :complex-float :complex-double)
  "The number types for which a LAPACK routine is defined.")

(defun lapack-type-p (type)
  (member type *lapack-types*))

(defun number-type (lapack-type)
  (ecase lapack-type
    ((:float :double) :real)
    ((:complex-double :complex-float) :complex)))

(defun base-type (lapack-type)
  (ecase lapack-type
    ((:float :complex-float) :float)
    ((:double :complex-double) :double)))

(defun convert-to-alien-type (type)
  "Converts @arg{x} to an alien type, if possible."
  (convert-type (if (lapack-type-p type)
		    (base-type type)
		    type)))

(defun lapack-available-p ()
  (and (member :blas *features*)
       (member :lapack *features*)))

(defun cl->lapack-type (type &optional (error-p t))
  "Converts a CL type to a LAPACK type, if possible."
  (and (lapack-available-p)
       (cond
         ((eql type 'double-float) :double)
         ((eql type 'single-float) :float)
         ((equal type '(complex double-float)) :complex-double)
         ((equal type '(complex single-float)) :complex-float)
         (t (when error-p (error "Unknown type"))))))

(defun ensure-lapack-type (type &optional (error-p t))
  (if (lapack-type-p type)
      type
      (cl->lapack-type type error-p)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; LAPACK functions
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; The following approach is not general enough for some set of
;;; BLAS/LAPACK routines, because function names can vary more than that:
;;; e.g. DSPGV and ZHEGV for symmetric and hermitian generalized eigenvalue
;;; problems.

(defparameter *underscore-p*
  #+allegro nil #-allegro t
  "Does the external name contain an underscore.")

(defun lapack-external-name (name type &optional (underscore-p *underscore-p*))
  "@arg{name} is either a string denoting the LAPACK name without the
number-type prefix, or a list of two strings denoting the names for the
real and complex routines without the number-type prefix.  This is useful
when we have separate but similar routines for the real and complex case
like, e.g., the LAPACK routines SPGV and HEGV."
  (assert (member type *lapack-types*))
  (if (atom name)
      (format nil "~A~A~A"
	      (ecase type
		(:float "s")
		(:double "d")
		(:complex-float "c")
		(:complex-double "z"))
	      name
              (if underscore-p "_" ""))
      (format nil "~A~A"
              (nth (position type *lapack-types*) name)
              (if underscore-p "_" ""))))

(defun lapack-lisp-name (name type)
  (intern
   (string-upcase (lapack-external-name name type t))
   "FL.LAPACK"))

(defun lapack (name type)
  "Ensures the CL binding to the specified LAPACK function."
  (and (lapack-available-p)
       (let* ((type (ensure-lapack-type type))
              (symbol (lapack-lisp-name name type)))
         (unless (fboundp symbol)
           (create-lapack-function name type))
         (values (symbol-function symbol) symbol))))

(defun remove-all-lapack-functions (name)
  "Removes the CL bindings to the LAPACK functions categorized by
@arg{name}."
  (loop for type in *lapack-types*
     for symbol = (lapack-lisp-name name type)
     when (fboundp symbol) do
       (fmakunbound symbol)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; LAPACK templates
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defparameter *lapack-templates* (make-hash-table :test 'equal)
  "The table of all LAPACK templates.")

(defun insert-lapack-template (name return-value argument-specifications)
  (setf (gethash name *lapack-templates*)
	(cons return-value argument-specifications)))
       
(defmacro define-lapack-template (name return-value &body argument-specifications)
  "Defines a LAPACK template.  The syntax is similar to the SBCL alien
interface, but keywords are used for specifying types.  Furthermore, the
type of the matrix elements is not yet specified, but can be refered to
as :element-type, and its base-type as :base-type.  There is also a :select
clause which allows to include arguments only for specific
situations.  :select can dispatch
on :element-type (:float, :double, :complex-float, :complex-double),
:base-type (:float, :double),
:number-type (:real, :complex)"
  `(progn
     (remove-all-lapack-functions ,name)
     (insert-lapack-template ,name ',return-value ',argument-specifications)))

(defun translate-types (spec type)
  (let ((base-type (base-type type))
	(number-type (number-type type)))
    (map-tree (lambda (x)
		(cond ((eq x :element-type) type)
		      ((eq x :base-type) base-type)
		      ((eq x :number-type) number-type)
		      (t x)))
	      spec)))

(defun extract-applicable-spec (arg-spec type)
  "Handles the @arg{case}-construct and converts :element-type to
@arg{type}."
  ;; convert :element-type and :base-type
  (setq arg-spec (translate-types arg-spec type))
  ;; handle select-statement
  (destructuring-bind (first second . others)
      arg-spec
    (if (eq first :select)
	(loop for (spec clause) in others
	   when (if (atom spec)
		    (eq second spec)
		    (member second spec))
	   do (return clause))
	arg-spec)))

(defun check-spec (arg-spec)
  "Checks if arg-spec makes sense for the implementation.  Should probably
be done after replacing types for the implementation?"
  (destructuring-bind (name type mode) arg-spec
    (declare (ignore name))
    (when (or (and (consp type) (not (eq mode :in))))
      (error "Bad argument specification.")))
  arg-spec
  )

(defun create-lapack-function (name type)
  (let ((template (gethash name *lapack-templates*)))
    (unless template
      (error "No LAPACK template is defined"))
    (unless (member type *lapack-types*)
      (error "No recognized LAPACK type"))
    (destructuring-bind (return-value . argument-specifications)
	template
      (let* ((c-name (lapack-external-name name type))
             (name (lapack-lisp-name name type))
	     (source
              `(fl.port:simplified-def-function
                ((,c-name ,name)
                 ,(intern (substitute #\- #\_ (string-upcase name))
                          "FL.LAPACK"))
                ,(loop for arg-spec in argument-specifications
                       for spec = (extract-applicable-spec arg-spec type)
                       when spec collect
                       (map-tree #'convert-to-alien-type
                                 (check-spec spec)))
                :returning ,(convert-to-alien-type (translate-types return-value type)))))
	(let ((*package* (find-package "FL.LAPACK")))
          (dbg :lapack "Generating:~%~S~%" source)
	  (eval source))))))

(defgeneric lapack-convert (arg)
  (:documentation "Convert argument for use in a LAPACK routine.")
  (:method (x) (error "Don't know to convert arg"))
  (:method ((x number)) x)
  (:method ((x vector)) (fl.port:vector-sap x)))

(defun call-lapack (routine &rest args)
  "Call the LAPACK routine @arg{routine} with the arguments @arg{args}.
NIL-arguments are discarded, arrays and standard-matrices are converted to
the necessary alien representation."
  ;; Furthermore, the routine is called
  ;; with the IEEE FP modes appropriately set for Fortran, and GC is disallowed
  ;; such that it does not interfere by moving the Lisp arrays after their
  ;; address has been calculated.
  (dbg :lapack "Calling with:~%~{~A~%~}~%" (remove nil args))
  (apply #'foreign-call routine
	 (loop for arg in args
               when arg collect (lapack-convert arg))))

(defun call-lapack-with-error-test (&rest args)
  "Wrapper for @function{call-lapack} which tests if the routine worked
satisfactorily."
  (let ((info (nth-value 1 (apply #'call-lapack args))))
    (assert (zerop info))))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Template definitions
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defparameter *nrm2* '("snrm2" "dnrm2" "scnrm2" "dznrm2"))
(define-lapack-template *nrm2* :base-type
  (n :int :copy)
  (x (* :element-type) :in)
  (incx :int :copy))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(define-lapack-template "gemm" :void
  (transa (* :char) :in)
  (transb (* :char) :in)
  (m :int :copy)
  (n :int :copy)
  (k :int :copy)
  (alphar :element-type :copy)
  (a (* :element-type) :in)
  (lda :int :copy)
  (b (* :element-type) :in)
  (ldb :int :copy)
  (beta :element-type :copy)
  (c (* :element-type) :in)
  (ldc :int :copy))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(define-lapack-template "ggev" :void
  (jobvl (* :char) :in)
  (jobvr (* :char) :in)
  (n :int :copy)
  (a (* :element-type) :in)
  (lda :int :copy)
  (b (* :element-type) :in)
  (ldb :int :copy)
  (alphar (* :element-type) :in)
  (:select :number-type
	   (:real (alphai (* :element-type) :in)))
  (beta (* :element-type) :in)
  (vl (* :element-type) :in)
  (ldvl :int :copy)
  (vr (* :element-type) :in)
  (ldvr :int :copy)
  (work (* :element-type) :in)
  (lwork :int :copy)
  (:select :number-type
	   (:complex (rwork (* :base-type) :in)))
  (info :int :in-out))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(define-lapack-template "gesv" :void
  (n :int :copy)
  (nrhs :int :copy)
  (a (* :element-type) :in)
  (lda :int :copy)
  (ipiv (* :int) :in)
  (b (* :element-type) :in)
  (ldb :int :copy)
  (info :int :in-out))

(define-lapack-template "getrf" :void
  (m :int :copy)
  (n :int :copy)
  (a (* :element-type) :in)
  (lda :int :copy)
  (ipiv (* :int) :in)
  (info :int :in-out))

(define-lapack-template "getrs" :void
  (trans (* :char) :in)
  (n :int :copy)
  (nrhs :int :copy)
  (a (* :element-type) :in)
  (lda :int :copy)
  (ipiv (* :int) :in)
  (b (* :element-type) :in)
  (ldb :int :copy)
  (info :int :in-out))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defparameter *hegv*
  '("sspgv" "dspgv" "chegv" "zhegv"))

(define-lapack-template *hegv* :void
  (itype :int :copy)
  (jobz (* :char) :in)
  (uplo (* :char) :in)
  (n :int :copy)
  (a (* :element-type) :in)
  (:select :number-type (:complex (lda :int :copy)))
  (b (* :element-type) :in)
  (:select :number-type (:complex (ldb :int :copy)))
  (w (* :base-type) :in)
  (:select :number-type (:real (z (* :element-type) :in)))
  (:select :number-type (:real (ldz :int :copy)))
  (work (* :element-type) :in)
  (lwork :int :copy)
  (:select :number-type (:complex (rwork (* :base-type) :in)))
  (info :int :in-out))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Testing
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


#+(or)
(progn  ; geht nicht
  (cffi::defcfun ("dnrm2_" dnrm2) :double
    (len :int)
    (x (:pointer :double))
    (inc :int))
  (cffi:with-foreign-object (x :double 10)
    (dotimes (i 10)
      (setf (cffi:mem-aref x :double i) (coerce i 'double-float)))
    (dnrm2 10 x 1)))

(defun test-lapack ()
  
  (lapack-lisp-name "gemm" :float)
  (lapack *hegv* :double)
  (lapack-external-name '("sspgv" "dspgv" "chegv" "zhegv") :complex-float)
  
  (check-spec '(jobvl (* :char) :in))

  (let ((x (constant-vector 10 1.0d0)))
    (call-lapack (lapack *nrm2* :double) 10 x 1))
  (let ((x (constant-vector 10 1.0f0)))
    (call-lapack (lapack *nrm2* :float) 10 x 1))
  (let ((x (constant-vector 10 #C(1.0d0 2.0d0))))
    (call-lapack (lapack *nrm2* :complex-double) 10 x 1))
  
  (let ((A (coerce #(2.0 0.0 0.0 0.0 2.0 0.0 0.0 0.0 2.0)
		   '(simple-array double-float (*))))
	(ipiv (make-array 3 :element-type '(unsigned-byte 32) :initial-element 0))
	(b (constant-vector 3 1.0d0)))
    (assert (equal (multiple-value-list (call-lapack (lapack "gesv" :double) 3 1 A 3 ipiv b 3 0))
		   '(nil 0)))
    b)

  (let ((A (coerce #(2.0 1.0 1.0 0.0 2.0 1.0 0.0 0.0 2.0)
		   '(simple-array double-float (*))))
	(ipiv (make-array 3 :element-type '(unsigned-byte 32) :initial-element 0)))
    (assert (equal (multiple-value-list (call-lapack (lapack "getrf" :double) 3 3 A 3 ipiv 0))
		   '(nil 0)))
    A)

  (let ((A (make-array 4 :initial-element 1.0 :element-type 'double-float))
        (B (make-array 4 :initial-element 2.0 :element-type 'double-float))
        (C (make-array 4 :initial-element 0.0 :element-type 'double-float)))
    (call-lapack (lapack "gemm" :double) "N" "N" 2 2 2 1.0 A 2 B 2 1.0 C 2)
    C)
  )

(fl.tests:adjoin-test 'test-lapack)
;;; (test-lapack)
