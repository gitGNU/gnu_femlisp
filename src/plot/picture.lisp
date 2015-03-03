;;; -*- mode: lisp; fill-column: 64; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; picture.lisp - picture definition, input/output
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;
;;; Copyright (C) 2003 Nicolas Neuss, University of Heidelberg.
;;; Copyright (C) 2015 Nicolas Neuss, University of Erlangen-Nuremberg.
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

(defpackage "FL.PICTURE"
  (:use "COMMON-LISP" "FL.MACROS" "FL.UTILITIES" "FL.MATLISP")
  (:export "PICTURE" "PICTURE-WIDTH" "PICTURE-HEIGHT" "PIC-REF"
           "MAKE-PICTURE" "UNSIGNED-BYTE-PICTURE"
           "WRITE-PICTURE-TO-STREAM" "WRITE-PICTURE-TO-FILE")
  (:documentation "This package provides routines for handling pictures.  A
picture is implemented as a full-tensor of rank 2 or 3."))

(in-package :fl.picture)

(defclass picture () ()
  (:documentation "Mixin for pictures."))

(defun picture (type)
  "Construct a picture with entries of @arg{type}."
  (fl.amop:find-programmatic-class
   (list (find-class 'picture) (find-class 'full-tensor) (store-vector type))
   (intern (format nil "~A" (list 'PICTURE type)))))

(defun picture-width (pic) (aref (dimensions pic) 0))
(defun picture-height (pic) (aref (dimensions pic) 1))
(defun picture-length (pic) (the fixnum (* (picture-width pic) (picture-height pic))))

(defmethod pic-ref ((pic picture) i j)
  (aref (store pic) (+ i (* j (picture-width pic)))))

(defmethod (setf pic-ref) (value (pic picture) i j)
  (setf (aref (store pic) (+ i (* j (picture-width pic))))
	value))

;;; Some constructors

(defun make-picture (width height &key (element-type 'single-float))
  (make-instance (picture element-type) :dimensions (fixnum-vec width height)))

(defmethod make-analog ((picture picture))
  (make-instance (picture (element-type picture))
		 :dimensions (dimensions picture)))

(defun map-picture (new-element-type func pic)
  (make-instance (picture new-element-type)
		 :dimensions (dimensions pic)
		 :store
		 (map `(simple-array ,new-element-type (*))
		      func (store pic))))

(defgeneric write-picture-to-stream (stream pic &key &allow-other-keys))

(defmethod write-picture-to-stream (stream (pic picture) &key (type :jpg))
  (let ((pic (unsigned-byte-picture pic)))
    (let ((nx (picture-width pic))
          (ny (picture-height pic)))
      #+cl-gd
      (cl-gd:with-image* (nx ny)
        (let ((palette (map 'vector (_ (cl-gd:allocate-color _ _ _))
                            (range< 0 256))))
          ;; set background (white) and make it transparent
          (setf (cl-gd:transparent-color) (aref palette 0))
          (loop for x below nx do
               (loop for y below ny do
                    (cl-gd:set-pixel x y :color (aref palette (pic-ref pic x y)))))
          (cl-gd:write-image-to-stream stream type)))
      #-cl-gd
      (error "Please ensure that the cl-gd library is available."))))

(defvar *plot-picture-file* (fl.start:femlisp-pathname "images/plot-picture.pgm"))

(defgeneric write-picture-to-file (pic &key &allow-other-keys))

(defmethod write-picture-to-file (pic &key (filename *plot-picture-file*) (type :jpg))
  (with-open-file (stream filename :direction :output :if-exists :supersede
                          :external-format :latin-1)
    (write-picture-to-stream stream pic :type type)))

(defun unsigned-byte-picture (pic)
  (if (subtypep (element-type pic) '(unsigned-byte 8))
      pic
      (map-picture '(unsigned-byte 8)
		   #'(lambda (val)
		       (let ((int-val (truncate val)))
			 (cond ((minusp int-val) 0)
			       ((> int-val 255) 255)
			       (t int-val))))
		   pic)))

(defgeneric write-pgm (picture &key &allow-other-keys))

(defmethod write-pgm ((pic picture) &key (filename *plot-picture-file*) (format 'P5))
  "PGM output for pictures."
  (assert (= (rank pic) 2) () "Only 2D pics are implemented.")
  (with-open-file (stream filename :direction :output :if-exists :supersede
                          :element-type #+sbcl :default #-sbcl '(unsigned-byte 8))
    ;; problems with Allegro and SBCL
    (format stream "~A~%~D ~D~%255~%" format (picture-width pic) (picture-height pic))
    (let ((byte-pic (unsigned-byte-picture pic)))
      (ecase format
	(P2 (loop for value across (store byte-pic) do
		 (format stream "~D " value)))
	(P5 (loop for value across (store byte-pic) do
		 (write-byte value stream)))))))

(defun load-pgm (&key (filename *plot-picture-file*))
  "PGM input for pictures.  Does not work for binary streams
because mixing ASCII and binary in files is non-standard."
  (with-open-file (stream filename :direction :input)
    (let* ((format (read stream))
	   (width (read stream))
	   (height (read stream))
	   (range (read stream))
	   (pic (make-picture width height :element-type '(unsigned-byte 8)))
	   (values (store pic)))
      (assert (= range 255))
      (case format
	(P2
	 (dotimes (i (* width height))
	   (setf (aref values i) (read stream))))
	(P5 
	 (dotimes (i (* width height))
	   (setf (aref values i) (read-byte stream)))))
      ;; return picture
      pic)))

(defun display-picture-running-p ()
  (let ((process (fl.port:run-program
                  (probe-file "/bin/ps")
                  '("-fwwu" "neuss")
                  :wait nil :input :stream :output :stream :error-output nil)))
    (unwind-protect
         (with-open-stream (stream (fl.port:process-output process))
           (loop for line = (read-line stream nil) while line
              thereis (and (search "display" line)
                           (search (namestring (truename *plot-picture-file*)) line))))
      (fl.port:process-close process))))

;;; (display-picture-running-p)

(defun ensure-display-running ()
  (unless (display-picture-running-p)
    (fl.port:run-program
     "/usr/bin/display"
     (list "-update" "1" (namestring (truename *plot-picture-file*)))
     :directory (make-pathname :directory (pathname-directory *plot-picture-file*))
     :wait nil)))

(defmethod fl.plot:plot ((pic picture) &key &allow-other-keys)
  (ensure-display-running)
  (write-pgm pic :filename *plot-picture-file*))

(defun test-picture ()
  (let* ((N 100)
	 (pic (make-instance (picture '(unsigned-byte 8))
			     :dimensions (vector N N))))
    (dotimes (i N)
      (dotimes (j N)
	(setf (pic-ref pic i j)
	      (if (< i j)
		  0
		  (if (> i 50)
		      150
		      255)))))
    (fl.plot:plot pic))
  (display-picture-running-p)
  ;;; (ensure-display-running)
  )
