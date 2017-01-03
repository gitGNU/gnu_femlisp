#-(and linux)
(error "CL-CPU-AFFINITY is Linux only.")

(defpackage :cl-cpu-affinity-system
  (:use :cl :asdf))

(in-package :cl-cpu-affinity-system)

(defclass c-so-source-file (c-source-file) ())

(defvar *gcc* "/usr/bin/gcc")
(defvar *gcc-options* '("-shared" "-fPIC"))

(defmethod output-files ((o compile-op) (c c-so-source-file))
  (list (make-pathname :type "so"
                       :defaults (component-pathname c))))

(defmethod perform ((o load-op) (c c-so-source-file))
  (destructuring-bind (so) (input-files o c)
    (funcall (read-from-string "cffi:load-foreign-library") so)))

(defmethod perform ((o compile-op) (c c-so-source-file))
  (let ((target (output-file o c))
        (source (component-pathname c)))
    (unless (zerop (run-shell-command
                    "~A ~A ~{~A ~}-o ~A"
                    *gcc*
                    source
                    *gcc-options*
                    target))
      (error 'operation-error :operation o :component c))))

(defsystem :cl-cpu-affinity
  :depends-on (:cffi)
  :components
  ((c-so-source-file "cpu-affinity-wrapper")
   (:file "package")
   (:file "cpu-affinity" :depends-on ("package" "cpu-affinity-wrapper")))
  )
