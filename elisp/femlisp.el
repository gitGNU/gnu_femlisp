;;;  -*- Mode: Emacs-Lisp -*-

(provide 'femlisp)

(require 'slime)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; SLIME enhancements
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun show-repl-maybe-start-slime ()
  (interactive)
  (if (slime-connected-p)
      (slime-display-output-buffer)
      (slime)))
(global-set-key [f9] 'show-repl-maybe-start-slime)

(defun slime-load-file-set-package (filename package)
  (let ((filename (slime-to-lisp-filename filename)))
    (slime-eval-async `(swank:load-file ,filename)
                      (lexical-let ((package package))
                        (lambda (ignored)
                          (slime-repl-set-package package))))))

(defun slime-start-and-load (filename &optional package)
  "Start Slime, if needed, load the current file and set the package."
  (interactive (list (expand-file-name (buffer-file-name))
                     (slime-find-buffer-package)))
  (cond ((slime-connected-p)
         (slime-load-file-set-package filename package))
        (t
         (slime-start-and-init (slime-lisp-options)
                               (slime-curry #'slime-start-and-load 
                                            filename package)))))

;; ;;; Hyperspec
;;(setq common-lisp-hyperspec-root "/usr/share/doc/hyperspec/")
;;(setq common-lisp-hyperspec-symbol-table nil)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Femlisp
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defvar *femlisp-root* nil
  "The location of the Femlisp directory.")

(defun femlisp ()
  "Start femlisp by loading the file start.lisp from the directory in
*femlisp-root*."
  (interactive)
  (slime-start-and-load (concat *femlisp-root* "start.lisp")
			"FL.APPLICATION")
  )

;;; Indenting for Femlisp's control constructs
(defun lisp-indent-defmethod (path state indent-point
                              sexp-column normal-indent)
  ;; Look for a method combination specifier...
  (let* ((combined (if (and (>= (car path) 3)
                            (null (cdr path)))
                       (save-excursion
                         (goto-char (car (cdr state)))
                         (forward-char)
                         (forward-sexp)
                         (forward-sexp)
                         (forward-sexp)
                         (backward-sexp)
                         (if (looking-at ":")
                             t
                             nil))
                       nil))
         (method (if combined
                     '(4 4 (&whole 4 &rest 1) &body)
                     '(4 (&whole 4 &rest 1) &body))))
    (funcall (function lisp-indent-259)
             method
             path state indent-point sexp-column normal-indent)))

(put 'defmethod 'common-lisp-indent-function 'lisp-indent-defmethod)

;;; Indentation and coloring for some Femlisp macros

(setq lisp-indent-function 'common-lisp-indent-function)

(put 'aif 'common-lisp-indent-function (get 'if 'common-lisp-indent-function))
(put 'awhen 'common-lisp-indent-function (get 'when 'common-lisp-indent-function))
(put 'whereas 'common-lisp-indent-function (get 'when 'common-lisp-indent-function))
(put 'dorows 'common-lisp-indent-function (get 'dotimes 'common-lisp-indent-function))
(put 'dorow 'common-lisp-indent-function (get 'dotimes 'common-lisp-indent-function))
(put 'doskel 'common-lisp-indent-function (get 'dotimes 'common-lisp-indent-function))
(put 'dotensor 'common-lisp-indent-function (get 'dotimes 'common-lisp-indent-function))
(put 'dovec 'common-lisp-indent-function (get 'dotimes 'common-lisp-indent-function))
(put 'dohash 'common-lisp-indent-function (get 'dotimes 'common-lisp-indent-function))
(put 'for 'common-lisp-indent-function 1)
(put 'for< 'common-lisp-indent-function 1)
(put 'multi-for 'common-lisp-indent-function 1)
(put 'with-properties 'common-lisp-indent-function
     (get 'destructuring-bind 'common-lisp-indent-function))
(put 'with-items 'common-lisp-indent-function
     (get 'destructuring-bind 'common-lisp-indent-function))

(put 'whereas 'lisp-indent-function (get 'when 'lisp-indent-function))
(put 'lret 'common-lisp-indent-function (get 'let 'common-lisp-indent-function))
(put 'lret* 'common-lisp-indent-function (get 'let* 'common-lisp-indent-function))
(put 'with-properties 'lisp-indent-function
     (get 'destructuring-bind 'lisp-indent-function))
(put 'with-items 'lisp-indent-function
     (get 'destructuring-bind 'lisp-indent-function))

;;; fonts
(font-lock-add-keywords
 'lisp-mode
 '(("\\<\\(awhen\\|aif\\|lret\\|lret[*]\\|whereas\\)\\>" . font-lock-keyword-face)))





