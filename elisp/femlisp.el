;;;  -*- Mode: Emacs-Lisp -*-

(provide 'femlisp)

(require 'slime)

(defun femlisp ()
  "Start femlisp by loading the file start.lisp from the directory in
*femlisp-root*."
  (interactive)
  (slime-start-and-load (concat *femlisp-root* "start.lisp")
			"FL.APPLICATION")
  )

;;; Improved DEFMETHOD indenting
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
(put 'with-properties 'lisp-indent-function
     (get 'destructuring-bind 'lisp-indent-function))
(put 'with-items 'lisp-indent-function
     (get 'destructuring-bind 'lisp-indent-function))

;;; fonts
(font-lock-add-keywords
 'lisp-mode
 '(("\\<\\(awhen\\|whereas\\|aif\\)\\>" . font-lock-keyword-face)))


;; Old stuff
;; (setq swank::*start-swank-in-background* t)

;; (add-hook 'lisp-mode-hook (lambda () (slime-mode t)))
;; (add-hook 'inferior-lisp-mode-hook (lambda () (inferior-slime-mode t)))

;; ;;; Hyperspec
(setq common-lisp-hyperspec-root "/usr/share/doc/hyperspec/")
;;(setq common-lisp-hyperspec-symbol-table nil)





