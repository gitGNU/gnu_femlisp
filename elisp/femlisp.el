;;;  -*- Mode: Emacs-Lisp -*-

(provide 'femlisp)

(require 'slime)

(defun femlisp ()
  "Start femlisp by loading the file start.lisp from the directory in
*femlisp-root*."
  (interactive)
  (unless (slime-connected-p) (slime))
  (let ((startfile (concat *femlisp-root* "start.lisp")))
    (slime-display-output-buffer)
    (slime-interactive-eval (format "(load %S)" startfile))
    (slime-repl-set-package "FL.APPLICATION")))



