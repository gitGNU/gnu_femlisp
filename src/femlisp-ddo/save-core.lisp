(in-package :fl.application)

(sb-ext:save-lisp-and-die
 #+ddo "femlisp-ddo" #-ddo "femlisp"
 :executable t)
