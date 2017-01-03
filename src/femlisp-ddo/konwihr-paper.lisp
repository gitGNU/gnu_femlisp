(in-package :fl.application)

#|

;;; Compiling

sbcl --eval "(asdf:req :femlisp-save-core)"

;;; Running (while commenting out different parts of konwihr-paper.lisp)

;;; On Laptop
likwid-perfctr -f -g MEM -C 1 femlisp --dynamic-space-size 8000 --load "konwihr-paper.lisp"

;;; On Sultana
likwid-perfctr -f -g MEM -C 1,5,9,13,17,21,25,29,33,37,41,45 femlisp --dynamic-space-size 240000 --load "konwihr-paper.lisp"

Alternative:
taskset -c 1,5,9,13,17,21,25,29,33,37,41,45  femlisp --dynamic-space-size 240000 --load "konwihr-paper.lisp"
|#

;;; <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
;;; Preliminaries:

(dbg-on :log-times)

;;; status output
(print (sb-cpu-affinity:get-cpu-affinity-mask))
(terpri)
(force-output)

;; (new-kernel)
;; (print (pwork (_ (princ-to-string (sb-cpu-affinity:get-cpu-affinity-mask)))))
;; (terpri)
;; (force-output)


#|
;;; >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

;;; <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
;;; Erste Rechnung auf Sultana

;;; umgekehrte Ordnung weil sonst die Affinitaet nicht mehr hochgesetzt werden kann
(elahom-performance-test '(12 11 10 9 8 7 6 5 4 3 2 1) :levels 3)

;;; >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
|#

#|
;;; <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
;;; Zweite Rechnung auf Sultana

(elahom-performance-test '(12) :levels 4)

;;; >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
|#


;;; <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
;;; MPI-parallele Rechnung auf Sultana
;;; sh sultana-worker.sh

(fl.application::elahom-performance-calculation :levels 2)

#|
(in-package :ddo-femlisp-test)

;;; (lfarm:end-kernel)
;;; (load "../connect-to-mpi-workers.lisp")
(ddo (princ-to-string (sb-cpu-affinity:get-cpu-affinity-mask)))
(ddo (dbg-on :communication)
     (dbg-on :distribute)
     (dbg-on :log-times))
(ddo (time (loop repeat 100000000 for x from 0.0d0 by 1.0 sum x)))
(ddo (fl.application::elahom-performance-calculation :levels 2))
|#

#| old stuff

(dbg-on :local-mat-pool)
(dbg-on :log-times)
(loop for nr-threads in '(1 1 2 3 4 5 6 7 8 9 10 11 12) do
  (sb-ext:gc :full t)
  (new-kernel nr-threads)
  (format t "~&*** ~D kernels ***~%" nr-threads)
  (time (pwork (lambda () (elahom-performance-calculation :levels 2 :output 1)))))

|#

#|
1 89
2 105
3 130
4 155
5 180
6 212
Unterbrochen mit:
 ** On entry to DGEMM  parameter number  1 had an illegal value
Zweite Rechnung -> Unterbrechung schon bei 5 workern

-> speedups (1.0f0 1.6952381f0 2.0538461f0 2.2967741f0 2.4722223f0 2.518868f0)
|#

#| ;;; interactive MPI
;;; (lfarm:end-kernel)
;;; (load "../connect-to-mpi-workers.lisp")
(ddo (princ-to-string (sb-cpu-affinity:get-cpu-affinity-mask)))
(ddo (dbg-on :communication)
     (dbg-on :distribute)
     (dbg-on :log-times))
(ddo (time (loop repeat 100000000 for x from 0.0d0 by 1.0 sum x)))
(ddo (fl.application::elahom-performance-calculation :levels 2))

|#
