;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; feeval.lisp - Evaluation of FE functions
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

(in-package :discretization)

(defmethod fe-value ((asv <ansatz-space-vector>) pos)
  "Evaluates a finite element function at point pos."
  (declare (optimize (debug 3)))
  (let* ((cell (find-cell-from-position (mesh asv) pos))
	 (local (global->local cell pos))
	 (fe (get-fe (fe-class asv) cell))
	 (f-values (get-local-from-global-vec cell fe asv))
	 (nr-comps (nr-of-components fe))
	 (components (components fe))
	 (result (make-array nr-comps :initial-element nil)))
    (dotimes (comp nr-comps result)
      (let ((shape-vals
	     (ensure-matlisp (map 'double-vec (rcurry #'evaluate local)
				  (fe-basis (aref components comp)))
			     :row)))
	(setf (aref result comp)
	      (m* shape-vals (if (vectorp f-values)
				 (aref f-values comp)
				 f-values)))))))

(defmethod fe-gradient ((asv <ansatz-space-vector>) pos)
  "Evaluates a finite element function at point pos."
  (declare (type <ansatz-space-vector> asv))
  (let* ((cell (find-cell-from-position (mesh asv) pos))
	 (local (global->local cell pos))
	 (fe (get-fe (fe-class asv) cell))
	 (f-values (get-local-from-global-vec cell fe asv))
	 (result nil))
    (loop with Dphi^-1 = (m/ (local->Dglobal cell local))
	  for shape in (fe-basis fe)
	  and i from 0 do
	  (m-incf result
		  (scal (vec-ref f-values i)
			(m* (ensure-matlisp (evaluate-gradient shape local) :row)
			    Dphi^-1))))
    result))

(defmethod cell-integrate (cell x &key (initial-value 0.0d0) (combiner #'+) (key #'identity))
  (let* ((fe-class (fe-class x))
	 (fe (get-fe fe-class cell))
	 (qrule (quadrature-rule fe-class fe))
	 (x-values (get-local-from-global-vec cell fe x))
	 (result initial-value))
    (loop for ip in (integration-points qrule)
	  for shape-vals in (ip-values fe qrule) ; (n-basis x 1)-matrix
	  for shape-vals-transposed =  (if (typep fe '<vector-fe>)
					   (vector-map #'transpose shape-vals)
					   (transpose shape-vals))
	  for xip = (if (typep fe '<vector-fe>)
			(map 'simple-vector #'m* shape-vals-transposed x-values)
			(m* shape-vals-transposed x-values))
	  do
	  (setq result
		(funcall combiner
			 (scal (* (ip-weight ip)
				  (area-of-span (local->Dglobal cell (ip-coords ip))))
			       (funcall key xip))
			 result)))
    result))

(defmethod fe-integrate ((asv <ansatz-space-vector>) &key cells skeleton
			 (initial-value 0.0d0) (combiner #'+) (key #'identity))
  "Integrates a finite element function over the domain.  key is a
transformer function, as always (e.g. #'abs if you want the L1-norm)."
  (let ((result initial-value))
    (flet ((accumulate (cell)
	     (setq result
		   (funcall combiner
			    (cell-integrate
			     cell asv :initial-value initial-value
			     :combiner combiner :key key)
			    result))))
      (cond (cells (mapc #'accumulate cells))
	    (t (doskel (cell (or skeleton (mesh asv)) :dimension :highest :where :surface)
		 (accumulate cell))))
      result)))

(defun update-cell-extreme-values (cell x minimum maximum initialized-p)
  (let* ((fe-class (fe-class x))
	 (fe (get-fe fe-class cell))
	 (x-values (get-local-from-global-vec cell fe x)))
    (dotimes (k (length x-values))
      (let ((comp-values (aref x-values k)))
	(dotimes (j (multiplicity x))
	  (symbol-macrolet ((min_kj (matrix-ref minimum k j))
			    (max_kj (matrix-ref maximum k j)))
	    (dotimes (i (nrows comp-values))
	      (let ((entry (matrix-ref comp-values i j)))
		(if initialized-p
		    (progn 
		      (setf min_kj (min min_kj entry))
		      (setf max_kj (max max_kj entry)))
		    (setf max_kj (setf min_kj entry)))))))))))

(defmethod fe-extreme-values ((asv <ansatz-space-vector>) &key cells skeleton)
  "Computes the extreme values of a finite element function over the domain
or some region."
  (let* ((fe-class (fe-class asv))
	 (nr-comps (nr-of-components fe-class))
	 (multiplicity (multiplicity asv))
	 (initialized-p nil)
	 (minimum (make-real-matrix nr-comps multiplicity))
	 (maximum (make-real-matrix nr-comps multiplicity)))
    (flet ((accumulate (cell)
	     (update-cell-extreme-values cell asv minimum maximum initialized-p)
	     (setq initialized-p t)))
      (cond (cells (mapc #'accumulate cells))
	    (t (doskel (cell (or skeleton (mesh asv)) :dimension :highest :where :surface)
		 (accumulate cell))))
      (list minimum maximum))))

