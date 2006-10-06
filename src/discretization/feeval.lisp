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

(in-package :fl.discretization)

(defmethod fe-local-value ((asv <ansatz-space-vector>) cell local-pos)
  "Evaluates a finite element function in CELL at LOCAL-POS."
  (let* ((fe (get-fe (fe-class asv) cell))
	 (f-values (get-local-from-global-vec cell fe asv))
	 (nr-comps (nr-of-components fe))
	 (components (components fe))
	 (result (make-array nr-comps :initial-element nil)))
    (dotimes (comp nr-comps result)
      (let ((shape-vals
	     (ensure-matlisp (map 'double-vec (rcurry #'evaluate local-pos)
				  (fe-basis (aref components comp)))
			     :row)))
	(setf (aref result comp)
	      (m* shape-vals (if (vectorp f-values)
				 (aref f-values comp)
				 f-values)))))))

(defmethod fe-value ((asv <ansatz-space-vector>) global-pos)
  "Evaluates a finite element function at GLOBAL-POS."
  (let ((cell (find-cell-from-position (mesh asv) global-pos)))
    (fe-local-value asv cell (global->local cell global-pos))))

(defmethod fe-local-gradient ((asv <ansatz-space-vector>) cell local-pos)
  "Evaluates a finite element gradient on cell at local coordinates local-pos."
  (let* ((fe (get-fe (fe-class asv) cell))
	 (f-values (get-local-from-global-vec cell fe asv))
	 (nr-comps (nr-of-components fe))
	 (components (components fe))
	 (result (make-array nr-comps :initial-element nil)))
    (dotimes (comp nr-comps result)
      (let ((Dphi^-1 (m/ (local->Dglobal cell local-pos)))
	    (f-values (if (vectorp f-values)
			  (aref f-values comp)
			  f-values)))
	(loop+ (i (shape (fe-basis (aref components comp)))) do
	   (m-incf (aref result comp)
		   (scal (vref f-values i)
			 (m* (ensure-matlisp
			      (coerce (evaluate-gradient shape local-pos) 'double-vec)
			      :row)
			     Dphi^-1))))))))

(defmethod fe-gradient ((asv <ansatz-space-vector>) pos)
  "Evaluates a finite element gradient at point pos."
  (let ((cell (find-cell-from-position (mesh asv) pos)))
    (fe-local-gradient asv cell (global->local cell pos))))

(defmethod cell-integrate (cell x &key initial-value (combiner #'m+!)
			   (key #'identity) coeff-func)
  "Integrates the ansatz-space vector @arg{x} on @arg{cell}.  If
@arg{coeff-fun} is set it should be a function which expects keyword
arguments @code{:solution} and @code{:global}."
  (dbg :feeval "Cell: ~A" cell)
  (let* ((fe (get-fe (fe-class x) cell))
	 (qrule (quadrature-rule fe))
	 (x-values (get-local-from-global-vec cell fe x))
	 (result initial-value))
    (loop
     for local across (integration-points qrule)
     and ip-weight across (integration-weights qrule)
     for shape-vals across (ip-values fe qrule) ; (n-basis x 1)-matrix
     for shape-vals-transposed =  (if (typep fe '<vector-fe>)
				      (vector-map #'transpose shape-vals)
				      (transpose shape-vals))
     for xip = (if (typep fe '<vector-fe>)
		   (map 'simple-vector #'m* shape-vals-transposed x-values)
		   (m* shape-vals-transposed x-values))
     for value = (if coeff-func
		     (evaluate coeff-func (list :solution xip :local local
						:global (local->global cell local)))
		     xip)
     for contribution = (scal (* ip-weight
				 (area-of-span (local->Dglobal cell local)))
			      (funcall key value))
     do (setq result
	      (if result
		  (funcall combiner contribution result)
		  contribution)))
    result))

(defmethod fe-integrate ((asv <ansatz-space-vector>) &key cells skeleton
			 initial-value (combiner #'m+!) (key #'identity))
  "Integrates a finite element function over the domain.  key is a
transformer function, as always (e.g. #'abs if you want the L1-norm)."
  (let ((result initial-value))
    (flet ((accumulate (cell)
	     (let ((contribution
		    (cell-integrate
		     cell asv :initial-value initial-value
		     :combiner combiner :key key)))
	       (setq result (if result
				(funcall combiner contribution result)
				contribution)))))
      (cond (cells (mapc #'accumulate cells))
	    (t (doskel (cell (or skeleton (mesh asv)) :dimension :highest :where :surface)
		 (accumulate cell))))
      result)))

(defun update-cell-extreme-values (cell x min/max)
  (let* ((fe (get-fe (fe-class x) cell))
	 (nr-comps (nr-of-components fe))
	 (multiplicity (multiplicity x))
	 (x-values (get-local-from-global-vec cell fe x)))
    (unless (vectorp x-values)
      (setf x-values (vector x-values)))
    (let ((minimum (or (car min/max) (make-array nr-comps :initial-element nil)))
	  (maximum (or (cdr min/max) (make-array nr-comps :initial-element nil))))
      (dotimes (k nr-comps)
	(let ((comp-values (aref x-values k)))
	  (unless (aref minimum k)
	    (setf (aref minimum k) (matrix-slice comp-values :nrows 1)))
	  (unless (aref maximum k)
	    (setf (aref maximum k) (matrix-slice comp-values :nrows 1)))
	  (dotimes (j multiplicity)
	    (symbol-macrolet ((min_kj (vref (aref minimum k) j))
			      (max_kj (vref (aref maximum k) j)))
	      (dotimes (i (nrows comp-values))
		(let ((entry (mref comp-values i j)))
		  (setf min_kj (min min_kj entry))
		  (setf max_kj (max max_kj entry))))))))
      (setf (car min/max) minimum (cdr min/max) maximum)
      min/max)))

(defmethod fe-extreme-values ((asv <ansatz-space-vector>) &key cells skeleton)
  "Computes the extreme values of a finite element function over the domain
or some region.  The result is a pair, the car being the minimum values and
the cdr the maximum values.  Each part is a vector of the size of the
number of components with entries being matrices of size 1 times
multiplicity."
  (let ((min/max (cons nil nil)))
    (cond (cells
	   (loop for cell in cells do
		(update-cell-extreme-values cell asv min/max)))
	  (t (doskel (cell (or skeleton (mesh asv)) :dimension :highest :where :surface)
	       (update-cell-extreme-values cell asv min/max))))
    min/max))

