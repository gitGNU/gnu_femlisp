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

;;; Evaluation interface

(defgeneric fe-local-value (asv cell local-pos)
  (:documentation
   "Evaluates a FE ansatz-space-vector in @arg{cell} at @arg{local-pos}."))
(defgeneric fe-value (asv global-pos)
  (:documentation
   "Evaluates a FE ansatz-space-vector at @arg{global-pos}."))
(defgeneric fe-local-gradient (asv cell local-pos)
  (:documentation
   "Evaluates the gradient of a FE ansatz-space-vector on @arg{cell} at @arg{local-pos}."))
(defgeneric fe-gradient (asv global-pos)
  (:documentation
   "Evaluates the gradient of the FE ansatz-space-vector @arg{asv} at @arg{global-pos}."))
(defgeneric cell-integrate (cell object &key &allow-other-keys)
  (:documentation "Integrates @arg{object} on @arg{cell}."))
(defgeneric fe-integrate (asv
                          &key cells skeleton initial-value combiner key coeff-func)
  (:documentation
  "Integrates a finite element function over the domain.  key is a
transformer function, as always (e.g. #'abs if you want the L1-norm)."))
(defgeneric fe-extreme-values (asv &key cells skeleton component)
  (:documentation
  "Computes the extreme values of a finite element function over the domain
or some region.  The result is a pair, the car being the minimum values and
the cdr the maximum values.  Each part is a matrix of the format ncomps x
multiplicity."))

;;; Implementation

(defmethod fe-local-value ((asv <ansatz-space-vector>) cell local-pos)
  "Evaluates a finite element function in @arg{cell} at @arg{local-pos}."
  (vector-map #'m*-tn
              (ip-values-at-point (get-fe (ansatz-space asv) cell) local-pos)
              (get-local-from-global-vec cell asv)))

(defmethod fe-value ((asv <ansatz-space-vector>) global-pos)
  "Evaluates a finite element function at @arg{global-pos}."
  (let ((cell (loop for threshold in '(nil 1.0e-12 1.0e-10 1.0e-8)
                    for cell = (let ((fl.mesh::*inside-threshold* threshold))
                                 (find-cell-from-position (mesh asv) global-pos))
                    when cell do (return cell))))
    (unless cell
      (error "No suitable cell found which may be due to round-off errors."))
    (fe-local-value asv cell (global->local cell global-pos))))

(defmethod fe-local-gradient ((asv <ansatz-space-vector>) cell local-pos)
  "Evaluates a finite element gradient on @arg{cell} at @arg{local-pos}."
  (vector-map
   (lambda (local-gradients component-values)
     (m*-tn (m/ (local->Dglobal cell local-pos))
            (m*-tn local-gradients component-values)))
   (ip-gradients-at-point (get-fe (ansatz-space asv) cell) local-pos)
   (get-local-from-global-vec cell asv)))

(defmethod fe-gradient ((asv <ansatz-space-vector>) pos)
  "Evaluates a finite element gradient at point pos."
  (let ((cell (find-cell-from-position (mesh asv) pos)))
    (fe-local-gradient asv cell (global->local cell pos))))

(defmethod cell-integrate (cell (func function)
                           &key initial-value (combiner #'m+!)
                           (quadrature-order *quadrature-order*))
  "Integrates @arg{func} over @arg{cell}."
  (dbg :feeval "Cell: ~A" cell)
  (lret* ((qrule (gauss-rule (mapcar #'dimension (factor-simplices cell))
                             (ceiling (/ (+ quadrature-order 1) 2))))
          (result initial-value))
    (loop
       for local across (integration-points qrule)
       and ip-weight across (integration-weights qrule)
       for value = (evaluate func (local->global cell local))
       for contribution = (scal (* ip-weight
                                   (area-of-span (local->Dglobal cell local)))
                                value)
       do (setq result
                (if result
                    (funcall combiner contribution result)
                    contribution)))))

(defun cell-volume (cell)
  (cell-integrate cell (_ 1.0) :quadrature-order (or *quadrature-order* 10)))

(defmethod cell-integrate (cell (x <ansatz-space-vector>)
                           &key initial-value (combiner #'m+!)
			   (key #'identity) coeff-func)
  "Integrates the ansatz-space vector @arg{x} on @arg{cell}.  If
@arg{coeff-fun} is set it should be a function which expects keyword
arguments @code{:solution} and @code{:global}."
  (dbg :feeval "Cell: ~A" cell)
  (lret* ((fe (get-fe (ansatz-space x) cell))
          (qrule (quadrature-rule fe))
          (x-values (get-local-from-global-vec cell x))
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
		     (evaluate coeff-func (list :solution (list xip) :local local
						:global (local->global cell local)))
		     xip)
     for contribution = (scal (* ip-weight
				 (area-of-span (local->Dglobal cell local)))
			      (funcall key value))
     do (setq result
	      (if result
		  (funcall combiner contribution result)
		  contribution)))))

(defmethod fe-integrate ((asv <ansatz-space-vector>) &key cells skeleton
			 initial-value (combiner #'m+!) (key #'identity) coeff-func)
  "Integrates a finite element function over the domain.  key is a
transformer function, as always (e.g. #'abs if you want the L1-norm)."
  (lret ((result initial-value))
    (flet ((accumulate (cell)
	     (let ((contribution
		    (cell-integrate
		     cell asv :initial-value initial-value
		     :combiner combiner :key key :coeff-func coeff-func)))
	       (setq result (if result
				(funcall combiner contribution result)
				contribution)))))
      (cond (cells (mapc #'accumulate cells))
	    (t (doskel (cell (or skeleton (mesh asv)) :dimension :highest :where :surface)
		 (accumulate cell)))))))

(defun compute-energy (asv &key cells skeleton)
  (lret* ((as (ansatz-space asv))
          (result (zeros (multiplicity as))))
    (flet ((accumulate (cell)
	     (let ((local-vec (get-local-from-global-vec cell asv))
                   (local-mat (getf (assemble-cell cell as :matrix t) :local-mat)))
               (for-each-key
                (lambda (i j)
                  (m+! (m* (transpose (vref local-vec i))
                           (m* (mref local-mat i j)
                               (vref local-vec j)))
                       result))
                local-mat))))
      (cond (cells (mapc #'accumulate cells))
	    (t (doskel (cell (or skeleton (mesh asv)) :dimension :highest :where :surface)
		 (accumulate cell)))))))

(defun update-cell-extreme-values (cell x min/max component)
  "Computes the extreme values of @arg{x} on @arg{cell}.  Note that the
current implementation works correctly only for Lagrange finite elements,
and even for those only approximate extrema are obtained."
  (multiple-value-bind (pos length)
      (component-position
       (components-of-cell cell (hierarchical-mesh x) (problem x))
       component)
    (when pos
      (let ((x-values (get-local-from-global-vec cell x))
	    (multiplicity (multiplicity x)))
	(symbol-macrolet ((minimum (car min/max))
			  (maximum (cdr min/max)))
	  (unless minimum
	    (setq minimum (zeros length multiplicity))
	    (setq maximum (zeros length multiplicity))
	    (dotimes (k length)
	      (let ((comp-values (aref x-values (+ pos k))))
		(dotimes (j multiplicity)
		  (setf (mref minimum k j) (mref comp-values 0 j)
			(mref maximum k j) (mref comp-values 0 j))))))
	  (dotimes (k length)
	    (let ((comp-values (aref x-values (+ pos k))))
	      (dotimes (j multiplicity)
		(symbol-macrolet ((min_kj (mref minimum k j))
				  (max_kj (mref maximum k j)))
		    (dotimes (i (nrows comp-values))
		      (let ((entry (mref comp-values i j)))
			(setf min_kj (min min_kj entry))
			(setf max_kj (max max_kj entry))))))))))))
  min/max)

(defmethod fe-extreme-values ((asv <ansatz-space-vector>) &key cells skeleton component)
  "Computes the extreme values of a finite element function over the domain
or some region.  The result is a pair, the car being the minimum values and
the cdr the maximum values.  Each part is a matrix of the format ncomps x
multiplicity."
  (unless component
    (setq component (first (components (problem asv))))
    (when (listp component)
      (setq component (first component))))
  (let ((min/max (cons nil nil)))
    (cond (cells
	   (loop for cell in cells do
		(update-cell-extreme-values cell asv min/max component)))
	  (t (doskel (cell (or skeleton (mesh asv)) :where :surface)
	       (update-cell-extreme-values cell asv min/max component))))
    min/max))

