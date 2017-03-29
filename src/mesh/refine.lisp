
;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; refine.lisp
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

(in-package :fl.mesh)

(file-documentation
 "This module provides definitions for refinement rules and functions for
refinements of skeletons.  The actual form of refinement information
depends on the cell type and is calculated separately in the files
vertex.lisp, simplex.lisp, and product-cell.lisp.")

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; refinement-rules
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass refinement-rule ()
  ((names :reader names :initarg :names
	  :documentation "Names identifying the rule.")
   (refcell :reader reference-cell :initarg :reference-cell :documentation
    "Reference cell for this refinement rule.")
   (boundary-refinement-rules :reader boundary-refinement-rules
    :documentation "Refinement rules for the sides required by this rule.")
   (refinement-info :reader refinement-info
                    :initarg :refinement-info :documentation
		     "Vector of refinement information for the children."))
  (:documentation "Rule for refining reference cells.  Those rules are
stored in the refine-info slot of the cell class."))

(defclass <anisotropic-rule-mixin> ()
  ((types :initarg :types))
  (:documentation "Mixin for an anisotropic refrule.  This type of rule
 includes copy and regular refinement as special cases."))

(defun anisotropic-rule-p (rule)
  (typep rule '<anisotropic-rule-mixin>))

(defun anisotropic-name (types)
  (intern (format nil "ANISOTROPIC-~{~:[F~;T~]~}" types)
          (find-package "KEYWORD")))

(defmethod print-object :after ((refrule refinement-rule) stream)
  (format stream "{~A}" (names refrule)))

(defun get-refinement-rule (cell id)
  "Finds the refinement rule for @arg{cell} defined by the @arg{id}.  This
@arg{id} can be a number (position of the rule, T (meaning 0), or some
symbol which is contained in the names of some rule.  Two values are
returned: the rule and its position in the refinement-rule vector."
  (let ((refrules (refinement-rules cell)))
    (whereas ((pos (cond ((and (numberp id) (< id (length refrules))) id)
			 ((eq id t) 0)
			 ((typep id 'refinement-rule) (position id refrules))
			 (t (position-if #'(lambda (rule) (member id (names rule)))
					 refrules)))))
      (values (aref refrules pos) pos))))

(defun (setf get-refinement-rule) (rule cell id)
  (with-cell-information (refinement-rules)
    cell
    (let ((pos (nth-value 1 (get-refinement-rule cell id))))
      (cond
	(pos
	 (unless (if (member :regular (names rule))
		     (zerop pos)
		     (plusp pos))
	     (error "Regular refinement should be at position zero."))
	 (setf (aref refinement-rules pos) rule))
	(t (when (and (numberp id) (not (= id (length refinement-rules))))
	     (error "The refinement-rule vector should be filled successively."))
	   (when (and (member :regular (names rule))
		      (plusp (length refinement-rules)))
	     (error "Regular refinement should be at position zero."))
	   (setf refinement-rules
		 (concatenate 'vector refinement-rules (vector rule)))))))
  rule)

(defun add-new-refinement-rule (rule)
  "Adds a new @arg{rule} to the refinement rules of @arg{cell}."
  (let ((refcell (reference-cell rule)))
    (assert (notany (_ (intersection (names rule) (names _)))
                    (refinement-rules refcell)))
    ;; ensure that the regular rule is already there
    (assert (get-refinement-rule refcell :regular))
    ;; add new rule 
    (setf (get-refinement-rule refcell (names rule))
          rule)))

(defun rule-position (id cell)
  (if (numberp id)
      id
      (let ((rules (refinement-rules cell)))
	(if (typep id 'refinement-rule)
	    (position id rules)
	    (position-if #'(lambda (rule)
			     (member id (names rule)))
			 rules)))))
      
(definline get-subcell-children (cell skel)
  "Returns a vector filled with all children of the subcells of
@arg{cell}."
  (map 'vector #'(lambda (cell) (children cell skel)) (subcells cell)))

(defun refine-cell-interior (refinfo cell subcell-refinements)
  "Refines @arg{cell} by filling the first position of the vector
@arg{subcell-refinements} with freshly generated children.  The other
positions should already be filled with the refinements of the cell's
boundary.  Returns a vector of children."
  (declare (type simple-vector subcell-refinements refinfo))
  ;; (declare (optimize speed))  ; ((speed (safety 1)))
  (lret ((my-refinement (make-array (length refinfo) :initial-element nil)))
    (setf (aref subcell-refinements 0) my-refinement)
    (loop for (refchild . vec) across refinfo ; ! was child-refcell
       and n of-type (integer 0 (#.array-dimension-limit)) from 0 do
       (symbol-macrolet ((new-child (aref my-refinement n)))
         (cond
           ((vertex-p refchild)
            (setf new-child (make-vertex (local->global cell vec))))
           (t (setf new-child (make-instance (class-of refchild)))
              (setf (slot-value new-child 'boundary)
                    (map-cell-vec #'(lambda (k&j)
                                      (aref (the simple-vector
                                                 (aref subcell-refinements (car k&j))) (cdr k&j)))
                                  vec))
              (when (mapped-p refchild)
                (setf (slot-value new-child 'mapping)
                      (compose-2 (cell-mapping cell) (cell-mapping refchild))))))))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; refinement information in a skeleton
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun refined-p (cell skeleton)
  (getf (skel-ref skeleton cell) 'CHILDREN))

(defun children (cell skeleton)
  (the (or null (simple-array <cell> (*)))
    (whereas ((refinement (getf (skel-ref skeleton cell) 'CHILDREN)))
      (if (consp refinement)
	  (cdr refinement)
	  refinement))))

(defun refinement (cell skeleton)
  "Returns the refinement of @arg{cell} in @arg{skeleton} as two values:
the rule and the children."
  (whereas ((refinement (getf (skel-ref skeleton cell) 'CHILDREN)))
    (if (consp refinement)
	(values (car refinement) (cdr refinement))
	(values (get-refinement-rule (reference-cell cell) 0) refinement))))

(defun refinement-rule (cell skel)
  "Returns the refinement rule of @arg{cell} in @arg{skel}."
  (nth-value 0 (refinement cell skel)))

(defun (setf children) (child-vec cell skeleton)
  (setf (getf (skel-ref skeleton cell) 'CHILDREN)
	child-vec))

(defun parent (cell skeleton)
  (the (or null <cell>)
    (getf (skel-ref skeleton cell) 'PARENT)))

(defun (setf parent) (parent cell skeleton)
  (setf (getf (skel-ref skeleton cell) 'PARENT)
	parent))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; generation of refinement information
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun skeleton->refinement-rule (names refcell skeleton)
  "Setup the refinement rule from a refined skeleton."
  (declare (optimize debug))
  (lret ((rule (make-instance 'refinement-rule :names names :reference-cell refcell)))
    (with-slots (boundary-refinement-rules refinement-info)
	rule
      ;; sanity checks
      (assert (reference-cell-p refcell))
      (let ((cells ()))
	(doskel (cell skeleton)
	  (assert (and (not (mapped-p cell))
		       (not (identified-p cell skeleton))))
	  (push cell cells))
	(assert (null (set-difference cells (coerce (subcells refcell) 'list)))))
      ;; internal setup of the refinement rule from the refined skeleton
      (setf boundary-refinement-rules
	    (map 'vector (rcurry #'refinement-rule skeleton) (boundary refcell)))
      (setf refinement-info
	    (let ((subcells (subcells refcell)))
	      (map 'vector
		   #'(lambda (child)
		       (cons (?2 (reference-cell child) child)
			     (if (vertex-p child)
				 (let* ((p1 (vertex-position child))
                                        (p2 (global->local refcell p1)))
				   (unless (mzerop (m- p1 p2) 1.0e-13)
				     (error "Should not happen"))
				   p2)
				 (map 'vector
				      #'(lambda (side)
					  (loop for k from 0 and subcell across subcells
					     for j = (position side (children subcell skeleton))
					     when j return (cons k j)))
				      (boundary child)))))
		   (children refcell skeleton)))))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; <child-info> - old format for refinement information
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defstruct (<child-info> (:conc-name child-))
  "This structure contains information about the children, their factor
simplices, and how to find their boundary sides.  The slots are:

class:  the child's cell-class

barycentric-corners:
   A list of the child's corners in barycentric coordinates

boundary-paths:
   For each side of the child this list contains a path in the form (i_1,
..., i_l, j): here i_1, ..., i_l are positions in subsequent boundary lists
starting from the boundary of the child's parent, and j is the position in
the refine-info vector of the boundary's parent.

transform-A, transform-b:  determine the transformation mapping for the child"
  
  (reference-cell (required-argument) :type <cell>)
  (barycentric-corners () :type list)
  (boundary-paths () :type list)
  (transform-A nil)
  (transform-b nil))

(deftype child-info-vec () '(simple-array <child-info> (*)))

(defun refine-info->refinement-rule (refcell refine-info)
  "This generates a refinement-rule from the older refine-info data."
  (let* ((types (make-list (length (factor-simplices refcell)) :initial-element t))
         (rule (make-instance
                (fl.amop:find-programmatic-class '(refinement-rule <anisotropic-rule-mixin>))
                :reference-cell refcell
                :names (if (vertex-p refcell)
                           '(:copy :regular :anisotropic-)
                           (list :regular (anisotropic-name types)))
                :types types)))
    (with-slots (boundary-refinement-rules refinement-info)
	rule
      (setf boundary-refinement-rules
	    (map 'vector
		 #'(lambda (side)
		     (aref (refinement-rules side) 0))
		 (boundary refcell)))
      (setf refinement-info
	    (let ((subcells (subcells refcell)))
	      (labels ((find-subcell (cell path)
			 (if (null path)
			     cell
			     (find-subcell (aref (boundary cell) (car path))
					   (cdr path)))))
		(vector-map
		 #'(lambda (child-info)
		     (let ((child-refcell (child-reference-cell child-info)))
		       (cons child-refcell
			     (if (vertex-p child-refcell)
				 (make-double-vec (dimension refcell) 0.5)
				 (map 'vector
				      #'(lambda (path)
					  (cons
					   (position (find-subcell refcell (butlast path))
						     subcells)
					   (car (last path))))
				      (child-boundary-paths child-info))))))
		 refine-info)))))
    rule))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; refinement-skeleton generation
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defparameter *refcell-refinement-memoize-depth* T
  "Depth up to which we store the refinements of the reference cell.  The
levels 0 and 1 are always memoized, becuause these are needed for
refinement operations.  NIL indicates no additional memoization beyond
levels 0 and 1, T indicates infinite depth.")

(defvar *refcell-refinement-table* (make-hash-table :test 'equal)
  "Table where refinements of reference cells are memoized.")

(defun induced-refinement-of-subcell-refcells (refcell rule subcell)
  "A refinement rule on refcell induces a refinement of all its subparts."
  (if (eql refcell subcell)
      rule
      (loop for side across (boundary refcell)
         for rule across (boundary-refinement-rules
                          (get-refinement-rule refcell rule))
           thereis (induced-refinement-of-subcell-refcells
                    side rule subcell))))

(defun refcell-refinement-skeleton (refcell &optional (level 1) (rule 0) reinit)
  "Returns an LEVEL times refined skeleton of REFCELL.  It is partially
memoized, see the documentation of *REFCELL-REFINEMENT-MEMOIZE-DEPTH*."
  (assert (reference-cell-p refcell))
  (setf rule (get-refinement-rule refcell rule))
  (assert rule)
  (when reinit
    ;; reinitialize the rule for refcell
    (loop for key being each hash-key of *refcell-refinement-table*
       when (and (eql (first key) refcell)
                 (eql (third key) rule))
       do (remhash key *refcell-refinement-table*)))
  (let ((key (list refcell level rule)))
    (or (gethash key *refcell-refinement-table*)
	(let ((result
	       (if (zerop level)
                   (skeleton refcell)
                   ;; the following makes sense at the moment only for special cases
                   ;; regular or copy refinement.
                   (refine (refcell-refinement-skeleton refcell (1- level) rule)
                           :indicator (constantly rule)
                           :decouple nil :highest t))))
	  (dbg :refine "Generating new refinement skeleton.~%Key=~A" key)
	  (when (or (<= level 1)	; is needed for refinement
		    (and *refcell-refinement-memoize-depth*
			 (or (eq *refcell-refinement-memoize-depth* T)
			     (<= level *refcell-refinement-memoize-depth*))))
	    (setf (gethash key *refcell-refinement-table*)
		  result))
	  result))))

(defun subcell-children (cell skeleton)
  "Returns a vector of all children of the subcells of @arg{cell} in
@arg{skeleton}."
  (coerce (loop for subcell across (subcells cell) nconcing
	       (loop for child across (children subcell skeleton)
		  collect child))
	  'vector))

(defun inner-refcell-children (refcell rule)
  "Returns the children of refcell."
  (children refcell (refcell-refinement-skeleton refcell 0 rule)))


(defun refcell-children (refcell rule)
  "Returns the children for refcell and subcells."
  (subcell-children refcell (refcell-refinement-skeleton refcell 0 rule)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; refine-info generation
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmethod generate-refine-info :before (refcell)
  (assert (reference-cell-p refcell)))

(defmethod generate-refine-info :after (refcell)
  "Generates the transformation mappings."
  ;; we have to ensure that the refined skeleton exists.  This refinement
  ;; of the standard cells needs the incomplete refine-info.
  ;; TODO: put in cell class information
  ;; generate regular refinement-rule from refine-info
  (declare (optimize debug))
  (with-cell-information (refinement-rules)
      refcell
    (setf refinement-rules
	  (vector (refine-info->refinement-rule
		   refcell (refine-info refcell)))))
  (assert (get-refinement-rule refcell 0))
  (assert (get-refinement-rule refcell t))
  (check (refcell-refinement-skeleton refcell 1 0))
  (let ((skel0 (refcell-refinement-skeleton refcell 0 0)))
    (doskel (cell skel0)
      (assert (children cell skel0))))
  (assert (inner-refcell-children refcell 0))
  (loop for child-info across (refine-info refcell)
	and child across (inner-refcell-children refcell 0) do
	(let ((corner (car (corners child))))
	  (setf (child-transform-b child-info) corner)
	  (assert corner)
	  (setf (child-transform-A child-info)
		(and (plusp (dimension child))
		     (l2Dg child (make-double-vec (dimension child))))))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Refinement of skeletons
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defgeneric refine-cell! (rule cell skel refined-skel refined-region)
  (:documentation "This local refinement function is used to assemble a
global refinement mapping which is represented by the skeleton
`refinement'.  It needs and ensures that the boundary of `cell' is
already refined.  An existing refinement of `cell' is simply kept."))

(defmethod refine-cell! (rule (cell <cell>) (skel <skeleton>)
			 (refined-skel <skeleton>) refined-region)
  (setf rule (get-refinement-rule cell rule))
  (with-slots (boundary-refinement-rules refinement-info) rule
    (unless (refined-p cell skel)
      ;; first ensure that the boundary is already refined
      (loop for side across (boundary cell)
         and side-rule across boundary-refinement-rules
	 unless (refined-p side skel) do
	 (refine-cell! side-rule side skel refined-skel refined-region))
      ;; then refine the interior
      (let* ((refine-info (refine-info cell))
	     (nr-of-children (length refine-info))
	     (children-vector
	      (refine-cell-interior
	       refinement-info cell (get-subcell-children cell skel))))
	;; set refinement information
	(setf (children cell skel)
	      (if (member :regular (names rule))
		  children-vector
		  (cons rule children-vector)))
	;; handle mappings
	(when (typep cell '<mapped-cell>)
	  (assert (member :regular (names rule)))
	  (dotimes (i nr-of-children)
	    (let ((child (aref children-vector i))
		  (child-info (aref refine-info i)))
	      (unless (vertex-p child)
		(change-class
		 child (mapped-cell-class (class-of child)
					  (typep cell '<distorted-cell>))
		 :mapping
		 (if (typep cell '<distorted-cell>) ; #+scl (scl-test cell)
		     (mapping cell)
		     (transform-function (mapping cell)
                                         :domain-transform
                                         (list (child-transform-A child-info)
                                               (child-transform-b child-info)))))))))
	  ;; put the pair cell/children-vector in the refined-refion
	(when refined-region
	  (setf (skel-ref refined-region cell) children-vector))
	;; finally, insert the children in the refined skeleton.
	(loop for child across children-vector do
	     (setf (parent child refined-skel) cell)))))
  nil)

(defmethod refine-cell! :after (rule (cell <cell>) (skel <skeleton>)
				(refined-skel <skeleton>) refined-region)
  "This after method ensures the refinement of identified cells and
identifies the children."
  (when (identified-p cell skel)
    (let ((identified-cells (identified-cells cell skel)))
      ;; ensure refinement of identified cells; this is a recursive call
      (when (loop with every-p = t
	       for id-cell in identified-cells
	       unless (refined-p id-cell skel) do
		 (setq every-p nil)
		 (refine-cell! rule id-cell skel refined-skel refined-region)
	       finally (return every-p))
	;; all identified cells are refined, now set identification for
	;; all children
	(loop for i from 0 below (length (children cell skel)) do
	   (identify (mapcar #'(lambda (id-cell) (aref (children id-cell skel) i))
			     identified-cells)
		     refined-skel))))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Skeleton refinement
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defgeneric do-refinement! (skel refined-skel task
                            &key refined-region)
  (:documentation "@arg{task} contains an association list of a cell/rule combination.
Note that the task is assumed to be consistent, which means that this function should
usually not be called directly by the user.  Its purpose is to allow for parallelization
to kick in, which needs to ensure consistent refinement of interfaces.")
  (:method ((skel <skeleton>) (refined-skel <skeleton>) (task list)
            &key refined-region)
    (loop for (cell . rule) in task do
      (unless (refined-p cell skel)
        (refine-cell! rule cell skel refined-skel refined-region)))))

(defgeneric update-refinement! (skel refined-skel &key region indicator refined highest)
  (:documentation "Low-level interface to refinement.  @arg{region} is
either NIL or a skeleton giving the region to be refined.  @arg{indicator}
is a function returning the refinement rule for each cell.  In the skeleton
@arg{refined} the actually refined cells can be collected, and
@arg{highest} determines if the indicator gives refinement for every cell
or only the highest-dimensional cells."))

(defmethod update-refinement! ((skel <skeleton>) (refined-skel <skeleton>)
			       &key region indicator refined highest)
  (assert (or region  ; h-mesh refinement
	      (not (eq skel refined-skel))))
  (ensure region skel)
  ;; check that refinement rules fit
  (let ((override-table (make-hash-table :test 'eq)))
    (doskel (cell region :direction :down)
      (whereas ((rule (aand (or (gethash cell override-table)
                                (funcall indicator cell))
                            (get-refinement-rule cell it))))
        (loop for side across (boundary cell)
           and side-rule across (boundary-refinement-rules rule)
           for indicated-side-rule =
             (get-refinement-rule side (funcall indicator side))
           do
             (unless (eq side-rule indicated-side-rule)
               (if indicated-side-rule
                   (error "Rule mismatch:~%cell ~A~%side ~A:~%~
                   induced:   ~A~%indicated: ~A"
                   cell side side-rule indicated-side-rule)
                   ;; we blame the inconsistency on the error indicator
                   ;; that marks only top-level cells
                   (setf (gethash side override-table) side-rule))
               ))))
    (let ((task ()))
      (doskel (cell region :direction :down :dimension (and highest :highest))
        (unless (refined-p cell skel)
          (whereas ((rule (or (gethash cell override-table)
                              (funcall indicator cell))))
            (let ((rule (get-refinement-rule cell rule)))
              (if rule
                  (push (cons cell rule) task)
                  (error "Refinement rule not found."))))))
      (do-refinement! skel refined-skel task :refined-region refined))))


(defgeneric refine (skel &key indicator &allow-other-keys)
  (:documentation "Refines @arg{skel} either locally or globally depending
on the @function{indicator}."))

(defparameter *constant-refinement-type*
  :regular
  "Standard refinement is regular refinement which has been implemented for
  all cell types in @femlisp{}.")

(defmethod refine ((skel <skeleton>)
                   &key (indicator (constantly *constant-refinement-type*))
                   highest (decouple t))
  "Refines the cells of @arg{skel} indicated by @arg{indicator}.  Returns
two values: the first is the refined skeleton, the second is the refinement
which is a skeleton for the old mesh referencing the refinement vectors.
This refinement algorithm usually makes sense only for global refinements.
Local refinements should be done with hierarchical-mesh structures."
  (let* ((refined-skel (make-analog skel))
	 (refinement (make-instance '<skeleton> :dimension (dimension skel))))
    (update-refinement! skel refined-skel :region skel :indicator indicator
			:refined refinement :highest highest)
    (when decouple
      (doskel (cell refined-skel)
	(remf (skel-ref refined-skel cell) :PARENT)))
    (values refined-skel refinement)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Tests
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(in-suite mesh-suite)

(test refine
  (is (= 4 (length (refinement-rules (n-cube 2)))))
  )
